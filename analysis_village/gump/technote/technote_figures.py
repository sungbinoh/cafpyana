#!/usr/bin/env python3
"""Figure inventory / prune / install helper for an Overleaf (LaTeX) technote.

Stdlib only. Every subcommand takes --technote DIR (the Overleaf git checkout;
main file main.tex unless --main is given) and works on paths relative to it.

  list     every image referenced by a (non-commented) \\includegraphics in the
           compiled document tree, with \\foreach loops expanded; --table adds
           the .tex file, section, figure label and caption per reference
  unused   image files on disk that no compiled \\includegraphics references
           (--keep-commented also protects images named in commented-out lines
           or in .tex files that are never \\input, so they can be re-enabled)
  check    validate a manifest (TSV: <source rel. to --plotbase> <dest rel. to
           --technote>): every source exists, every dest is referenced, and
           list the referenced images the manifest does not cover
  install  copy every manifest source onto its dest (--dry-run to preview)

Typical:
  technote_figures.py list   --technote ../6973acd1e6f673f7a8b495e7 --table
  technote_figures.py unused --technote ... --keep-commented
  technote_figures.py check  --technote ... --plotbase ../plots-gumple-... --manifest figure_manifest.tsv
  technote_figures.py install --technote ... --plotbase ... --manifest figure_manifest.tsv
"""

import argparse
import itertools
import os
import re
import shutil
import sys

IMG_EXTS = (".pdf", ".png", ".jpg", ".jpeg", ".eps")

RE_INCLUDE = re.compile(r"\\(?:input|include)\s*\{([^}]*)\}")
RE_GRAPHICS = re.compile(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]*)\}")
RE_FOREACH = re.compile(r"\\foreach\s+((?:\\[A-Za-z@]+\s*/?\s*)+)\s+in\s*\{")
RE_SECTION = re.compile(r"\\(section|subsection|subsubsection)\*?\s*\{")
RE_LABEL = re.compile(r"\\label\s*\{([^}]*)\}")
RE_CAPTION = re.compile(r"\\caption\s*(?:\[[^\]]*\])?\s*\{")
RE_MACRO = re.compile(r"\\[A-Za-z@]+")


# ----------------------------------------------------------------------------
# tex reading
# ----------------------------------------------------------------------------
def strip_comments(text):
    """Blank everything from an unescaped % to end of line, keeping line
    structure (so character offsets stay meaningful)."""
    out = []
    for line in text.split("\n"):
        i = 0
        while True:
            j = line.find("%", i)
            if j < 0:
                break
            if j > 0 and line[j - 1] == "\\":
                i = j + 1
                continue
            line = line[:j] + " " * (len(line) - j)
            break
        out.append(line)
    return "\n".join(out)


def balanced(text, start):
    """text[start] == '{'; return index one past the matching '}'."""
    depth = 0
    i = start
    n = len(text)
    while i < n:
        c = text[i]
        if c == "\\":
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise ValueError("unbalanced braces from offset %d" % start)


def resolve_tex(root, name):
    p = os.path.join(root, name)
    if os.path.isfile(p):
        return p
    if os.path.isfile(p + ".tex"):
        return p + ".tex"
    return None


def tex_tree(root, main):
    """Ordered list of .tex files reachable from main via \\input/\\include
    (comments stripped before following)."""
    order = []
    seen = set()

    def visit(path):
        if path in seen:
            return
        seen.add(path)
        order.append(path)
        text = strip_comments(open(path, encoding="utf-8", errors="replace").read())
        for m in RE_INCLUDE.finditer(text):
            child = resolve_tex(root, m.group(1).strip())
            if child:
                visit(child)
            else:
                print("WARNING: cannot resolve \\input{%s} in %s"
                      % (m.group(1), os.path.relpath(path, root)), file=sys.stderr)

    visit(os.path.join(root, main))
    return order


# ----------------------------------------------------------------------------
# \foreach expansion
# ----------------------------------------------------------------------------
def foreach_loops(text):
    """[(body_start, body_end, [varnames], [value tuples])] for every
    \\foreach \\a[/\\b...] in {list}{body} in text."""
    loops = []
    for m in RE_FOREACH.finditer(text):
        names = [v.strip() for v in re.findall(r"\\[A-Za-z@]+", m.group(1))]
        list_start = m.end() - 1
        list_end = balanced(text, list_start)
        items = [s.strip() for s in text[list_start + 1:list_end - 1].split(",") if s.strip()]
        values = []
        for it in items:
            parts = [p.strip() for p in it.split("/")] if len(names) > 1 else [it]
            if len(parts) != len(names):
                parts = (parts + [""] * len(names))[:len(names)]
            values.append(tuple(parts))
        # body: the next '{' after the list
        k = list_end
        while k < len(text) and text[k] in " \t\r\n":
            k += 1
        if k >= len(text) or text[k] != "{":
            continue
        body_end = balanced(text, k)
        loops.append((k, body_end, names, values))
    return loops


def expand_macros(path, pos, loops):
    """All concrete paths for an \\includegraphics argument at offset pos,
    substituting enclosing \\foreach variables (cross product)."""
    enclosing = [l for l in loops if l[0] < pos < l[1]]
    if not RE_MACRO.search(path) or not enclosing:
        return [path]
    out = []
    for combo in itertools.product(*[l[3] for l in enclosing]):
        p = path
        subs = []
        for (_, _, names, _), vals in zip(enclosing, combo):
            subs += list(zip(names, vals))
        # longest macro names first so \Detector wins over \Det
        for name, val in sorted(subs, key=lambda nv: -len(nv[0])):
            p = re.sub(re.escape(name) + r"(?![A-Za-z@])", lambda _m: val, p)
        out.append(p)
    return out


def resolve_image(root, rel):
    rel = rel.strip()
    p = os.path.join(root, rel)
    if os.path.isfile(p):
        return rel
    base, ext = os.path.splitext(rel)
    if not ext:
        for e in IMG_EXTS:
            if os.path.isfile(os.path.join(root, rel + e)):
                return rel + e
    return None


# ----------------------------------------------------------------------------
# inventory
# ----------------------------------------------------------------------------
def context_for(text, pos):
    """(section title, figure label, caption) for an \\includegraphics at pos:
    the last \\section*/\\subsection* before it, and the \\label / \\caption
    inside the enclosing figure environment (searched forward to \\end{figure})."""
    sec = ""
    for m in RE_SECTION.finditer(text, 0, pos):
        end = balanced(text, m.end() - 1)
        sec = text[m.end():end - 1].strip()
    fig_end = text.find("\\end{figure", pos)
    if fig_end < 0:
        fig_end = min(len(text), pos + 2000)
    fig_start = text.rfind("\\begin{figure", 0, pos)
    if fig_start < 0:
        fig_start = max(0, pos - 2000)
    block = text[fig_start:fig_end]
    lab = RE_LABEL.search(block)
    cap = RE_CAPTION.search(block)
    caption = ""
    if cap:
        cend = balanced(block, cap.end() - 1)
        caption = " ".join(block[cap.end():cend - 1].split())
    return sec, (lab.group(1).strip() if lab else ""), caption


def inventory(root, main):
    """List of dicts: path (as written, expanded), resolved (on-disk rel path
    or None), tex, section, label, caption."""
    refs = []
    for tex in tex_tree(root, main):
        raw = open(tex, encoding="utf-8", errors="replace").read()
        text = strip_comments(raw)
        loops = foreach_loops(text)
        for m in RE_GRAPHICS.finditer(text):
            sec, label, caption = context_for(text, m.start())
            for p in expand_macros(m.group(1), m.start(), loops):
                refs.append(dict(path=p.strip(), resolved=resolve_image(root, p),
                                 tex=os.path.relpath(tex, root), section=sec,
                                 label=label, caption=caption))
    return refs


def all_tex_files(root):
    out = []
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d != ".git"]
        for f in fn:
            if f.endswith(".tex"):
                out.append(os.path.join(dp, f))
    return out


def protected_paths(root, main):
    """Images named by ANY \\includegraphics anywhere (commented lines and
    non-\\input .tex files included), for --keep-commented. Macro-bearing paths
    are expanded against the loops of the raw file."""
    used = set()
    for tex in all_tex_files(root):
        raw = open(tex, encoding="utf-8", errors="replace").read()
        # raw text: keep comments, but still need balanced-brace parsing
        try:
            loops = foreach_loops(raw)
        except ValueError:
            loops = []
        for m in RE_GRAPHICS.finditer(raw):
            for p in expand_macros(m.group(1), m.start(), loops):
                r = resolve_image(root, p)
                if r:
                    used.add(os.path.normpath(r))
    return used


def images_on_disk(root):
    out = []
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d != ".git"]
        for f in fn:
            if f.lower().endswith(IMG_EXTS):
                out.append(os.path.relpath(os.path.join(dp, f), root))
    return sorted(out)


def read_manifest(path):
    rows = []
    for ln, line in enumerate(open(path, encoding="utf-8"), 1):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split("\t")
        if len(parts) != 2:
            raise SystemExit("%s:%d: expected 2 tab-separated fields" % (path, ln))
        rows.append((parts[0].strip(), parts[1].strip()))
    return rows


# ----------------------------------------------------------------------------
# subcommands
# ----------------------------------------------------------------------------
def cmd_list(args):
    refs = inventory(args.technote, args.main)
    uniq = sorted({os.path.normpath(r["resolved"] or r["path"]) for r in refs})
    missing = sorted({r["path"] for r in refs if r["resolved"] is None})
    if args.table:
        print("\t".join(["image", "tex", "section", "label", "caption"]))
        for r in refs:
            print("\t".join([r["resolved"] or ("MISSING:" + r["path"]), r["tex"],
                             r["section"], r["label"], r["caption"]]))
    else:
        for p in uniq:
            print(p)
    print("# %d references, %d unique images, %d missing on disk"
          % (len(refs), len(uniq), len(missing)), file=sys.stderr)
    for p in missing:
        print("# MISSING: %s" % p, file=sys.stderr)
    return 1 if missing else 0


def cmd_unused(args):
    refs = inventory(args.technote, args.main)
    used = {os.path.normpath(r["resolved"]) for r in refs if r["resolved"]}
    prot = protected_paths(args.technote, args.main) if args.keep_commented else set()
    n_prot = 0
    total = 0
    for rel in images_on_disk(args.technote):
        n = os.path.normpath(rel)
        if n in used:
            continue
        if n in prot:
            n_prot += 1
            continue
        total += os.path.getsize(os.path.join(args.technote, rel))
        print(rel)
    print("# used: %d; protected (commented/orphan refs): %d; unused: %.1f MB"
          % (len(used), n_prot, total / 1e6), file=sys.stderr)
    return 0


def cmd_check(args):
    refs = inventory(args.technote, args.main)
    used = {os.path.normpath(r["resolved"] or r["path"]) for r in refs}
    rows = read_manifest(args.manifest)
    bad = 0
    dests = set()
    for src, dst in rows:
        sp = os.path.join(args.plotbase, src)
        if not os.path.isfile(sp):
            print("MISSING SOURCE: %s" % sp)
            bad += 1
        if os.path.normpath(dst) not in used:
            print("DEST NOT REFERENCED BY TECHNOTE: %s" % dst)
            bad += 1
        if dst in dests:
            print("DUPLICATE DEST: %s" % dst)
            bad += 1
        dests.add(os.path.normpath(dst))
    uncovered = sorted(u for u in used if u not in dests)
    print("# manifest rows: %d, problems: %d" % (len(rows), bad))
    print("# referenced images NOT covered by the manifest (%d) -- external / "
          "not regenerated:" % len(uncovered))
    for u in uncovered:
        print("  %s" % u)
    return 1 if bad else 0


def cmd_install(args):
    rows = read_manifest(args.manifest)
    n = 0
    for src, dst in rows:
        sp = os.path.join(args.plotbase, src)
        dp = os.path.join(args.technote, dst)
        if not os.path.isfile(sp):
            print("SKIP (missing source): %s" % sp)
            continue
        if args.dry_run:
            print("%s -> %s" % (sp, dp))
        else:
            os.makedirs(os.path.dirname(dp), exist_ok=True)
            shutil.copy2(sp, dp)
        n += 1
    print("# %s %d files" % ("would install" if args.dry_run else "installed", n))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--technote", required=True, help="Overleaf checkout dir")
        p.add_argument("--main", default="main.tex", help="main file (default main.tex)")

    p = sub.add_parser("list"); common(p)
    p.add_argument("--table", action="store_true",
                   help="TSV with tex file / section / label / caption per reference")
    p.set_defaults(fn=cmd_list)

    p = sub.add_parser("unused"); common(p)
    p.add_argument("--keep-commented", action="store_true",
                   help="also protect images named in commented-out lines or in "
                        ".tex files that are never \\input")
    p.set_defaults(fn=cmd_unused)

    p = sub.add_parser("check"); common(p)
    p.add_argument("--plotbase", required=True)
    p.add_argument("--manifest", required=True)
    p.set_defaults(fn=cmd_check)

    p = sub.add_parser("install"); common(p)
    p.add_argument("--plotbase", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(fn=cmd_install)

    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
