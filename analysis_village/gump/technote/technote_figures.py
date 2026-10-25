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

Paths are expanded through \\foreach loops (brace-aware: `a/{text, with, commas}`
items stay whole), zero-argument macros (\\newcommand / \\providecommand / \\def,
e.g. a figure-directory macro, also in \\input paths) and parameterized
\\newcommand macros whose body holds an \\includegraphics (expanded at each call
site; a macro that is never called contributes nothing, like a commented line).

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


def split_top(s, sep):
    """Split s at sep only where the brace depth is 0."""
    parts, depth, cur = [], 0, []
    i = 0
    while i < len(s):
        c = s[i]
        if c == "\\" and i + 1 < len(s):
            cur.append(s[i:i + 2])
            i += 2
            continue
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
        if c == sep and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    parts.append("".join(cur))
    return parts


def strip_outer(s):
    """'{abc}' -> 'abc' when the braces enclose the whole string."""
    s = s.strip()
    if s.startswith("{"):
        try:
            if balanced(s, 0) == len(s):
                return s[1:-1].strip()
        except ValueError:
            pass
    return s


# ----------------------------------------------------------------------------
# macros
# ----------------------------------------------------------------------------
RE_NEWCMD = re.compile(r"\\(newcommand|renewcommand|providecommand)\*?\s*"
                       r"(?:\{\s*(\\[A-Za-z@]+)\s*\}|(\\[A-Za-z@]+))\s*"
                       r"(?:\[(\d)\])?\s*(?:\[[^\]]*\])?\s*\{")
RE_DEF = re.compile(r"\\def\s*(\\[A-Za-z@]+)\s*\{")


class Macros:
    """Zero-argument macros (name -> replacement) and parameterized macros
    (name -> (nargs, body, defining file, body offset))."""

    def __init__(self):
        self.simple = {}
        self.param = {}

    def collect(self, text, path):
        """Record the definitions in text; return [(start, end)] spans of the
        parameterized definitions (their bodies are expanded at call sites)."""
        spans = []
        for m in RE_NEWCMD.finditer(text):
            kind, name = m.group(1), m.group(2) or m.group(3)
            try:
                end = balanced(text, m.end() - 1)
            except ValueError:
                continue
            body = text[m.end():end - 1]
            if m.group(4) and int(m.group(4)) > 0:
                self.param[name] = (int(m.group(4)), body, path, m.end())
                spans.append((m.start(), end))
            elif kind != "providecommand" or name not in self.simple:
                self.simple[name] = body.strip()
        for m in RE_DEF.finditer(text):
            try:
                end = balanced(text, m.end() - 1)
            except ValueError:
                continue
            self.simple[m.group(1)] = text[m.end():end - 1].strip()
        return spans

    def expand(self, s):
        """Substitute zero-argument macros (repeatedly, longest name first)."""
        for _ in range(5):
            new = s
            for name in sorted(self.simple, key=len, reverse=True):
                new = re.sub(re.escape(name) + r"(?![A-Za-z@])",
                             lambda _m, v=self.simple[name]: v, new)
            if new == s:
                break
            s = new
        return s

    def calls(self, text):
        """[(name, pos, [args])] for every call of a parameterized macro."""
        out = []
        if not self.param:
            return out
        pat = re.compile(r"(%s)(?![A-Za-z@])"
                         % "|".join(re.escape(n) for n in sorted(self.param, key=len, reverse=True)))
        for m in pat.finditer(text):
            nargs = self.param[m.group(1)][0]
            k, args = m.end(), []
            while len(args) < nargs:
                while k < len(text) and text[k] in " \t\r\n":
                    k += 1
                if k < len(text) and text[k] == "[":  # optional argument: skip
                    k = text.find("]", k) + 1
                    continue
                if k >= len(text) or text[k] != "{":
                    break
                try:
                    end = balanced(text, k)
                except ValueError:
                    break
                args.append(text[k + 1:end - 1])
                k = end
            if len(args) == nargs:
                out.append((m.group(1), m.start(), args))
        return out

    def instantiate(self, name, args):
        body = self.param[name][1]
        return re.sub(r"#(\d)", lambda m: args[int(m.group(1)) - 1]
                      if 0 < int(m.group(1)) <= len(args) else m.group(0), body)


def in_spans(pos, spans):
    return any(a <= pos < b for a, b in spans)


def resolve_tex(root, name):
    p = os.path.join(root, name)
    if os.path.isfile(p):
        return p
    if os.path.isfile(p + ".tex"):
        return p + ".tex"
    return None


def tex_tree(root, main, macros=None):
    """Ordered list of .tex files reachable from main via \\input/\\include
    (comments stripped before following). Macro definitions met on the way are
    recorded in macros (a Macros), so \\input{\\FigDir/x} resolves."""
    order = []
    seen = set()
    macros = Macros() if macros is None else macros

    def visit(path):
        if path in seen:
            return
        seen.add(path)
        order.append(path)
        text = strip_comments(open(path, encoding="utf-8", errors="replace").read())
        macros.collect(text, path)
        for m in RE_INCLUDE.finditer(text):
            child = resolve_tex(root, macros.expand(m.group(1).strip()))
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
        items = [s.strip() for s in split_top(text[list_start + 1:list_end - 1], ",") if s.strip()]
        values = []
        for it in items:
            parts = [strip_outer(p) for p in split_top(it, "/")] if len(names) > 1 else [strip_outer(it)]
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


def graphics_in(text, macros, skip_spans=()):
    """[(pos, path, ctx_text, ctx_pos)] for every \\includegraphics in text,
    \\foreach and macro expanded. Graphics inside skip_spans (parameterized
    macro definitions) are skipped; they are produced at each call site instead,
    with the instantiated body as the caption/label context."""
    out = []
    loops = foreach_loops(text)
    for m in RE_GRAPHICS.finditer(text):
        if in_spans(m.start(), skip_spans):
            continue
        for p in expand_macros(m.group(1), m.start(), loops):
            out.append((m.start(), macros.expand(p).strip(), text, m.start()))
    for name, pos, args in macros.calls(text):
        if in_spans(pos, skip_spans):
            continue
        body = macros.instantiate(name, args)
        body_loops = foreach_loops(body)
        for m in RE_GRAPHICS.finditer(body):
            for p in expand_macros(m.group(1), m.start(), body_loops):
                for q in expand_macros(p, pos, loops):  # call inside a \foreach
                    out.append((pos, macros.expand(q).strip(), body, m.start()))
    return out


def inventory(root, main):
    """List of dicts: path (as written, expanded), resolved (on-disk rel path
    or None), tex, section, label, caption."""
    refs = []
    macros = Macros()
    texs = tex_tree(root, main, macros)
    for tex in texs:
        text = strip_comments(open(tex, encoding="utf-8", errors="replace").read())
        spans = Macros().collect(text, tex)  # parameterized definitions in this file
        for pos, p, ctext, cpos in graphics_in(text, macros, spans):
            sec = context_for(text, pos)[0]
            _, label, caption = context_for(ctext, cpos)
            refs.append(dict(path=p, resolved=resolve_image(root, p),
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
    macros = Macros()
    texs = all_tex_files(root)
    raws = {}
    for tex in texs:  # all definitions first: orphan files may call each other's macros
        raws[tex] = open(tex, encoding="utf-8", errors="replace").read()
        try:
            macros.collect(raws[tex], tex)
        except ValueError:
            pass
    for tex in texs:
        raw = raws[tex]
        # raw text: keep comments, but still need balanced-brace parsing.
        # Graphics inside macro definitions count too (an uncalled macro is
        # like a commented-out line).
        try:
            refs = graphics_in(raw, macros)
        except ValueError:
            refs = [(m.start(), macros.expand(m.group(1)), raw, m.start())
                    for m in RE_GRAPHICS.finditer(raw)]
        for _, p, _, _ in refs:
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
