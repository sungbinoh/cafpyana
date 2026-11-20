#!/usr/bin/env python3
"""
match_fom.py -- rewrite cafpyana dataframes with a corrected BNB FOM.

Reads one or more `.df` files and writes, for each, a new `.df` that is a copy
of every key except the spill-info keys `bnb_0`, `bnb_1`, ... which gain the
recomputed FOM. The input is never modified.

The corrected values come from the beam-quality ROOT files written by
`bnb_fom_caf.py` (`beamqual_*_fom.root`). Those carry, per spill,
    FOM        the production value reproduced bug-for-bug  (used as a check)
    FOM_best   the corrected value: NaN where no FOM can be formed, else in
               [0, 1]; no -999, no 100+x "nominal width" encoding

Columns added to each bnb_* frame (the production FOM is left untouched):
    FOM_corrected           the corrected FOM
    FOM_corrected_matched   bool, spill found in the ROOT files
    fom_source              width used          (see fom_sources in the ROOT file)
    fom_status              status bits         (see fom_status_bits)

    # what is in here, and is the schema understood?
    ./match_fom.py mydata.df --root beamqual_run2_fom.root --inspect

    # report the match without writing anything
    ./match_fom.py mydata.df --root beamqual_run2_fom.root --dry-run

    # mydata.df -> mydata_fomfix.df
    ./match_fom.py mydata.df --root beamqual_run2_fom.root

    # many inputs into one directory, all three beam-quality files
    ./match_fom.py run*/*.df --root beamqual_*_fom.root --outdir fixed/

Join key
    (spill_time_sec, spill_time_nsec) as one int64 nanosecond stamp -- the same
    integers sbncode wrote into the spill record, so the match is exact and
    there is no tolerance to tune. In the run2 file the stamp is unique for
    5,345,417 of 5,345,446 spills; the 29 collisions are duplicate spill
    records, resolved by preferring the row with a finite FOM_best, then first
    seen. run/subrun are cross-checked when present, since a disagreement there
    means the timestamp joined to the wrong spill.

Verification
    The ROOT `FOM` column reproduces the production value bug-for-bug, so it
    must equal the FOM already in the dataframe. The report prints that
    agreement -- below ~100% the join is wrong, not the correction large.

Needs pandas + pytables + h5py (.df I/O) and uproot (.root).
"""
import argparse
import os
import re
import sys

import numpy as np

# cafpyana / BNBSpillInfo spellings, most specific first
SEC_CANDIDATES = ["spill_time_sec", "spill_time_s", "spill_time_secs", "time_sec", "sec"]
NSEC_CANDIDATES = ["spill_time_nsec", "spill_time_ns", "spill_time_nsecs", "time_nsec", "nsec"]
FLOAT_TIME_CANDIDATES = ["spill_time", "spill_time_f", "time"]
FOM_CANDIDATES = ["FOM", "fom"]

BNB_KEY_RE = re.compile(r"^/?bnb_(\d+)$")
NS = 1_000_000_000


def die(msg):
    print(f"match_fom: {msg}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------- column helpers
def flat_columns(df):
    """{lowercased name: column key}, tolerating cafpyana's MultiIndex columns.
    Indexes the first and last non-empty level and their dotted join, so
    'spill_time_sec', ('spill_time_sec', '') and ('bnbinfo', 'spill_time_sec')
    are all reachable."""
    out = {}
    for c in df.columns:
        if isinstance(c, tuple):
            parts = [str(x) for x in c if str(x) != ""]
            if not parts:
                continue
            for name in {parts[0], parts[-1], ".".join(parts)}:
                out.setdefault(name.lower(), c)
        else:
            out.setdefault(str(c).lower(), c)
    return out


def pick(cols, candidates):
    for name in candidates:
        if name.lower() in cols:
            return cols[name.lower()]
    return None


def ns_stamp(sec, nsec):
    return np.asarray(sec).astype(np.int64) * NS + np.asarray(nsec).astype(np.int64)


# ---------------------------------------------------------------- ROOT side
def load_lookup(paths, tree, fom_col, extra_cols):
    """Spill lookup from the beam-quality ROOT files.
    Returns (sorted unique int64 stamp, {column: array aligned with it}, meta)."""
    try:
        import uproot
    except ImportError:
        die("uproot is required to read the ROOT files (pip install uproot)")

    want = list(dict.fromkeys(
        ["run", "subrun", "spill_time_sec", "spill_time_nsec", "FOM", fom_col]
        + [c for c in extra_cols if c]))
    chunks, meta = [], {"n_raw": 0}
    for p in paths:
        f = uproot.open(p)
        trees = [k.split(";")[0] for k in f.keys()]
        if tree not in trees:
            die(f"{p}: no tree '{tree}' (has {trees})")
        t = f[tree]
        missing = [c for c in want if c not in set(t.keys())]
        if missing:
            die(f"{p}: tree '{tree}' lacks {missing}")
        chunks.append(t.arrays(want, library="np"))
        meta["n_raw"] += t.num_entries
        print(f"  {t.num_entries:>10,} spills  {os.path.basename(p)}:{tree}")

    cat = {k: np.concatenate([c[k] for c in chunks]) for k in want}
    stamp = ns_stamp(cat["spill_time_sec"], cat["spill_time_nsec"])
    # duplicate stamps: prefer a finite corrected FOM, then first seen
    order = np.lexsort((~np.isfinite(cat[fom_col]), stamp))
    srt = stamp[order]
    keep = np.concatenate(([True], srt[1:] != srt[:-1]))
    meta["n_dup"] = int((~keep).sum())
    sel = order[keep]
    meta["n_unique"] = int(keep.sum())
    return stamp[sel], {k: v[sel] for k, v in cat.items()}, meta


# ---------------------------------------------------------------- dataframe side
def bnb_keys(store):
    ks = [(int(BNB_KEY_RE.match(k).group(1)), k) for k in store.keys() if BNB_KEY_RE.match(k)]
    return [k for _, k in sorted(ks)]


def frame_stamp(df, cols, tol):
    sec, nsec = pick(cols, SEC_CANDIDATES), pick(cols, NSEC_CANDIDATES)
    if sec is not None and nsec is not None:
        return ns_stamp(df[sec], df[nsec]), f"{sec} + {nsec} (exact)"
    ft = pick(cols, FLOAT_TIME_CANDIDATES)
    if ft is not None:
        if tol <= 0:
            die(f"only a float time column ('{ft}') is present; rerun with --tolerance "
                f"(e.g. 1e-6) to allow nearest-time matching")
        return np.rint(np.asarray(df[ft], dtype=np.float64) * NS).astype(np.int64), f"{ft} (float)"
    die(f"no spill-time column found. Looked for {SEC_CANDIDATES + FLOAT_TIME_CANDIDATES}; "
        f"frame has {sorted(cols)[:40]}")


def match(stamp_df, stamp_lut, tol):
    """Lookup index per df row, -1 where unmatched."""
    if len(stamp_lut) == 0:
        return np.full(len(stamp_df), -1)
    idx = np.searchsorted(stamp_lut, stamp_df)
    if tol <= 0:
        ok = idx < len(stamp_lut)
        safe = np.where(ok, idx, 0)
        return np.where(ok & (stamp_lut[safe] == stamp_df), safe, -1)
    tol_ns = int(round(tol * NS))
    lo, hi = np.clip(idx - 1, 0, len(stamp_lut) - 1), np.clip(idx, 0, len(stamp_lut) - 1)
    dlo, dhi = np.abs(stamp_lut[lo] - stamp_df), np.abs(stamp_lut[hi] - stamp_df)
    best, dist = np.where(dlo <= dhi, lo, hi), np.minimum(dlo, dhi)
    return np.where(dist <= tol_ns, best, -1)


def report(key, df, stamp_df, how, idx, values, cols, fom_col, new_col):
    n, m = len(df), idx >= 0
    print(f"\n  {key}: {n:,} spill rows, key = {how}")
    print(f"    matched                {m.sum():,} / {n:,}  ({m.mean() if n else 0:.4%})")
    if n and not m.all():
        for i in np.where(~m)[0][:3]:
            print(f"      unmatched row {i}: {stamp_df[i] // NS}.{stamp_df[i] % NS:09d}")
    if not m.any():
        return

    for name in ("run", "subrun"):
        c = pick(cols, [name])
        if c is None:
            continue
        agree = np.asarray(df[c])[m].astype(np.int64) == values[name][idx[m]].astype(np.int64)
        print(f"    {name} agreement        {agree.mean():.4%}"
              f"{'' if agree.all() else '   <-- WRONG JOIN'}")

    fc = pick(cols, FOM_CANDIDATES)
    lhs = np.asarray(df[fc], dtype=np.float64)[m] if fc is not None else None
    if lhs is not None:
        ref = values["FOM"][idx[m]].astype(np.float64)
        both = np.isfinite(lhs) & np.isfinite(ref)
        if both.any():
            close = np.isclose(lhs[both], ref[both], atol=1e-4, rtol=0)
            print(f"    production FOM reproduced  {close.mean():.4%}"
                  f"{'' if close.mean() > 0.999 else '   <-- join or file mismatch'}")
    else:
        print("    (no FOM column in the frame; join verification skipped)")

    v = values[fom_col][idx[m]].astype(np.float64)
    print(f"    {new_col}: finite {np.isfinite(v).mean():.2%}, "
          f"in [0,1] {((v >= 0) & (v <= 1)).mean():.2%}")
    if lhs is not None:
        # the production FOM uses sentinels (-999 = none, 100+x = nominal width),
        # so a bare max|delta| is meaningless; split the categories out
        real = np.isfinite(lhs) & (lhs >= 0) & (lhs <= 1) & np.isfinite(v)
        if real.any():
            print(f"    both real FOMs ({real.sum():,} rows): changed "
                  f"{(~np.isclose(v[real], lhs[real], atol=1e-6)).mean():.2%}, "
                  f"max |delta| {np.nanmax(np.abs(v[real] - lhs[real])):.4g}")
        was999, was100 = np.isclose(lhs, -999.0, atol=1e-3), lhs > 100
        print(f"    was -999, now a FOM    {(was999 & np.isfinite(v)).sum():,}"
              f"   (still none {(was999 & ~np.isfinite(v)).sum():,})")
        print(f"    was 100+ (nominal)     {was100.sum():,}"
              f"   (now a real FOM {(was100 & np.isfinite(v)).sum():,})")
        print(f"    had a FOM, now NaN     "
              f"{(np.isfinite(lhs) & (lhs >= 0) & (lhs <= 1) & ~np.isfinite(v)).sum():,}")


def augment(df, idx, values, fom_col, new_col, extra_cols):
    out = df.copy()
    m = idx >= 0
    safe = np.where(m, idx, 0)
    nlev = getattr(out.columns, "nlevels", 1)
    cname = (lambda n: (n,) + ("",) * (nlev - 1)) if nlev > 1 else (lambda n: n)

    def take(name, fill=np.nan):
        return np.where(m, np.asarray(values[name])[safe].astype(np.float64), fill)

    out[cname(new_col)] = take(fom_col)
    out[cname(new_col + "_matched")] = m
    for c in extra_cols:
        out[cname(c)] = take(c, fill=-1)
    return out


# ---------------------------------------------------------------- output
def write_copy(src, dst, frames, formats):
    """dst = copy of every key in src except bnb_*, which come from `frames`.

    Non-bnb keys are copied at the HDF5 level (h5py), so memory stays flat no
    matter how large the event frames are and the stored layout is preserved
    exactly; only the bnb_* frames go through pandas.
    """
    import pandas as pd
    try:
        import h5py
    except ImportError:
        h5py = None

    n_copied = 0
    if h5py is not None:
        with h5py.File(src, "r") as fsrc, h5py.File(dst, "w") as fdst:
            for k, v in fsrc.attrs.items():    # pytables root attrs
                fdst.attrs[k] = v
            for name in fsrc:
                if BNB_KEY_RE.match(name):
                    continue
                fsrc.copy(name, fdst, name=name)
                n_copied += 1
        mode = "a"
    else:
        # fallback: copy through pandas. Correct, but every key is held in
        # memory one at a time, which hurts on multi-GB event frames.
        print("    (h5py not available: copying through pandas, needs more memory)")
        with pd.HDFStore(src, mode="r") as ssrc, pd.HDFStore(dst, mode="w") as sdst:
            for key in ssrc.keys():
                if BNB_KEY_RE.match(key):
                    continue
                fmt = getattr(ssrc.get_storer(key), "format_type", "fixed")
                sdst.put(key, ssrc[key], format=fmt)
                n_copied += 1
        mode = "a"
    with pd.HDFStore(dst, mode=mode) as store:
        for key, df in frames.items():
            store.put(key, df, format=formats.get(key, "fixed"))
    return n_copied


def out_path(inp, args):
    if args.output:
        return args.output
    base = os.path.basename(inp)
    if args.outdir:
        return os.path.join(args.outdir, base)
    stem, ext = os.path.splitext(base)
    return os.path.join(os.path.dirname(inp) or ".", f"{stem}{args.suffix}{ext}")


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dataframes", nargs="+", help="input cafpyana .df file(s)")
    ap.add_argument("--root", nargs="+", required=True,
                    help="beamqual_*_fom.root file(s) from bnb_fom_caf.py")
    ap.add_argument("-o", "--output", help="output path (one input only)")
    ap.add_argument("--outdir", help="write outputs into this directory, same basenames")
    ap.add_argument("--suffix", default="_fomfix",
                    help="default naming: <stem><suffix>.df beside the input")
    ap.add_argument("--tree", choices=["bnbinfo", "spill"], default="bnbinfo",
                    help="ROOT tree to join against (default bnbinfo: every spill)")
    ap.add_argument("--fom-column", default="FOM_best", help="corrected FOM in the ROOT file")
    ap.add_argument("--new-column", default="FOM_corrected", help="name of the added column")
    ap.add_argument("--extra-columns", default="fom_source,fom_status",
                    help="further ROOT columns to copy (comma separated, '' for none)")
    ap.add_argument("--tolerance", type=float, default=0.0,
                    help="[s] nearest-time matching window (default 0: exact stamp)")
    ap.add_argument("--inspect", action="store_true", help="print the bnb_* schema and exit")
    ap.add_argument("--dry-run", action="store_true", help="report the match, write nothing")
    ap.add_argument("--force", action="store_true", help="overwrite an existing output")
    a = ap.parse_args(argv)

    if a.output and len(a.dataframes) > 1:
        die("-o takes a single input; use --outdir for several")
    if a.output and a.outdir:
        die("-o and --outdir are mutually exclusive")
    extra = [c.strip() for c in a.extra_columns.split(",") if c.strip()]

    try:
        import pandas as pd
    except ImportError:
        die("pandas is required (with pytables and h5py for .df I/O)")

    if a.inspect:
        for path in a.dataframes:
            print(f"\n=== {path}")
            with pd.HDFStore(path, mode="r") as s:
                ks = bnb_keys(s)
                print(f"  bnb_* keys : {ks or 'NONE FOUND'}")
                print(f"  other keys : {[k for k in s.keys() if k not in ks]}")
                for k in ks[:2]:
                    df = s[k]
                    cols = flat_columns(df)
                    print(f"\n  --- {k}: {len(df):,} rows, index {df.index.names}")
                    print(f"      columns ({len(df.columns)}): {[str(c) for c in df.columns]}")
                    print(f"      time key : {pick(cols, SEC_CANDIDATES)} / "
                          f"{pick(cols, NSEC_CANDIDATES)}")
                    print(f"      FOM      : {pick(cols, FOM_CANDIDATES)}")
        return 0

    if a.outdir:
        os.makedirs(a.outdir, exist_ok=True)

    print(f"spill lookup, tree '{a.tree}':")
    stamp_lut, values, meta = load_lookup(a.root, a.tree, a.fom_column, extra)
    print(f"  {meta['n_unique']:,} unique spill stamps "
          f"({meta['n_dup']:,} duplicate stamps collapsed)")

    rc = 0
    for path in a.dataframes:
        dst = None if a.dry_run else out_path(path, a)
        print(f"\n=== {path}" + (f"\n    -> {dst}" if dst else "  (dry run)"))
        if dst:
            if os.path.abspath(dst) == os.path.abspath(path):
                die(f"output would overwrite the input ({dst})")
            if os.path.exists(dst) and not a.force:
                die(f"{dst} exists (use --force)")

        with pd.HDFStore(path, mode="r") as s:
            ks = bnb_keys(s)
            if not ks:
                print("  no bnb_* keys; skipped")
                rc = 1
                continue
            formats = {k: getattr(s.get_storer(k), "format_type", "fixed") for k in ks}
            frames = {k: s[k] for k in ks}

        fixed = {}
        for k in ks:
            df = frames[k]
            cols = flat_columns(df)
            stamp_df, how = frame_stamp(df, cols, a.tolerance)
            idx = match(stamp_df, stamp_lut, a.tolerance)
            report(k, df, stamp_df, how, idx, values, cols, a.fom_column, a.new_column)
            if not (idx >= 0).all():
                rc = 1
            fixed[k] = augment(df, idx, values, a.fom_column, a.new_column, extra)

        if dst:
            n_copied = write_copy(path, dst, fixed, formats)
            print(f"\n  wrote {dst}: {n_copied} key(s) copied, "
                  f"{len(fixed)} bnb_* key(s) corrected")
        else:
            print("\n  dry run: nothing written")

    if rc:
        print("\nsome rows did not match; see the report above", file=sys.stderr)
    return rc


if __name__ == "__main__":
    sys.exit(main())
