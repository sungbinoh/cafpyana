#!/usr/bin/env python3
"""Efficiency, purity and figures of merit of the NCDelta selection.

The two channels are reported separately and are NOT summed. They are tuned
to different operating points (see the "Tuning" section of cuts.py): 1g1p is
the pure channel and 1g0p the statistics channel, so a combined row is a
signal-weighted average of two deliberately different selections. It is
always worse than 1g1p and better than 1g0p, and is not a figure of merit
anything should be optimized against. Sum the raw mc_signal / mc_bkg columns
yourself if you need the total selected yield.

For each channel (1g1p, 1g0p):
  signal, background   selected slices, scaled to --pot
  mc_signal, mc_bkg    raw (unweighted) MC slice counts, for MC stat errors
  efficiency           selected signal slices / true signal neutrinos (mcnu)
  purity               signal / (signal + background)
  s_over_sqrt_b        signal / sqrt(background)
  s_over_sqrt_sb       signal / sqrt(signal + background)

Cut overrides come from a JSON file mapping a channel to CUTBASED keys, e.g.
  {"1g1p": {"trk_proton_like": false, "nu_score_min": 0.4},
   "1g0p": {"shw_dedx_min": 2.8}}
Keys not given keep their defaults from cuts.py; null disables a cut.

Use --split even/odd to evaluate on half of the files (by position in the
sorted input list). Tune on one half and check the result on the other.

Example:
  python evaluate.py -i data/ --pot 1e20 --config my_cuts.json --split odd
"""

import argparse
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cuts                            # noqa: E402
from dfio import expand_inputs, load   # noqa: E402


def channel_config(overrides, topo):
    """CUTBASED[topo] with overrides applied (None means the defaults)."""
    if not overrides or topo not in overrides:
        return None
    cfg = dict(cuts.CUTBASED[topo])
    unknown = set(overrides[topo]) - set(cfg)
    if unknown:
        raise KeyError("Unknown %s cut(s) in config: %s" % (topo, sorted(unknown)))
    cfg.update(overrides[topo])
    return cfg


def split_mask(index, split):
    """Files at even or odd positions in the input list."""
    if split == "all":
        return np.ones(len(index), dtype=bool)
    parity = np.asarray(index.get_level_values("file")) % 2
    return parity == (0 if split == "even" else 1)


def evaluate(slc, nu, scale=1.0, overrides=None):
    """DataFrame of selection metrics, one row per channel.

    The channels are kept separate on purpose; see the module docstring.

    `slc` and `nu` are the slice and true-neutrino frames from dfio.load;
    `scale` multiplies every slice and neutrino (e.g. target POT / MC POT).
    """
    cat = cuts.truth_category(slc)
    is_sig = (cat == 0).values
    n_true = scale * cuts.is_signal(nu).sum()

    rows = {}
    for topo in cuts.TOPOLOGIES:
        sel = cuts.selection(slc, topo, cutbased=channel_config(overrides, topo) or True).values
        mc_s = int((sel & is_sig).sum())
        mc_b = int((sel & ~is_sig).sum())
        rows[topo] = {"mc_signal": mc_s, "mc_bkg": mc_b}

    out = {}
    for name, r in rows.items():
        s, b = scale * r["mc_signal"], scale * r["mc_bkg"]
        out[name] = {
            "signal": s,
            "background": b,
            "mc_signal": r["mc_signal"],
            "mc_bkg": r["mc_bkg"],
            "efficiency": s / n_true if n_true > 0 else float("nan"),
            "purity": s / (s + b) if s + b > 0 else float("nan"),
            "s_over_sqrt_b": s / math.sqrt(b) if b > 0 else float("nan"),
            "s_over_sqrt_sb": s / math.sqrt(s + b) if s + b > 0 else float("nan"),
        }
    df = pd.DataFrame(out).T
    df.attrs["true_signal"] = n_true
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--inputs", nargs="+", required=True,
                    help="directories, .df files, glob patterns (quote them), or .list/.txt files")
    ap.add_argument("--pot", type=float, default=None,
                    help="scale to this POT (default: unscaled MC counts)")
    ap.add_argument("--config", default=None, help="JSON file of CUTBASED overrides per channel")
    ap.add_argument("--split", choices=("all", "even", "odd"), default="all",
                    help="evaluate on all files, or only even/odd-numbered files")
    ap.add_argument("--csv", default=None, help="also write the table to this CSV file")
    args = ap.parse_args()

    overrides = None
    if args.config:
        with open(args.config) as f:
            overrides = json.load(f)

    files = expand_inputs(args.inputs)
    frames, pot_by_file = load(files, keys=("slc", "mcnu"))
    slc, nu = frames["slc"], frames["mcnu"]
    if args.split != "all":
        slc = slc[split_mask(slc.index, args.split)]
        nu = nu[split_mask(nu.index, args.split)]
        pot_by_file = pot_by_file[split_mask(pot_by_file.index, args.split)]
        print("Split '%s': %d files, %d slices, %d true neutrinos"
              % (args.split, len(pot_by_file), len(slc), len(nu)))
    mc_pot = pot_by_file.sum()

    scale = 1.0
    if args.pot is not None:
        if mc_pot <= 0:
            sys.exit("MC POT is zero; cannot scale to --pot.")
        scale = args.pot / mc_pot

    table = evaluate(slc, nu, scale=scale, overrides=overrides)
    print("MC POT %.3g; true signal neutrinos (scaled): %.1f"
          % (mc_pot, table.attrs["true_signal"]))
    with pd.option_context("display.float_format", "{:.4g}".format):
        print(table.to_string())
    if args.csv:
        table.to_csv(args.csv)


if __name__ == "__main__":
    main()
