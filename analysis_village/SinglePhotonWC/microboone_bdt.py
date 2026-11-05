"""The four MicroBooNE single-photon TMVA BDTs (numu, other, ncpi0, nue),
scored on Wire-Cell candidate rows (make_wc_candidates.py columns).

The models are the original MicroBooNE GradBoost weight files,
single_photon_<name>_bdt_final.xml, not SBND-trained; they are not in git.
Pass their directory, or set SP_UBOONE_BDT_DIR.

Evaluation follows abhatfnal/sbnd-wirecell-single-photon
score_microboone_bdts.py, which agrees with TMVA::Reader to 5e-7:
EvaluateMVA = tanh(F), F = sum of boostWeight * leaf res, going left when
value <= cut. The score is MicroBooNE's Wire-Cell convention,
log10((1 + mva) / (1 - mva)) = 2 F / ln 10, which the selection thresholds
are in (the paper's nue score axis spans -8 to 4; tanh saturates at -1).
Inputs are read by their XML name from the candidate columns, except
reco_nuvtxY -> kine_nu_y_corr. A score exists only if shw_sp_filled == 1
and every input is finite (no defaults); otherwise it is NaN.

  df = df.join(score(df, bdt_dir))     # adds ub_numu, ub_other, ub_ncpi0, ub_nue
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd

NAMES = ["numu", "other", "ncpi0", "nue"]
EXPECTED_NVARS = {"numu": 73, "other": 146, "ncpi0": 45, "nue": 56}
COLUMN_OF = {"reco_nuvtxY": "kine_nu_y_corr"}
DEFAULT_DIR = os.environ.get("SP_UBOONE_BDT_DIR", "")


class TmvaGradBdt:
    """One TMVA GradBoost forest as flat node arrays, evaluated vectorised."""

    def __init__(self, xml_path):
        root = ET.parse(xml_path).getroot()
        self.variables = [v.get("Expression") for v in
                          sorted(root.find("Variables").findall("Variable"), key=lambda v: int(v.get("VarIndex")))]
        self.trees = []
        for bt in root.find("Weights").findall("BinaryTree"):
            ivar, cut, left, right, res = [], [], [], [], []

            def add(node):
                i = len(ivar)
                ivar.append(-1); cut.append(0.); left.append(-1); right.append(-1); res.append(0.)
                kids = {c.get("pos"): c for c in node if c.tag == "Node"}
                if not kids:
                    res[i] = float(node.get("res"))
                    return i
                assert node.get("cType") == "1", "only cType=1 (go right if value > cut) is handled"
                ivar[i], cut[i] = int(node.get("IVar")), float(node.get("Cut"))
                left[i] = add(kids["l"])
                right[i] = add(kids["r"])
                return i

            add(bt.find("Node"))
            self.trees.append((float(bt.get("boostWeight")), np.array(ivar), np.array(cut),
                               np.array(left), np.array(right), np.array(res)))

    def evaluate(self, X):
        """X: (n, nvar) array in XML variable order -> F (n,); EvaluateMVA is tanh(F)."""
        X = np.asarray(X, dtype=np.float64)
        rows = np.arange(len(X))
        raw = np.zeros(len(X))
        for bw, ivar, cut, left, right, res in self.trees:
            node = np.zeros(len(X), dtype=int)
            active = ivar[node] >= 0
            while active.any():
                n = node[active]
                go_left = X[rows[active], ivar[n]] <= cut[n]
                node[active] = np.where(go_left, left[n], right[n])
                active = ivar[node] >= 0
            raw += bw * res[node]
        return raw


_MODELS = {}


def models(bdt_dir=None):
    bdt_dir = bdt_dir or DEFAULT_DIR
    if not bdt_dir:
        raise ValueError("No BDT directory: pass bdt_dir or set SP_UBOONE_BDT_DIR")
    if bdt_dir not in _MODELS:
        M = {n: TmvaGradBdt(os.path.join(bdt_dir, "single_photon_%s_bdt_final.xml" % n)) for n in NAMES}
        for n, m in M.items():
            assert len(m.variables) == EXPECTED_NVARS[n], (n, len(m.variables))
        _MODELS[bdt_dir] = M
    return _MODELS[bdt_dir]


def input_columns(bdt_dir=None):
    """Candidate columns each BDT reads, in XML order."""
    return {n: [COLUMN_OF.get(v, v) for v in m.variables] for n, m in models(bdt_dir).items()}


def score(df, bdt_dir=None):
    """The four scores per candidate row, columns ub_<name>; NaN where not
    evaluable."""
    out = pd.DataFrame(index=df.index)
    filled = (df.shw_sp_filled == 1).values
    for n, m in models(bdt_dir).items():
        cols = input_columns(bdt_dir)[n]
        missing = [c for c in cols if c not in df.columns]
        if missing:
            raise KeyError("%s BDT inputs missing from the dataframe: %s" % (n, missing[:5]))
        X = df[cols].to_numpy(dtype=np.float64)
        ok = filled & np.isfinite(X).all(axis=1)
        s = np.full(len(df), np.nan)
        if ok.any():
            s[ok] = 2 * m.evaluate(X[ok]) / np.log(10)
        out["ub_" + n] = s
    return out
