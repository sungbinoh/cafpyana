#!/usr/bin/env python3
"""
Chunked post-processing for SBND CAFs.

This mirrors ``postprocess_cut3.py`` but processes matching HDF5 keys (e.g.
``mcnu_0``, ``evt_pand_0``) independently so that multi-terabyte inputs stay
memory-friendly. Each chunk is handled in its own worker, scaled using the
global POT/livetime totals, and written back with the same key suffix.

Example:
    python -m analysis_village.numuincl.postprocess_chunked \\
        --apply-cuts --spine --max-workers 8 --output-suffix postprocess_cut3
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Sequence, Tuple

import argparse
import gc
import json
import multiprocessing as mp
import os
import sys
from time import time

import h5py
import numpy as np
import pandas as pd

# -- Environment bootstrap ----------------------------------------------------

CAFPYANA_WD = os.environ.get(
    "CAFPYANA_WD", "/exp/sbnd/app/users/brindenc/develop/cafpyana"
)
os.environ["CAFPYANA_WD"] = CAFPYANA_WD

if CAFPYANA_WD not in sys.path:
    sys.path.insert(0, CAFPYANA_WD)
    sys.path.insert(0, os.path.join(CAFPYANA_WD, "pyanalib"))

SBNDANA_DIR = os.path.join(
    CAFPYANA_WD, "analysis_village", "numuincl", "sbnd"
)
if SBNDANA_DIR not in sys.path:
    sys.path.append(SBNDANA_DIR)

from sbnd.cafclasses.slice import CAFSlice  # type: ignore
from sbnd.cafclasses.interaction import CAFInteraction  # type: ignore
from sbnd.cafclasses.nu import NU  # type: ignore
from sbnd.constants import *  # noqa: F401,F403
from sbnd.numu.numu_constants import *  # noqa: F401,F403
from sbnd.detector.definitions import *  # noqa: F401,F403

from naming import PAND_CUTS, PAND_CUTS_CONT, SPINE_CUTS  # type: ignore  # noqa: F401

# -- Dataclasses --------------------------------------------------------------


def _default_key_list() -> List[str]:
    return []


@dataclass
class ChunkKeySet:
    suffix: str
    mcnu: List[str] = field(default_factory=_default_key_list)
    mc_hdr: List[str] = field(default_factory=_default_key_list)
    mc_pandora: List[str] = field(default_factory=_default_key_list)
    mc_spine: List[str] = field(default_factory=_default_key_list)
    offbeam_hdr: List[str] = field(default_factory=_default_key_list)
    offbeam_pandora: List[str] = field(default_factory=_default_key_list)
    offbeam_spine: List[str] = field(default_factory=_default_key_list)
    data_hdr: List[str] = field(default_factory=_default_key_list)
    data_pandora: List[str] = field(default_factory=_default_key_list)
    data_spine: List[str] = field(default_factory=_default_key_list)

    def has_any(self) -> bool:
        return any(
            len(getattr(self, attr))
            for attr in (
                "mcnu",
                "mc_pandora",
                "mc_spine",
                "offbeam_pandora",
                "offbeam_spine",
                "data_pandora",
                "data_spine",
            )
        )


@dataclass
class Totals:
    pot_mc: float = 0.0
    pot_data: float = 0.0
    livetime_offbeam: float = 0.0
    livetime_data: float = 0.0

    def as_dict(self) -> Dict[str, float]:
        return asdict(self)


@dataclass
class ProcessingConfig:
    apply_cuts: bool
    process_mcnu: bool
    process_pandora: bool
    process_spine: bool
    process_weights: bool
    livetime_data_override: Optional[float]
    max_workers: int
    output_suffix: str
    copy_headers: bool


@dataclass
class InputFiles:
    data_dir: str
    mc_fname: str
    offbeam_fname: str
    data_fname: str

    @property
    def mc_path(self) -> str:
        return os.path.join(self.data_dir, self.mc_fname)

    @property
    def offbeam_path(self) -> str:
        return os.path.join(self.data_dir, self.offbeam_fname)

    @property
    def data_path(self) -> str:
        return os.path.join(self.data_dir, self.data_fname)


@dataclass
class OutputPaths:
    mc: str
    data: str


@dataclass
class ChunkResult:
    suffix: str
    logs: List[str]
    duration: float


# -- Helpers -----------------------------------------------------------------


def human_gb(num_bytes: float) -> float:
    return num_bytes / (1024**3)


def extract_suffix(key: str) -> str:
    if "_" not in key:
        return "full"
    return key.split("_")[-1]


def categorize_key(key: str) -> str:
    if "hdr" in key:
        return "hdr"
    if "mcnu" in key:
        return "mcnu"
    if "pand" in key:
        return "pandora"
    if "trk" in key or ("evt" in key and "pand" not in key):
        return "spine"
    return "other"


def sum_column(path: str, key: str, column: str) -> float:
    if column is None:
        return 0.0
    try:
        df = pd.read_hdf(path, key=key, columns=[column])
    except (KeyError, ValueError):
        df = pd.read_hdf(path, key=key)
    if column not in df.columns:
        return 0.0
    return float(df[column].sum())


def discover_chunks(
    files: InputFiles,
    cfg: ProcessingConfig,
    max_mc_keys: Optional[int] = None,
    max_offbeam_keys: Optional[int] = None,
) -> Tuple[Totals, Dict[str, ChunkKeySet]]:
    totals = Totals()
    chunks: Dict[str, ChunkKeySet] = {}

    def ensure_chunk(key: str) -> ChunkKeySet:
        suffix = extract_suffix(key)
        if suffix not in chunks:
            chunks[suffix] = ChunkKeySet(suffix=suffix)
        return chunks[suffix]

    # MC
    with h5py.File(files.mc_path, "r") as handle:
        for idx, key in enumerate(handle.keys()):
            if max_mc_keys is not None and idx >= max_mc_keys:
                break
            cat = categorize_key(key)
            chunk = ensure_chunk(key)
            if cat == "pandora":
                chunk.mc_pandora.append(key)
            elif cat == "spine":
                chunk.mc_spine.append(key)
            elif cat == "mcnu":
                chunk.mcnu.append(key)
            elif cat == "hdr":
                chunk.mc_hdr.append(key)
                totals.pot_mc += sum_column(files.mc_path, key, "pot")

    # Offbeam
    with h5py.File(files.offbeam_path, "r") as handle:
        for idx, key in enumerate(handle.keys()):
            if max_offbeam_keys is not None and idx >= max_offbeam_keys:
                break
            cat = categorize_key(key)
            chunk = ensure_chunk(key)
            if cat == "pandora":
                chunk.offbeam_pandora.append(key)
            elif cat == "spine":
                chunk.offbeam_spine.append(key)
            elif cat == "hdr":
                chunk.offbeam_hdr.append(key)
                totals.livetime_offbeam += sum_column(
                    files.offbeam_path, key, "noffbeambnb"
                )

    # Data
    with h5py.File(files.data_path, "r") as handle:
        for key in handle.keys():
            cat = categorize_key(key)
            chunk = ensure_chunk(key)
            if cat == "pandora":
                chunk.data_pandora.append(key)
            elif cat == "spine":
                chunk.data_spine.append(key)
            elif cat == "hdr":
                chunk.data_hdr.append(key)
                totals.pot_data += sum_column(files.data_path, key, "pot")
                totals.livetime_data += sum_column(
                    files.data_path, key, "noffbeambnb"
                )

    if totals.livetime_data == 0.0 and cfg.livetime_data_override is not None:
        totals.livetime_data = cfg.livetime_data_override

    # Drop empty chunks
    chunks = {k: v for k, v in chunks.items() if v.has_any()}
    return totals, chunks


def write_hdf(df: pd.DataFrame, path: str, key: str, lock: mp.synchronize.Lock):
    if df is None:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with lock:
        df.to_hdf(path, key=key, mode="a", format="table")


def _log_save(logs: List[str], label: str, df: pd.DataFrame, suffix: str):
    if df is None:
        return
    size_gb = human_gb(df.memory_usage(deep=True).sum())
    logs.append(f"Saved {label}_{suffix}: {size_gb:.4f} GB")


def load_slice_sequence(
    path: str, keys: Sequence[str], **kwargs
) -> Optional[CAFSlice]:
    obj: Optional[CAFSlice] = None
    for key in keys:
        current = CAFSlice.load(path, key=key, **kwargs)
        if obj is None:
            obj = current
        else:
            obj.combine(current)
            del current
            gc.collect()
    return obj


def load_interaction_sequence(
    path: str, keys: Sequence[str], **kwargs
) -> Optional[CAFInteraction]:
    obj: Optional[CAFInteraction] = None
    for key in keys:
        current = CAFInteraction.load(path, key=key, **kwargs)
        if obj is None:
            obj = current
        else:
            obj.combine(current)
            del current
            gc.collect()
    return obj


def load_nu_sequence(path: str, keys: Sequence[str]) -> Optional[NU]:
    obj: Optional[NU] = None
    for key in keys:
        current = NU.load(path, key=key)
        if obj is None:
            obj = current
        else:
            obj.combine(current)
            del current
            gc.collect()
    return obj


def process_chunk(
    chunk: ChunkKeySet,
    files: InputFiles,
    cfg: ProcessingConfig,
    totals: Totals,
    outputs: OutputPaths,
    lock: mp.synchronize.Lock,
) -> ChunkResult:
    start = time()
    logs: List[str] = []

    mcnu_obj: Optional[NU] = None
    if cfg.process_mcnu and chunk.mcnu:
        mcnu_obj = load_nu_sequence(files.mc_path, chunk.mcnu)
        if mcnu_obj is not None:
            mcnu_obj.scale_to_pot(totals.pot_data, sample_pot=totals.pot_mc)
            mcnu_obj.add_fv()
            mcnu_obj.add_av()
            mcnu_obj.cut_muon(cut=False, min_ke=0.1)
            mcnu_obj.cut_fv(cut=False)
            mcnu_obj.cut_cosmic(cut=False)
            mcnu_obj.cut_cont(cut=False)

    slc: Optional[CAFSlice] = None
    slc_data: Optional[CAFSlice] = None
    slc_signal: Optional[CAFSlice] = None

    if cfg.process_pandora and chunk.mc_pandora:
        slc = load_slice_sequence(
            files.mc_path, chunk.mc_pandora, pot=totals.pot_mc
        )
        if slc is not None:
            slc.scale_to_pot(totals.pot_data, sample_pot=totals.pot_mc)
            gc.collect()

            if chunk.offbeam_pandora:
                offbeam = load_slice_sequence(
                    files.offbeam_path,
                    chunk.offbeam_pandora,
                    livetime=totals.livetime_offbeam,
                )
                if offbeam is not None:
                    offbeam.scale_to_livetime(
                        totals.livetime_data,
                        sample_livetime=totals.livetime_offbeam,
                    )
                    slc.combine(offbeam, duplicate_ok=True)
                    del offbeam
                    gc.collect()

            if chunk.data_pandora:
                slc_data = load_slice_sequence(
                    files.data_path, chunk.data_pandora
                )

            if mcnu_obj is not None:
                mcnu_obj = slc.set_mcnu_containment(mcnu_obj)
                mcnu_obj.add_event_type("pandora")

            slc_list = [obj for obj in (slc, slc_data) if obj is not None]
            for idx, sample in enumerate(slc_list):
                sample.clean(dummy_vals=[-9999, -999, 999, 9999, -5])
                sample.add_has_muon()
                sample.add_in_av()
                sample.add_in_fv()
                if idx == 0:
                    sample.add_event_type()
                    slc_signal = sample.copy()
                    is_signal = np.isin(sample.data.truth.event_type, [0, 1])
                    slc_signal.data = sample.data[is_signal]
                sample.cut_flashmatch(cut=cfg.apply_cuts)
                sample.cut_fv(cut=cfg.apply_cuts)
                sample.cut_muon(cut=cfg.apply_cuts, min_ke=0.1)
                sample.cut_cosmic(
                    cut=cfg.apply_cuts,
                    fmatch_score=320,
                    nu_score=0.5,
                    use_opt0=True,
                    use_isclearcosmic=False,
                )
                sample.cut_lowz(cut=cfg.apply_cuts, z_max=6, include_start=True)
                sample.cut_is_cont(cut=False)

    inter: Optional[CAFInteraction] = None
    inter_data: Optional[CAFInteraction] = None
    inter_signal: Optional[CAFInteraction] = None

    if cfg.process_spine and chunk.mc_spine:
        inter = load_interaction_sequence(
            files.mc_path, chunk.mc_spine, pot=totals.pot_mc
        )
        if inter is not None:
            inter.scale_to_pot(totals.pot_data, sample_pot=totals.pot_mc)
            gc.collect()

            if chunk.offbeam_spine:
                inter_offbeam = load_interaction_sequence(
                    files.offbeam_path,
                    chunk.offbeam_spine,
                    livetime=totals.livetime_offbeam,
                )
                if inter_offbeam is not None:
                    inter_offbeam.scale_to_livetime(
                        totals.livetime_data,
                        sample_livetime=totals.livetime_offbeam,
                    )
                    inter.combine(inter_offbeam, duplicate_ok=True)
                    del inter_offbeam
                    gc.collect()

            if chunk.data_spine:
                inter_data = load_interaction_sequence(
                    files.data_path, chunk.data_spine
                )

            if mcnu_obj is not None:
                mcnu_obj = inter.set_mcnu_containment(mcnu_obj)
                mcnu_obj.add_event_type("spine")

            inter_list = [obj for obj in (inter, inter_data) if obj is not None]
            for idx, sample in enumerate(inter_list):
                sample.clean(
                    dummy_vals=[-9999, -999, 999, 9999, -5, np.inf, -np.inf]
                )
                sample.add_in_fv()
                sample.add_in_av()
                if idx == 0:
                    sample.add_event_type()
                    inter_signal = sample.copy()
                    is_signal = np.isin(sample.data.truth.event_type, [0, 1])
                    inter_signal.data = sample.data[is_signal]
                sample.cut_time_contained(cut=cfg.apply_cuts)
                sample.cut_cosmic(cut=cfg.apply_cuts)
                sample.cut_fv(cut=cfg.apply_cuts)
                sample.cut_muon(cut=cfg.apply_cuts, min_ke=0.1)
                sample.cut_cosmic_score(cut=cfg.apply_cuts, score=102.35)
                sample.cut_lowz(cut=cfg.apply_cuts, z_max=6, include_start=True)
                sample.cut_start_dedx(cut=cfg.apply_cuts, dedx=4.17)
                sample.cut_is_cont(cut=False)

    if cfg.process_weights:
        if cfg.process_pandora and slc is not None and chunk.mc_pandora:
            slc.add_universe_weights(
                files.mc_path, keys=chunk.mc_pandora, duplicate_ok=False
            )
            if slc_signal is not None:
                slc_signal.add_universe_weights(
                    files.mc_path, keys=chunk.mc_pandora, duplicate_ok=False
                )
        if cfg.process_spine and inter is not None and chunk.mc_spine:
            inter.add_universe_weights(
                files.mc_path, keys=chunk.mc_spine, duplicate_ok=False
            )
            if inter_signal is not None:
                inter_signal.add_universe_weights(
                    files.mc_path, keys=chunk.mc_spine, duplicate_ok=False
                )

    if cfg.process_pandora and slc is not None:
        write_hdf(slc.data, outputs.mc, f"pandora_{chunk.suffix}", lock)
        _log_save(logs, "pandora", slc.data, chunk.suffix)
        if slc_signal is not None:
            write_hdf(
                slc_signal.data,
                outputs.mc,
                f"pandora_signal_{chunk.suffix}",
                lock,
            )
            _log_save(logs, "pandora_signal", slc_signal.data, chunk.suffix)
        if slc_data is not None:
            write_hdf(
                slc_data.data,
                outputs.data,
                f"pandora_{chunk.suffix}",
                lock,
            )
            _log_save(logs, "pandora_data", slc_data.data, chunk.suffix)

    if cfg.process_spine and inter is not None:
        write_hdf(inter.data, outputs.mc, f"spine_{chunk.suffix}", lock)
        _log_save(logs, "spine", inter.data, chunk.suffix)
        if inter_signal is not None:
            write_hdf(
                inter_signal.data,
                outputs.mc,
                f"spine_signal_{chunk.suffix}",
                lock,
            )
            _log_save(logs, "spine_signal", inter_signal.data, chunk.suffix)
        if inter_data is not None:
            write_hdf(
                inter_data.data,
                outputs.data,
                f"spine_{chunk.suffix}",
                lock,
            )
            _log_save(logs, "spine_data", inter_data.data, chunk.suffix)

    if cfg.process_mcnu and mcnu_obj is not None:
        write_hdf(mcnu_obj.data, outputs.mc, f"mcnu_{chunk.suffix}", lock)
        _log_save(logs, "mcnu", mcnu_obj.data, chunk.suffix)

    if cfg.copy_headers:
        for key in chunk.mc_hdr:
            df = pd.read_hdf(files.mc_path, key=key)
            write_hdf(df, outputs.mc, key, lock)
            _log_save(logs, key, df, chunk.suffix)
        for key in chunk.data_hdr:
            df = pd.read_hdf(files.data_path, key=key)
            write_hdf(df, outputs.data, key, lock)
            _log_save(logs, key, df, chunk.suffix)

    duration = time() - start
    return ChunkResult(suffix=chunk.suffix, logs=logs, duration=duration)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Chunked CAF post-processing with parallel execution."
    )
    parser.add_argument(
        "--data-dir",
        default="/pnfs/sbnd/persistent/users/brindenc/numuincl_analysis/v10_06_00/pandora",
        help="Base directory that holds the CAF HDF5 files.",
    )
    parser.add_argument(
        "--mc-fname",
        default="mc_syst/mc_tinypand_fullsyst_cut3.df",
        help="Relative MC HDF5 file name inside data-dir.",
    )
    parser.add_argument(
        "--offbeam-fname",
        default="data_offbeam2.df",
        help="Relative offbeam HDF5 file name.",
    )
    parser.add_argument(
        "--data-fname",
        default="data_dev2.df",
        help="Relative data HDF5 file name.",
    )
    parser.add_argument(
        "--output-suffix",
        default="postprocess_chunked",
        help="Tag appended to the output HDF5 file names.",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=max(1, mp.cpu_count() // 2),
        help="Number of worker processes.",
    )
    parser.add_argument(
        "--apply-cuts",
        action="store_true",
        help="Apply the standard selection cuts per chunk.",
    )
    parser.add_argument(
        "--no-mcnu",
        action="store_true",
        help="Skip MCNU processing.",
    )
    parser.add_argument(
        "--no-pandora",
        action="store_true",
        help="Skip Pandora slice processing.",
    )
    parser.add_argument(
        "--spine",
        action="store_true",
        help="Enable SPINE processing.",
    )
    parser.add_argument(
        "--weights",
        action="store_true",
        help="Re-attach universe weights per chunk.",
    )
    parser.add_argument(
        "--livetime-data",
        type=float,
        default=9.51e5,
        help="Override for data livetime if headers lack the value.",
    )
    parser.add_argument(
        "--metadata-out",
        default=None,
        help="Optional path to dump POT/livetime totals as JSON.",
    )
    parser.add_argument(
        "--copy-headers",
        action="store_true",
        help="Copy header keys into the output file alongside processed data.",
    )
    parser.add_argument(
        "--max-mc-keys",
        type=int,
        default=None,
        help="Limit MC keys to this many (debug/testing).",
    )
    parser.add_argument(
        "--max-offbeam-keys",
        type=int,
        default=None,
        help="Limit offbeam keys to this many (debug/testing).",
    )
    return parser.parse_args()


def main() -> None:
    global args  # pylint: disable=global-statement
    args = parse_args()

    files = InputFiles(
        data_dir=args.data_dir,
        mc_fname=args.mc_fname,
        offbeam_fname=args.offbeam_fname,
        data_fname=args.data_fname,
    )

    cfg = ProcessingConfig(
        apply_cuts=args.apply_cuts,
        process_mcnu=not args.no_mcnu,
        process_pandora=not args.no_pandora,
        process_spine=args.spine,
        process_weights=args.weights,
        livetime_data_override=args.livetime_data,
        max_workers=args.max_workers,
        output_suffix=args.output_suffix,
        copy_headers=args.copy_headers,
    )

    outputs = OutputPaths(
        mc=os.path.join(
            files.data_dir,
            f"{os.path.splitext(files.mc_fname)[0]}_{cfg.output_suffix}.df",
        ),
        data=os.path.join(
            files.data_dir,
            f"{os.path.splitext(files.data_fname)[0]}_{cfg.output_suffix}.df",
        ),
    )

    print("Scanning keys and computing global POT/livetime...")
    totals, chunks = discover_chunks(
        files,
        cfg,
        max_mc_keys=args.max_mc_keys,
        max_offbeam_keys=args.max_offbeam_keys,
    )
    print(f"Found {len(chunks)} chunk suffixes.")
    print(
        f"Global POT/LT: POT_MC={totals.pot_mc:.2e}, POT_DATA={totals.pot_data:.2e}, "
        f"LT_OFFBEAM={totals.livetime_offbeam:.2e}, LT_DATA={totals.livetime_data:.2e}"
    )
    if args.metadata_out:
        with open(args.metadata_out, "w", encoding="utf-8") as f:
            json.dump(totals.as_dict(), f, indent=2)

    suffixes = sorted(chunks.keys(), key=lambda x: (len(x), x))
    ctx = mp.get_context("spawn")
    lock = ctx.Lock()
    futures = []
    start = time()
    with ProcessPoolExecutor(max_workers=cfg.max_workers, mp_context=ctx) as pool:
        for suffix in suffixes:
            futures.append(
                pool.submit(
                    process_chunk,
                    chunks[suffix],
                    files,
                    cfg,
                    totals,
                    outputs,
                    lock,
                )
            )
        for future in as_completed(futures):
            result = future.result()
            for line in result.logs:
                print(line)
            print(f"Chunk {result.suffix} finished in {result.duration:.2f}s")
    total_time = time() - start
    print(f"All chunks complete in {total_time/60:.2f} minutes")


if __name__ == "__main__":
    main()

