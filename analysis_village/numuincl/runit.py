#!/usr/bin/env python3
"""
Submit run_df_maker jobs described in yamls/mcbnb.yaml.
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Optional
from sbnd.general.utils import convert_pnfs_to_xroot

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
NUMU_DIR = Path(__file__).resolve().parent
CONFIG_DIR = NUMU_DIR / "configs"
GENERATED_CONFIG_DIR = CONFIG_DIR / "generated"
DEFAULT_YAML = NUMU_DIR / "yamls" / "mcbnb.yaml"
RUN_DF_MAKER = REPO_ROOT / "run_df_maker.py"

FLAG_MAP = {
    "ncpu": "-ncpu",
    "ngrid": "-ngrid",
    "nfile": "-nfile",
    "split": "-split",
}


def main() -> None:
    os.system("source /exp/sbnd/app/users/brindenc/develop/cafpyana/setup.sh") # Source the setup script

    parser = argparse.ArgumentParser(description="Submit makedf jobs from YAML.")
    parser.add_argument(
        "-y", "--yaml", default=str(DEFAULT_YAML), help="Path to mcbnb.yaml (default: %(default)s)"
    )
    parser.add_argument("--only", nargs="+", help="Only run selected job names.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running.")
    args = parser.parse_args()

    yaml_path = Path(args.yaml).resolve()
    config = load_config(yaml_path)
    jobs = config["jobs"]
    input_dir = convert_pnfs_to_xroot(config.get("input_dir"))
    output_dir = convert_pnfs_to_xroot(config.get("output_dir"))
    selected = set(args.only) if args.only else None

    GENERATED_CONFIG_DIR.mkdir(parents=True, exist_ok=True)

    for job in jobs:
        name = job.get("name")
        if not name:
            print("Skipping unnamed job entry.")
            continue
        if not job.get("enabled", True):
            continue
        if selected and name not in selected:
            continue
        run_job(
            job,
            yaml_path=yaml_path,
            dry_run=args.dry_run,
            input_dir=input_dir,
            output_dir=output_dir,
        )

def load_config(path: Path) -> Dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"YAML file not found: {path}")
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    jobs = data.get("jobs")
    if not isinstance(jobs, list):
        raise SystemExit("mcbnb.yaml must define a top-level 'jobs' list.")
    return data


def file_exists(path: str) -> bool:
    """Check if a file exists, handling both local paths and XRootD URLs."""
    if path.startswith("root://"):
        # Extract host:port and path from XRootD URL
        # Format: root://host:port/path/to/file
        parts = path.split("/", 4)
        if len(parts) >= 4:
            host_port = parts[2]  # e.g., "fndcadoor.fnal.gov:1094"
            xrootd_path = "/" + parts[3]  # Get the path part
            # Use xrdfs to check if file exists
            result = subprocess.run(
                ["xrdfs", host_port, "stat", xrootd_path],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        return False
    else:
        # Local filesystem path
        return Path(path).exists()


def run_job(
    job: Dict[str, Any],
    yaml_path: Path,
    dry_run: bool,
    input_dir: Optional[str],
    output_dir: Optional[str],
) -> None:
    slug = slugify(job["name"])
    cfg_path = GENERATED_CONFIG_DIR / f"{slug}.py"
    cfg_path.write_text(render_config(job))

    # Get output path to check if file already exists
    output = job.get("output")
    if not output:
        raise SystemExit(f"Job '{job['name']}' missing output field.")
    output = output.format(name=job["name"], slug=slug)
    output = resolve_path(output, output_dir)
    output = output + ".df"  # run_df_maker adds .df extension
    
    # Check if output file already exists
    if file_exists(output):
        print("=" * 80)
        print(f"[{job['name']}] SKIPPED - Output file already exists:")
        print(f"  {output}")
        print("=" * 80)
        return

    cmd = build_command(job, cfg_path, slug, input_dir, output_dir)
    print("=" * 80)
    print(f"[{job['name']}] {' '.join(shlex.quote(part) for part in cmd)}")
    print(f"  Config: {cfg_path}")
    if dry_run:
        _print_settings(job)
        print("  Config preview:")
        with cfg_path.open() as handle:
            for i, line in enumerate(handle):
                print(f"    {line.rstrip()}")
                if i >= 40:
                    print("    ...")
                    break
        return

    env = os.environ.copy()
    env.setdefault("MAKEDF1MUX_YAML", str(yaml_path))
    env.setdefault("MAKEDF1MUX_JOB", job["name"])
    subprocess.run(cmd, check=True, env=env)


def render_config(job: Dict[str, Any]) -> str:
    settings = job.get("settings", {})
    base_config = job.get("base_config", "analysis_village.numuincl.configs.base")
    dfs_override = job.get("dfs")
    names_override = job.get("names")
    lines = [
        "from importlib import import_module",
        "from analysis_village.numuincl import makedf1muX as maker",
        "",
        f"# Auto-generated for job '{job['name']}'",
    ]
    for key, value in settings.items():
        literal = python_literal(value)
        if key == "UPDATE_RECOMB":
            lines.append(f"maker.set_update_recomb({literal})")
        else:
            lines.append(f"maker.{key} = {literal}")
    if settings:
        lines.append("maker.apply_setting_dependencies()")
    lines.extend(
        [
            "",
            f"_base = import_module('{module_name(base_config)}')",
        ]
    )
    if dfs_override is not None:
        lines.append("DFS = [")
        for entry in dfs_override:
            lines.append(f"    {format_callable(entry)},")
        lines.append("]")
    else:
        lines.append("DFS = _base.DFS")
    if names_override is not None:
        lines.append(f"NAMES = {names_override!r}")
    else:
        lines.append("NAMES = _base.NAMES")
    # Only set PREPROCESS if INCLUDE_WEIGHTS is True
    # lines.append("if maker.INCLUDE_WEIGHTS:")
    # lines.append("    from preprocess import Script")
    # lines.append('    PREPROCESS = [Script("/exp/sbnd/app/users/gputnam/Ar23-knobs/update_reweight_anywhere.sh")]')
    # lines.append("else:")
    # lines.append("    PREPROCESS = []")
    return "\n".join(lines) + "\n"


def build_command(
    job: Dict[str, Any],
    cfg_path: Path,
    slug: str,
    input_dir: Optional[str],
    output_dir: Optional[str],
) -> list[str]:
    output = job.get("output")
    if not output:
        raise SystemExit(f"Job '{job['name']}' missing output field.")
    output = output.format(name=job["name"], slug=slug)
    output = resolve_path(output, output_dir)
    cmd = [sys.executable, str(RUN_DF_MAKER), "-c", str(cfg_path), "-o", output]

    if job.get("input_list"):
        input_list = resolve_path(job["input_list"], input_dir)
        cmd.extend(["-l", input_list])
    elif job.get("input_files"):
        inputs = job["input_files"]
        if isinstance(inputs, (list, tuple)):
            resolved = [resolve_path(item, input_dir) for item in inputs]
            input_str = ",".join(resolved)
        else:
            input_str = resolve_path(inputs, input_dir)
        cmd.extend(["-i", input_str])
    else:
        raise SystemExit(f"Job '{job['name']}' requires input_list or input_files.")

    args = job.get("args", {})
    for key, flag in FLAG_MAP.items():
        if key in args:
            cmd.extend([flag, str(args[key])])
    if extra := args.get("extra"):
        cmd.extend(extra)
    return cmd


def slugify(name: str) -> str:
    chars = []
    for ch in name.strip():
        if ch.isalnum():
            chars.append(ch.lower())
        elif ch in ("-", "_"):
            chars.append(ch)
        else:
            chars.append("-")
    slug = "".join(chars).strip("-")
    return slug or "job"


def module_name(value: str) -> str:
    path = Path(value)
    if path.suffix == ".py":
        rel = path if path.is_absolute() else (REPO_ROOT / value)
        rel = rel.resolve().relative_to(REPO_ROOT)
        return ".".join(rel.with_suffix("").parts)
    return value


def python_literal(value: Any) -> str:
    if isinstance(value, bool):
        return "True" if value else "False"
    if isinstance(value, str):
        return repr(value)
    return str(value)


def resolve_path(value: str, base: Optional[str]) -> str:
    if not value:
        return value
    # If value is already an XRootD URL, return it as-is
    if value.startswith("root://"):
        return value
    # If base is an XRootD URL, append the relative path to it
    if base and base.startswith("root://"):
        # Remove trailing slash from base if present, then append value
        base_clean = base.rstrip("/")
        value_clean = value.lstrip("/")
        return f"{base_clean}/{value_clean}"
    # Normal filesystem path handling
    path = Path(value)
    if path.is_absolute() or not base:
        return str(path)
    return str(Path(base) / path)


def _print_settings(job: Dict[str, Any]) -> None:
    settings = job.get("settings", {})
    if not settings:
        return
    print("  makedf1muX overrides:")
    for key, value in settings.items():
        print(f"    {key} = {value!r}")


def format_callable(path: str) -> str:
    module, attr = path.rsplit(".", 1)
    return f"getattr(import_module('{module}'), '{attr}')"


if __name__ == "__main__":
    main()

