#!/usr/bin/env python3
"""Install benchmark reference-data symlinks from a dataset inventory.

The current workflow links files from an already-populated inventory directory:

    by_dataset/<dataset>/reference_inputs.tsv
    by_dataset/<dataset>/reference_files/<file>

The `install` command also knows how to unpack local
`<dataset>.ref_data.tar.gz` payloads before linking. Download support can be
added later by extending `ensure_reference_payloads`.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_INVENTORY_DIR = Path(
    "/home/unix/bhaas/projects/LRAA_PAPER_Analyses/"
    "benchmark_reference_inventory/by_dataset"
)

REGIME_DIRS = {
    "DENOVO": "DENOVO_ID",
    "QUANT_ONLY": "QUANT_ONLY",
    "REF_GUIDED": "REF_Guided",
}


@dataclass(frozen=True, order=True)
class LinkAction:
    source: Path
    dest: Path
    dataset: str
    regime: str
    role: str


@dataclass
class LinkStats:
    created: int = 0
    unchanged: int = 0
    replaced: int = 0
    missing_source: int = 0
    missing_target_dir: int = 0
    conflicts: int = 0
    would_create: int = 0
    would_replace: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create benchmark reference_data symlinks from the reference inventory."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Benchmark repo root. Defaults to the current directory.",
    )
    parser.add_argument(
        "--inventory-dir",
        type=Path,
        default=DEFAULT_INVENTORY_DIR,
        help="Directory containing by-dataset reference inventory payloads.",
    )
    parser.add_argument(
        "--dataset",
        action="append",
        help="Dataset to process. May be given multiple times. Defaults to all datasets.",
    )
    parser.add_argument(
        "--regime",
        choices=sorted(REGIME_DIRS),
        action="append",
        help="Inventory regime to process. May be given multiple times. Defaults to all regimes.",
    )

    subparsers = parser.add_subparsers(dest="command")
    link_parser = subparsers.add_parser("link", help="Create symlinks from an existing inventory.")
    add_link_options(link_parser)

    install_parser = subparsers.add_parser(
        "install",
        help="Unpack local reference tarballs if needed, then create symlinks.",
    )
    add_link_options(install_parser)
    install_parser.add_argument(
        "--force-unpack",
        action="store_true",
        help="Re-extract local <dataset>.ref_data.tar.gz payloads even when dataset directories exist.",
    )

    args = parser.parse_args()
    if args.command is None:
        args.command = "link"
        args.dry_run = False
        args.force = False
        args.create_missing_dirs = False
        args.skip_raw_proxy_quants = False
    return args


def add_link_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report planned symlink changes without modifying the filesystem.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace existing non-matching files or symlinks at destination paths.",
    )
    parser.add_argument(
        "--create-missing-dirs",
        action="store_true",
        help="Create missing reference_data/raw_prog_results directories instead of reporting them as skipped.",
    )
    parser.add_argument(
        "--skip-raw-proxy-quants",
        action="store_true",
        help="Do not link proxy quantification files into sample raw_prog_results directories.",
    )


def selected_dataset_dirs(inventory_dir: Path, dataset_names: set[str] | None) -> list[Path]:
    if dataset_names:
        return [inventory_dir / dataset for dataset in sorted(dataset_names)]

    return sorted(
        path
        for path in inventory_dir.iterdir()
        if path.is_dir() and (path / "reference_inputs.tsv").exists()
    )


def load_actions(
    inventory_dir: Path,
    repo_root: Path,
    dataset_names: set[str] | None,
    regimes: set[str] | None,
    skip_raw_proxy_quants: bool,
) -> list[LinkAction]:
    actions: set[LinkAction] = set()
    for dataset_dir in selected_dataset_dirs(inventory_dir, dataset_names):
        tsv_path = dataset_dir / "reference_inputs.tsv"
        if not tsv_path.exists():
            raise FileNotFoundError(f"Missing inventory TSV: {tsv_path}")

        with tsv_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                regime = row["regime"]
                if regimes and regime not in regimes:
                    continue
                if regime not in REGIME_DIRS:
                    raise ValueError(f"Unsupported regime {regime!r} in {tsv_path}")
                if row.get("file_type") != "file":
                    continue

                rel_path = clean_relative_path(row["path"], tsv_path)
                source = dataset_dir / "reference_files" / rel_path
                target_dir = target_reference_dir(repo_root, row)
                dest = target_dir / rel_path.name
                actions.add(
                    LinkAction(
                        source=source,
                        dest=dest,
                        dataset=row["target_dataset"],
                        regime=regime,
                        role=row["role"],
                    )
                )
                if not skip_raw_proxy_quants and is_proxy_quant(row):
                    actions.add(
                        LinkAction(
                            source=source,
                            dest=target_raw_prog_results_dir(repo_root, row) / rel_path.name,
                            dataset=row["target_dataset"],
                            regime=regime,
                            role=f"{row['role']}:raw_proxy_quant",
                        )
                    )

    return sorted(actions)


def clean_relative_path(path_value: str, tsv_path: Path) -> Path:
    rel_path = Path(path_value)
    if rel_path.is_absolute() or ".." in rel_path.parts:
        raise ValueError(f"Unsafe inventory path {path_value!r} in {tsv_path}")
    return rel_path


def target_reference_dir(repo_root: Path, row: dict[str, str]) -> Path:
    regime_dir = REGIME_DIRS[row["regime"]]
    dataset = row["target_dataset"]

    if dataset == "MORFs":
        return repo_root / regime_dir / dataset / row["sample"] / "reference_data"
    return repo_root / regime_dir / dataset / "reference_data"


def target_raw_prog_results_dir(repo_root: Path, row: dict[str, str]) -> Path:
    regime_dir = REGIME_DIRS[row["regime"]]
    dataset = row["target_dataset"]
    base_dir = repo_root / regime_dir / dataset

    if dataset == "MORFs":
        sample_dir = base_dir / row["sample"]
    elif row.get("subdataset"):
        sample_dir = base_dir / row["subdataset"] / row["sample"]
    else:
        sample_dir = base_dir / row["sample"]

    return sample_dir / "raw_prog_results"


def is_proxy_quant(row: dict[str, str]) -> bool:
    return (
        row.get("reference_category") == "proxy-quant"
        and row.get("input_type") == "quantification"
        and row.get("file_type") == "file"
    )


def ensure_reference_payloads(args: argparse.Namespace) -> None:
    """Ensure local inventory directories exist before linking.

    This is the intended future hook for downloading Zenodo tarballs. For now it
    supports unpacking tarballs already present beside the dataset directories.
    """

    dataset_names = set(args.dataset) if args.dataset else None
    dataset_dirs = selected_dataset_dirs(args.inventory_dir, dataset_names)
    if dataset_names:
        missing = sorted(dataset_names - {path.name for path in dataset_dirs if path.exists()})
        dataset_dirs.extend(args.inventory_dir / dataset for dataset in missing)

    for dataset_dir in dataset_dirs:
        reference_files_dir = dataset_dir / "reference_files"
        if reference_files_dir.exists() and not args.force_unpack:
            continue

        tarball = args.inventory_dir / f"{dataset_dir.name}.ref_data.tar.gz"
        if not tarball.exists():
            if not dataset_dir.exists():
                raise FileNotFoundError(
                    f"Missing dataset inventory and local tarball for {dataset_dir.name}: {tarball}"
                )
            continue

        print(f"Unpacking {tarball} into {args.inventory_dir}")
        safe_extract_tarball(tarball, args.inventory_dir)


def safe_extract_tarball(tarball: Path, dest_dir: Path) -> None:
    dest_dir = dest_dir.resolve()
    with tarfile.open(tarball, "r:*") as tar:
        for member in tar.getmembers():
            member_path = (dest_dir / member.name).resolve()
            if os.path.commonpath([dest_dir, member_path]) != str(dest_dir):
                raise ValueError(f"Refusing unsafe tar member path: {member.name}")
        tar.extractall(dest_dir)


def apply_actions(actions: Iterable[LinkAction], args: argparse.Namespace) -> LinkStats:
    stats = LinkStats()
    for action in actions:
        stats = apply_action(action, args, stats)
    return stats


def apply_action(action: LinkAction, args: argparse.Namespace, stats: LinkStats) -> LinkStats:
    if not action.source.exists():
        stats.missing_source += 1
        print(f"MISSING source: {action.source}", file=sys.stderr)
        return stats

    if not action.dest.parent.exists():
        if args.create_missing_dirs:
            if args.dry_run:
                print(f"Would create directory: {action.dest.parent}")
            else:
                action.dest.parent.mkdir(parents=True, exist_ok=True)
        else:
            stats.missing_target_dir += 1
            print(f"MISSING target directory: {action.dest.parent}", file=sys.stderr)
            return stats

    if action.dest.is_symlink() and same_symlink_target(action.dest, action.source):
        stats.unchanged += 1
        return stats

    if path_lexists(action.dest):
        if not args.force:
            stats.conflicts += 1
            print(f"CONFLICT destination exists: {action.dest}", file=sys.stderr)
            return stats

        if action.dest.is_dir() and not action.dest.is_symlink():
            stats.conflicts += 1
            print(f"CONFLICT destination is a directory: {action.dest}", file=sys.stderr)
            return stats

        if args.dry_run:
            stats.would_replace += 1
            print(f"Would replace: {action.dest} -> {action.source}")
            return stats

        action.dest.unlink()
        action.dest.symlink_to(action.source)
        stats.replaced += 1
        print(f"Replaced: {action.dest} -> {action.source}")
        return stats

    if args.dry_run:
        stats.would_create += 1
        print(f"Would link: {action.dest} -> {action.source}")
        return stats

    action.dest.symlink_to(action.source)
    stats.created += 1
    print(f"Linked: {action.dest} -> {action.source}")
    return stats


def path_lexists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def same_symlink_target(link_path: Path, source: Path) -> bool:
    try:
        return link_path.resolve(strict=True) == source.resolve(strict=True)
    except FileNotFoundError:
        return False


def print_summary(stats: LinkStats) -> None:
    print(
        "Summary: "
        f"created={stats.created} "
        f"unchanged={stats.unchanged} "
        f"replaced={stats.replaced} "
        f"would_create={stats.would_create} "
        f"would_replace={stats.would_replace} "
        f"missing_source={stats.missing_source} "
        f"missing_target_dir={stats.missing_target_dir} "
        f"conflicts={stats.conflicts}"
    )


def main() -> int:
    args = parse_args()
    args.repo_root = args.repo_root.resolve()
    args.inventory_dir = args.inventory_dir.resolve()
    dataset_names = set(args.dataset) if args.dataset else None
    regimes = set(args.regime) if args.regime else None

    if args.command == "install":
        ensure_reference_payloads(args)

    actions = load_actions(
        args.inventory_dir,
        args.repo_root,
        dataset_names,
        regimes,
        args.skip_raw_proxy_quants,
    )
    stats = apply_actions(actions, args)
    print_summary(stats)

    if stats.missing_source or stats.missing_target_dir or stats.conflicts:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
