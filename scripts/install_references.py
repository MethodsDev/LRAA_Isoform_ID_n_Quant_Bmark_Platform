#!/usr/bin/env python3
"""Install benchmark reference-data symlinks from Zenodo.

Reference archives are downloaded from the pinned Zenodo record into a local
cache, unpacked, and linked into the benchmark tree.  The unpacked inventory
has this layout:

    by_dataset/<dataset>/reference_inputs.tsv
    by_dataset/<dataset>/reference_files/<file>

The ``link`` command remains available for an already-populated inventory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
import tarfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_ZENODO_RECORD = "20719930"
DEFAULT_CACHE_DIRNAME = ".benchmark_inputs/reference_inventory/by_dataset"
ZENODO_API_ROOT = "https://zenodo.org/api/records"
REFERENCE_ARCHIVE_SUFFIX = ".ref_data.tar.gz"

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


@dataclass(frozen=True)
class RemoteFile:
    name: str
    url: str
    size: int
    checksum: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download benchmark references from Zenodo and link them into the repository."
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
        help=(
            "Cache containing by-dataset reference inventory payloads. "
            f"Defaults to <repo-root>/{DEFAULT_CACHE_DIRNAME}."
        ),
    )
    parser.add_argument(
        "--zenodo-record",
        default=DEFAULT_ZENODO_RECORD,
        help=f"Zenodo record ID used by install (default: {DEFAULT_ZENODO_RECORD}).",
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
        help="Download reference archives from Zenodo, unpack them, and create symlinks.",
    )
    add_link_options(install_parser)
    install_parser.add_argument(
        "--force-unpack",
        action="store_true",
        help=(
            "Re-extract cached <dataset>.ref_data.tar.gz archives even when "
            "dataset directories exist."
        ),
    )
    install_parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download selected archives again even when valid cached copies exist.",
    )

    args = parser.parse_args()
    if args.command is None:
        args.command = "install"
        args.dry_run = False
        args.force = False
        args.create_missing_dirs = False
        args.skip_raw_proxy_quants = False
        args.force_unpack = False
        args.force_download = False
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
        help=(
            "Create missing reference_data/raw_prog_results directories instead "
            "of reporting them as skipped."
        ),
    )
    parser.add_argument(
        "--skip-raw-proxy-quants",
        action="store_true",
        help="Do not link proxy quantification files into sample raw_prog_results directories.",
    )


def selected_dataset_dirs(inventory_dir: Path, dataset_names: set[str] | None) -> list[Path]:
    if dataset_names:
        return [inventory_dir / dataset for dataset in sorted(dataset_names)]

    if not inventory_dir.is_dir():
        raise FileNotFoundError(
            f"Reference inventory does not exist: {inventory_dir}. Run the install command first."
        )

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
    """Download and unpack the selected Zenodo reference archives."""

    args.inventory_dir.mkdir(parents=True, exist_ok=True)
    remote_files = fetch_zenodo_reference_files(args.zenodo_record)
    available_datasets = set(remote_files)
    dataset_names = set(args.dataset) if args.dataset else available_datasets
    unknown_datasets = dataset_names - available_datasets
    if unknown_datasets:
        unknown = ", ".join(sorted(unknown_datasets))
        available = ", ".join(sorted(available_datasets))
        raise ValueError(
            f"Dataset(s) not found in Zenodo record: {unknown}. Available: {available}"
        )

    for dataset_name in sorted(dataset_names):
        dataset_dir = args.inventory_dir / dataset_name
        reference_files_dir = dataset_dir / "reference_files"
        if (
            reference_files_dir.is_dir()
            and (dataset_dir / "reference_inputs.tsv").is_file()
            and not args.force_unpack
            and not args.force_download
        ):
            print(f"Using unpacked reference data: {dataset_dir}")
            continue

        remote_file = remote_files[dataset_name]
        tarball = args.inventory_dir / remote_file.name
        if args.force_download or not checksum_matches(
            tarball, remote_file.checksum, remote_file.size
        ):
            download_remote_file(remote_file, tarball, args.force_download)
        else:
            print(f"Using verified archive: {tarball}")

        print(f"Unpacking {tarball} into {args.inventory_dir}")
        safe_extract_tarball(tarball, args.inventory_dir)
        if not reference_files_dir.is_dir() or not (dataset_dir / "reference_inputs.tsv").is_file():
            raise ValueError(
                f"Archive {tarball.name} did not contain the expected {dataset_name} inventory"
            )


def fetch_zenodo_reference_files(record_id: str) -> dict[str, RemoteFile]:
    """Return the per-dataset reference archives published in a Zenodo record."""

    return fetch_zenodo_files(record_id, REFERENCE_ARCHIVE_SUFFIX)


def fetch_zenodo_files(record_id: str, suffix: str) -> dict[str, RemoteFile]:
    """Return files with a given suffix from a published Zenodo record."""

    if not record_id.isdigit():
        raise ValueError(f"Zenodo record ID must be numeric: {record_id!r}")

    url = f"{ZENODO_API_ROOT}/{record_id}"
    print(f"Reading Zenodo record metadata: {url}")
    request = urllib.request.Request(url, headers={"User-Agent": "LRAA-benchmark-installer/1"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            record = json.load(response)
    except (urllib.error.URLError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read Zenodo record {record_id}: {exc}") from exc

    remote_files: dict[str, RemoteFile] = {}
    for entry in record.get("files", []):
        name = entry.get("key", "")
        if not name.endswith(suffix):
            continue
        key = name[: -len(suffix)]
        links = entry.get("links", {})
        content_url = links.get("content") or links.get("self")
        checksum = entry.get("checksum")
        size = entry.get("size")
        if not content_url or not checksum or not isinstance(size, int):
            raise ValueError(f"Zenodo metadata is incomplete for {name}")
        remote_files[key] = RemoteFile(name, content_url, size, checksum)

    if not remote_files:
        raise ValueError(f"Zenodo record {record_id} contains no *{suffix} files")
    return remote_files


def checksum_matches(path: Path, checksum: str, expected_size: int) -> bool:
    if not path.is_file() or path.stat().st_size != expected_size:
        return False

    algorithm, separator, expected_digest = checksum.partition(":")
    if not separator:
        raise ValueError(f"Invalid checksum from Zenodo: {checksum!r}")
    try:
        digest = hashlib.new(algorithm)
    except ValueError as exc:
        raise ValueError(f"Unsupported checksum algorithm from Zenodo: {algorithm}") from exc

    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().lower() == expected_digest.lower()


def download_remote_file(remote_file: RemoteFile, destination: Path, restart: bool) -> None:
    """Download to a partial file, resuming it when Zenodo supports byte ranges."""

    partial = destination.with_name(f"{destination.name}.part")
    if restart and partial.exists():
        partial.unlink()
    start_at = partial.stat().st_size if partial.exists() else 0
    if start_at > remote_file.size:
        partial.unlink()
        start_at = 0
    elif start_at == remote_file.size:
        if checksum_matches(partial, remote_file.checksum, remote_file.size):
            os.replace(partial, destination)
            print(f"Downloaded and verified: {destination}")
            return
        partial.unlink()
        start_at = 0

    headers = {"User-Agent": "LRAA-benchmark-installer/1"}
    if start_at:
        headers["Range"] = f"bytes={start_at}-"
        print(f"Resuming {remote_file.name} at {format_bytes(start_at)}")
    else:
        print(f"Downloading {remote_file.name} ({format_bytes(remote_file.size)})")

    request = urllib.request.Request(remote_file.url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            response_status = getattr(response, "status", response.getcode())
            append = start_at > 0 and response_status == 206
            if start_at and not append:
                start_at = 0
            mode = "ab" if append else "wb"
            downloaded = start_at
            next_report = downloaded + 256 * 1024 * 1024
            with partial.open(mode) as handle:
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    downloaded += len(chunk)
                    if downloaded >= next_report:
                        print(
                            f"  {remote_file.name}: {format_bytes(downloaded)} "
                            f"of {format_bytes(remote_file.size)}"
                        )
                        next_report = downloaded + 256 * 1024 * 1024
    except urllib.error.URLError as exc:
        raise RuntimeError(
            f"Download failed for {remote_file.name}; rerun to resume it: {exc}"
        ) from exc

    if partial.stat().st_size != remote_file.size:
        raise RuntimeError(
            f"Incomplete download for {remote_file.name}: got {partial.stat().st_size} "
            f"bytes, expected {remote_file.size}; rerun to resume it"
        )
    if not checksum_matches(partial, remote_file.checksum, remote_file.size):
        raise RuntimeError(
            f"Checksum mismatch for {remote_file.name}; remove {partial} and try again"
        )

    os.replace(partial, destination)
    print(f"Downloaded and verified: {destination}")


def format_bytes(size: int) -> str:
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def safe_extract_tarball(tarball: Path, dest_dir: Path) -> None:
    dest_dir = dest_dir.resolve()
    with tarfile.open(tarball, "r:*") as tar:
        for member in tar.getmembers():
            member_path = (dest_dir / member.name).resolve()
            if os.path.commonpath([dest_dir, member_path]) != str(dest_dir):
                raise ValueError(f"Refusing unsafe tar member path: {member.name}")
        if hasattr(tarfile, "data_filter"):
            tar.extractall(dest_dir, filter="data")
        else:
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
        create_relative_symlink(action.dest, action.source)
        stats.replaced += 1
        print(f"Replaced: {action.dest} -> {action.source}")
        return stats

    if args.dry_run:
        stats.would_create += 1
        print(f"Would link: {action.dest} -> {action.source}")
        return stats

    create_relative_symlink(action.dest, action.source)
    stats.created += 1
    print(f"Linked: {action.dest} -> {action.source}")
    return stats


def path_lexists(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def create_relative_symlink(dest: Path, source: Path) -> None:
    dest.symlink_to(os.path.relpath(source, start=dest.parent))


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
    if args.inventory_dir is None:
        args.inventory_dir = args.repo_root / DEFAULT_CACHE_DIRNAME
    args.inventory_dir = args.inventory_dir.resolve()
    dataset_names = set(args.dataset) if args.dataset else None
    regimes = set(args.regime) if args.regime else None

    try:
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
    except (OSError, RuntimeError, ValueError, tarfile.TarError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if stats.missing_source or stats.missing_target_dir or stats.conflicts:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
