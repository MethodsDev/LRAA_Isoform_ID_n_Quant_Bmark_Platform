#!/usr/bin/env python3
"""Download deposited benchmark predictions from Zenodo and link them locally."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path

from install_references import (
    DEFAULT_ZENODO_RECORD,
    RemoteFile,
    checksum_matches,
    create_relative_symlink,
    download_remote_file,
    fetch_zenodo_files,
    path_lexists,
    safe_extract_tarball,
)


PREDICTION_ARCHIVE_SUFFIX = ".predictions.tar.gz"
DEFAULT_CACHE_DIRNAME = ".benchmark_inputs/predictions"
ARCHIVE_ROOTS = {
    "DENOVO": "DENOVO_ID",
    "QUANT_ONLY": "QUANT_ONLY",
    "REF_GUIDED": "REF_Guided",
}


@dataclass
class LinkStats:
    created: int = 0
    unchanged: int = 0
    replaced: int = 0
    conflicts: int = 0
    would_create: int = 0
    would_replace: int = 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download deposited benchmark predictions from Zenodo."
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path.cwd(),
        help="Benchmark repository root. Defaults to the current directory.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help=f"Download cache. Defaults to <repo-root>/{DEFAULT_CACHE_DIRNAME}.",
    )
    parser.add_argument(
        "--zenodo-record",
        default=DEFAULT_ZENODO_RECORD,
        help=f"Zenodo record ID (default: {DEFAULT_ZENODO_RECORD}).",
    )
    parser.add_argument(
        "--regime",
        action="append",
        choices=sorted(ARCHIVE_ROOTS),
        help="Regime to install. May be repeated; defaults to all regimes.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Download and unpack inputs, but only report link changes.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace conflicting files or links in raw_prog_results.",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Download selected archives again even when cached copies are valid.",
    )
    parser.add_argument(
        "--force-unpack",
        action="store_true",
        help="Re-extract selected archives even when an unpacked cache exists.",
    )
    return parser.parse_args()


def unpack_prediction_archive(
    remote_file: RemoteFile,
    archive_key: str,
    archive_root: str,
    cache_dir: Path,
    force_download: bool,
    force_unpack: bool,
) -> Path:
    archive = cache_dir / remote_file.name
    if force_download or not checksum_matches(archive, remote_file.checksum, remote_file.size):
        download_remote_file(remote_file, archive, force_download)
    else:
        print(f"Using verified archive: {archive}")

    digest = remote_file.checksum.replace(":", "-")
    unpacked_dir = cache_dir / "unpacked" / f"{archive_key}.{digest}"
    marker = unpacked_dir / ".complete"
    if marker.is_file() and not force_unpack:
        print(f"Using unpacked predictions: {unpacked_dir}")
        return unpacked_dir

    unpacked_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=f".{archive_key}.", dir=unpacked_dir.parent)
    )
    try:
        print(f"Unpacking {archive} into prediction cache")
        safe_extract_tarball(archive, temporary_dir)
        validate_prediction_tree(temporary_dir, archive_root, remote_file.name)
        (temporary_dir / ".complete").write_text(f"{remote_file.checksum}\n")
        if unpacked_dir.exists():
            shutil.rmtree(unpacked_dir)
        os.replace(temporary_dir, unpacked_dir)
    except Exception:
        shutil.rmtree(temporary_dir, ignore_errors=True)
        raise
    return unpacked_dir


def validate_prediction_tree(unpacked_dir: Path, archive_root: str, archive_name: str) -> None:
    root = unpacked_dir / archive_root
    if not root.is_dir():
        raise ValueError(f"{archive_name} does not contain the expected {archive_root}/ root")

    prediction_files = list(iter_prediction_files(unpacked_dir, archive_root))
    if not prediction_files:
        raise ValueError(f"{archive_name} contains no raw_prog_results files")

    unexpected = [
        path
        for path in root.rglob("*")
        if path.is_file() and "raw_prog_results" not in path.relative_to(root).parts
    ]
    if unexpected:
        raise ValueError(
            f"Unexpected file outside raw_prog_results in {archive_name}: {unexpected[0]}"
        )


def iter_prediction_files(unpacked_dir: Path, archive_root: str):
    root = unpacked_dir / archive_root
    for raw_results_dir in sorted(root.rglob("raw_prog_results")):
        if raw_results_dir.is_dir():
            yield from sorted(path for path in raw_results_dir.rglob("*") if path.is_file())


def link_predictions(
    unpacked_dir: Path,
    archive_root: str,
    repo_root: Path,
    force: bool,
    dry_run: bool,
    stats: LinkStats,
) -> None:
    for source in iter_prediction_files(unpacked_dir, archive_root):
        destination = repo_root / source.relative_to(unpacked_dir)
        destination.parent.mkdir(parents=True, exist_ok=True)

        if destination.is_symlink() and destination.resolve() == source.resolve():
            stats.unchanged += 1
            continue
        if path_lexists(destination):
            if destination.is_file() and files_equal(destination, source):
                stats.unchanged += 1
                continue
            if not force or (destination.is_dir() and not destination.is_symlink()):
                stats.conflicts += 1
                print(f"CONFLICT destination exists: {destination}", file=sys.stderr)
                continue
            if dry_run:
                stats.would_replace += 1
                print(f"Would replace: {destination} -> {source}")
                continue
            destination.unlink()
            create_relative_symlink(destination, source)
            stats.replaced += 1
            print(f"Replaced: {destination} -> {source}")
            continue

        if dry_run:
            stats.would_create += 1
            print(f"Would link: {destination} -> {source}")
        else:
            create_relative_symlink(destination, source)
            stats.created += 1
            print(f"Linked: {destination} -> {source}")


def files_equal(left: Path, right: Path) -> bool:
    if left.stat().st_size != right.stat().st_size:
        return False
    with left.open("rb") as left_handle, right.open("rb") as right_handle:
        while True:
            left_chunk = left_handle.read(1024 * 1024)
            right_chunk = right_handle.read(1024 * 1024)
            if left_chunk != right_chunk:
                return False
            if not left_chunk:
                return True


def main() -> int:
    args = parse_args()
    args.repo_root = args.repo_root.resolve()
    args.cache_dir = (
        args.cache_dir.resolve()
        if args.cache_dir is not None
        else args.repo_root / DEFAULT_CACHE_DIRNAME
    )
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    try:
        remote_files = fetch_zenodo_files(args.zenodo_record, PREDICTION_ARCHIVE_SUFFIX)
        regimes = set(args.regime) if args.regime else set(ARCHIVE_ROOTS)
        missing = regimes - set(remote_files)
        if missing:
            raise ValueError(
                f"Prediction archive(s) missing from Zenodo: {', '.join(sorted(missing))}"
            )

        stats = LinkStats()
        for regime in sorted(regimes):
            unpacked_dir = unpack_prediction_archive(
                remote_files[regime],
                regime,
                ARCHIVE_ROOTS[regime],
                args.cache_dir,
                args.force_download,
                args.force_unpack,
            )
            link_predictions(
                unpacked_dir,
                ARCHIVE_ROOTS[regime],
                args.repo_root,
                args.force,
                args.dry_run,
                stats,
            )
    except (OSError, RuntimeError, ValueError, tarfile.TarError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        "Summary: "
        f"created={stats.created} unchanged={stats.unchanged} replaced={stats.replaced} "
        f"would_create={stats.would_create} would_replace={stats.would_replace} "
        f"conflicts={stats.conflicts}"
    )
    return 1 if stats.conflicts else 0


if __name__ == "__main__":
    raise SystemExit(main())
