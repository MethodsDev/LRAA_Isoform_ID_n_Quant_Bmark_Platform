#!/usr/bin/env python3

import argparse
import glob
import os
import re
import shutil
import sys


DEBUG = False


def parse_args():
    parser = argparse.ArgumentParser(
        description="Copy comparator result files into their matching raw_prog_results directories."
    )
    parser.add_argument(
        "program_dirs",
        nargs="*",
        default=["*"],
        help=(
            "Program output directories to process (e.g. Flair IsoQuant), "
            "or '*' for all directories under --input-root."
        ),
    )
    parser.add_argument(
        "--input-root",
        default=".",
        help="Root dir containing program output directories. Default: current dir.",
    )
    parser.add_argument(
        "--base-project-out-dir",
        default="../",
        help=(
            "Base project output dir containing CellLines/SIRVs/... "
            "Default: ../ relative to current dir."
        ),
    )
    parser.add_argument(
        "--registry",
        default=None,
        help=(
            "Path to tool_registry.yaml. Used to warn loudly on "
            "deposited filenames that no registry entry recognizes "
            "(typically a brand-new tool/version that needs an entry "
            "added). If omitted, walks up from --base-project-out-dir "
            "looking for tool_registry.yaml; deposit proceeds with no "
            "warnings if none is found."
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Dry-run mode; report copies without writing files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    global DEBUG
    DEBUG = args.debug

    input_root = os.path.abspath(args.input_root)
    base_project_out_dir = os.path.abspath(args.base_project_out_dir)
    program_dirs = resolve_program_dirs(input_root, args.program_dirs)
    registry_entries = _try_load_registry(args.registry, base_project_out_dir)

    input_files = []
    for program_dir in program_dirs:
        input_files.extend(
            f
            for f in glob.glob(os.path.join(program_dir, "*"))
            if os.path.isfile(f)
        )
    input_files = sorted(input_files)

    if not input_files:
        print("No files found to process.")
        return 0

    unrecognized_by_registry = []

    for filename in input_files:
        base_fname = os.path.basename(filename)

        if base_fname in {"gs.files.list", "notes"}:
            continue

        if (
            deposit_cellline_or_sirv(filename, base_project_out_dir)
            or deposit_isoseqsim(filename, base_project_out_dir)
            or deposit_morf2(filename, base_project_out_dir)
            or deposit_sgnex(filename, base_project_out_dir)
        ):
            if registry_entries is not None and not _matches_any_entry(base_fname, registry_entries):
                unrecognized_by_registry.append(base_fname)
            continue

        print(f"************ Error, not able to process file: {filename}")

    if unrecognized_by_registry:
        print(
            "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
            f"WARNING: {len(unrecognized_by_registry)} deposited file(s) "
            "matched no entry in the tool registry. They will be ignored "
            "by bmark_nb_runner.py until you add a matching entry to "
            "tool_registry.yaml:"
        )
        for f in unrecognized_by_registry:
            print(f"  - {f}")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

    print("\n\nDone.\n")
    return 0


def _try_load_registry(explicit_path, base_project_out_dir):
    """Load the tool registry if findable; return None if not.

    The registry is OPTIONAL for deposit -- the script still runs without
    one. It's only used to surface 'this filename has no parser entry'
    warnings.
    """
    if explicit_path:
        path = explicit_path
    else:
        path = _find_registry_walking_up(base_project_out_dir)
        if path is None:
            return None

    # Lazy-import yaml + QuantParser so missing PyYAML doesn't break
    # the deposit script for callers that don't pass --registry.
    try:
        import yaml  # noqa: F401
    except ImportError:
        print(
            f"NOTE: PyYAML not installed; skipping registry-aware warnings. "
            f"Install pyyaml to enable them."
        )
        return None

    iso_pylib = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "..", "..", "benchmarking", "pylib",
        )
    )
    if iso_pylib not in sys.path:
        sys.path.insert(0, iso_pylib)
    try:
        import QuantParser
    except Exception as e:
        print(
            "NOTE: Could not import QuantParser; skipping registry-aware "
            f"warnings. Details: {e}"
        )
        return None

    try:
        return QuantParser.load_registry(path)
    except Exception as e:
        print(
            f"NOTE: Could not load registry from {path}; skipping "
            f"registry-aware warnings. Details: {e}"
        )
        return None


def _find_registry_walking_up(start_dir):
    d = os.path.abspath(start_dir)
    while True:
        candidate = os.path.join(d, "tool_registry.yaml")
        if os.path.isfile(candidate):
            return candidate
        parent = os.path.dirname(d)
        if parent == d:
            return None
        d = parent


def _matches_any_entry(basename, entries):
    for e in entries:
        if e.get("quant_pattern") and re.search(e["quant_pattern"], basename):
            return True
        if e.get("gtf_pattern") and re.search(e["gtf_pattern"], basename):
            return True
    return False


def resolve_program_dirs(input_root, program_dir_patterns):
    found_dirs = []
    seen = set()

    for pattern in program_dir_patterns:
        if os.path.isabs(pattern):
            glob_pattern = pattern
        else:
            glob_pattern = os.path.join(input_root, pattern)

        matches = sorted(glob.glob(glob_pattern))
        for match in matches:
            if os.path.isdir(match):
                abspath = os.path.abspath(match)
                if abspath not in seen:
                    seen.add(abspath)
                    found_dirs.append(abspath)

    if not found_dirs:
        raise RuntimeError(
            "No matching program directories found for: {}".format(
                ", ".join(program_dir_patterns)
            )
        )

    return found_dirs


def deposit_morf2(filename, proj_base_dir):
    base_fname = os.path.basename(filename)

    m = re.match(r"(morf2_[^\.]+)", base_fname)
    if m is None:
        return False

    morf2_sample_name = m.group(1)
    dest_loc = os.path.join(
        proj_base_dir, "MORFs", morf2_sample_name, "raw_prog_results"
    )

    assert_dest_exists(dest_loc)
    copy_file(filename, dest_loc)
    return True


def deposit_isoseqsim(filename, proj_base_dir):
    base_fname = os.path.basename(filename)

    m = re.search(r"^([^_]+)_isoseqsim_(e\d{3})\.", base_fname)
    if m is None:
        return False

    orgname = m.group(1)
    error_rate = m.group(2)

    dest_loc = os.path.join(
        proj_base_dir,
        f"{orgname}_sim",
        f"{orgname}_isoseqsim_{error_rate}",
        "raw_prog_results",
    )

    assert_dest_exists(dest_loc)
    copy_file(filename, dest_loc)
    return True


def deposit_cellline_or_sirv(filename, proj_base_dir):
    base_fname = os.path.basename(filename)

    m = re.search(r"CL_([^_]+)_E(\d)_([^\.]+)\.", base_fname)
    if m is None:
        return False

    cell_line = m.group(1)
    erep = m.group(2)
    sirv_or_human = m.group(3)

    if sirv_or_human == "human":
        sample_name = f"CL_{cell_line}_E{erep}_human"
        dest_loc = os.path.join(
            proj_base_dir,
            "CellLines",
            sample_name,
            "raw_prog_results",
        )
    elif sirv_or_human == "sirv":
        sample_name = f"CL_{cell_line}_E{erep}_sirv"
        dest_loc = os.path.join(
            proj_base_dir,
            "SIRVs",
            f"e{erep}_sirv",
            sample_name,
            "raw_prog_results",
        )
    else:
        raise RuntimeError(f"Error, not recognizing {filename} as human or sirv")

    assert_dest_exists(dest_loc)
    copy_file(filename, dest_loc)
    return True


def deposit_sgnex(filename, proj_base_dir):
    base_fname = os.path.basename(filename)

    m = re.search(r"^(SGNex_[^_]+)_cDNAStranded", base_fname)
    if m is None:
        return False

    sgnex_sample_name = m.group(1)

    dest_loc = os.path.join(
        proj_base_dir,
        "SG-NEx",
        sgnex_sample_name,
        "raw_prog_results",
    )

    assert_dest_exists(dest_loc)
    copy_file(filename, dest_loc)
    return True


def assert_dest_exists(dest_loc):
    if not os.path.exists(dest_loc):
        raise RuntimeError(f"Error, destination path does not exist: {dest_loc}")


def copy_file(target, dest_dir):
    out_file = os.path.join(dest_dir, os.path.basename(target))
    if DEBUG:
        print(f"-[DEBUGMODE] would copy {target} -> {out_file}")
    else:
        print(f"-copying {target} -> {out_file}")
        shutil.copyfile(target, out_file)


if __name__ == "__main__":
    sys.exit(main())
