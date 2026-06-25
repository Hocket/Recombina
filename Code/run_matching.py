"""
run_matching.py

Main workflow for the ViralRecombinant pipeline.

Examples
--------
Generate PED/INFO only (no pipeline):
    python run_matching.py --alignment alignment.fasta --pedinfo-only

Run full pipeline (with temp intermediate files, auto-cleanup):
    python run_matching.py --alignment alignment.fasta --haploview-jar path/to/haploview.jar

Run full pipeline and keep intermediate files:
    python run_matching.py --alignment alignment.fasta --haploview-jar path/to/haploview.jar --keep-intermediate

Show help:
    python run_matching.py --help

Author: Alex Poyer
"""

import argparse
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
import os
import time
import threading

import matching
import tree_utils
import fasta_to_ped as ftp


# --- Helper utilities -----------------------------------------------------
def predict_total_pairs_from_info(info_file_path: Path) -> int:
    """Estimate upper bound of pairwise comparisons from an INFO file.

    Counts unique positions from the INFO file (column 2) and returns
    n*(n-1)//2 as an upper bound of pairwise comparisons.
    """
    positions = set()
    try:
        with open(info_file_path, "r") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) >= 2:
                    positions.add(parts[1])
    except Exception:
        return 0
    n = len(positions)
    return n * (n - 1) // 2


def monitor_ld_file(ld_path: Path, predicted_total: int, proc: subprocess.Popen, interval: float = 1.0):
    """Tail an LD file while `proc` runs and display progress.

    Uses `tqdm` if available; otherwise prints a single-line status updated
    in-place using carriage return. The function returns after `proc`
    exits and the file stabilizes.
    """
    try:
        from tqdm import tqdm
        has_tqdm = True
    except Exception:
        has_tqdm = False

    last_size = 0
    last_count = 0

    if has_tqdm:
        pbar = tqdm(total=predicted_total if predicted_total > 0 else None, unit="pairs", desc="Haploview LD")
    else:
        print("Monitoring LD file (press Ctrl-C to stop)...")

    # loop until process ends and file stabilizes
    while True:
        if ld_path.exists():
            try:
                size = os.path.getsize(ld_path)
            except OSError:
                size = last_size

            if size > last_size:
                # read only the new bytes and count newlines
                try:
                    with open(ld_path, "rb") as fh:
                        fh.seek(last_size)
                        data = fh.read(size - last_size)
                except Exception:
                    data = b""
                new_lines = data.count(b"\n")
                last_size = size
            else:
                new_lines = 0

            last_count += new_lines
            if has_tqdm:
                if new_lines:
                    pbar.update(new_lines)
            else:
                pct = (last_count / predicted_total * 100) if predicted_total else 0
                print(f"\rLD lines: {last_count} / {predicted_total} ({pct:.1f}%)", end="", flush=True)

        # If process finished and file hasn't grown for a short time, exit
        if proc.poll() is not None:
            # allow final writes
            time.sleep(1.5)
            if ld_path.exists():
                # read any remaining bytes
                try:
                    size = os.path.getsize(ld_path)
                    if size > last_size:
                        with open(ld_path, "rb") as fh:
                            fh.seek(last_size)
                            data = fh.read()
                        last_count += data.count(b"\n")
                except Exception:
                    pass
            break

        time.sleep(interval)

    if has_tqdm:
        pbar.close()
    else:
        print()


def make_run_dir(base_dir: Path, stem: str) -> Path:
    """Create a unique run directory under base_dir using stem.

    If `base_dir/stem` exists, append an incrementing suffix `_1`, `_2`, ...
    to avoid clobbering previous runs.
    """
    base_dir.mkdir(parents=True, exist_ok=True)
    run_dir = base_dir / stem
    if not run_dir.exists():
        run_dir.mkdir(parents=True, exist_ok=False)
        return run_dir
    idx = 1
    while True:
        candidate = base_dir / f"{stem}_{idx}"
        if not candidate.exists():
            candidate.mkdir(parents=True, exist_ok=False)
            return candidate
        idx += 1


ROOT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = ROOT_DIR / "Data"
INPUT_DIR = DATA_DIR / "InputFiles"
OUTPUT_DIR = DATA_DIR / "OutputFiles"
IQTREE_OUT_DIR = DATA_DIR / "IQTree_out"

OUTPUT_PREFIX = "run_matching"
CONTREE = IQTREE_OUT_DIR / f"{OUTPUT_PREFIX}.contree"


def main(
    alignment_file: str | Path,
    haploview_jar: str | Path = None,
    output_excel: str | Path = OUTPUT_DIR / "matches_output.xlsx",
    keep_intermediate: bool = False,
    pedinfo_only: bool = False,
    skip_iqtree: bool = False,
) -> None:
    """
    Full pipeline:
        1. (Optional) Generate PED + INFO from FASTA if not supplied.
        2. Run IQ-TREE to build the phylogenetic tree.
        3. Match recombinant SNP pairs per sample.
        4. Summarize and export to Excel + colored Nexus tree.
    """

    alignment_file = Path(alignment_file)

    # --- Validate alignment ---------------------------------------------------
    if not alignment_file.exists():
        print(f"ERROR: alignment file not found: {alignment_file}", file=sys.stderr)
        sys.exit(1)

    # --- Step 1: PED/INFO generation ------------------------------------------
    print("\n[Step 1] Generating PED and INFO files from FASTA alignment...")
    # Create a per-run output directory to contain all outputs for this alignment
    run_dir = make_run_dir(OUTPUT_DIR, alignment_file.stem)
    run_iqtree_dir = run_dir / "IQTree_out"
    run_iqtree_dir.mkdir(parents=True, exist_ok=True)

    # Place PED/INFO/LD in the run directory so users inspect live outputs there.
    if keep_intermediate:
        ped_path = run_dir / f"{alignment_file.stem}_ped.txt"
        info_path = run_dir / f"{alignment_file.stem}_info.txt"
        recombinant_path = run_dir / f"{alignment_file.stem}_haploview_ld.txt"
    else:
        tempdir = tempfile.TemporaryDirectory()
        ped_path = Path(tempdir.name) / f"{alignment_file.stem}_ped.txt"
        info_path = Path(tempdir.name) / f"{alignment_file.stem}_info.txt"
        recombinant_path = (
            Path(tempdir.name) / f"{alignment_file.stem}_haploview_ld.txt"
        )

    ftp.fasta_to_ped(
        alignment_file,
        ped_path,
        info_path,
        verbose=True,
    )

    print(f"\tPED  → {ped_path}")
    print(f"\tINFO → {info_path}")

    if pedinfo_only:
        print(
            "PED and INFO file generation complete. Exiting as requested by --pedinfo-only flag."
        )
        if not keep_intermediate and "tempdir" in locals():
            tempdir.cleanup()
        return

    # --- Step 2: Run Haploview to generate recombinant file -------------------
    print("\n[Step 2] Running Haploview to generate recombinant file...")

    # Haploview files placed in run_dir so the LD file Haploview writes
    # is created directly in the run-specific output folder (live-updating file).
    haploview_ped = run_dir / ped_path.name
    haploview_info = run_dir / info_path.name

    # If the ped/info are not already in the run-specific directory, copy them
    # there so Haploview will write its LD output into the run folder (live file).
    copied_to_output = False
    if ped_path.parent != run_dir:
        # Only copy when source and destination are different to avoid SameFileError
        if ped_path.resolve() != haploview_ped.resolve():
            shutil.copy(str(ped_path), str(haploview_ped))
        if info_path.resolve() != haploview_info.resolve():
            shutil.copy(str(info_path), str(haploview_info))
        copied_to_output = True
    else:
        haploview_ped = ped_path
        haploview_info = info_path

    haploview_cmd = [
        "java",
        "-jar",
        str(haploview_jar),
        "-nogui",
        "-pedfile",
        str(haploview_ped),
        "-info",
        str(haploview_info),
        "-dprime",
    ]

    print("Running:", " ".join(haploview_cmd))

    # Start Haploview as subprocess and monitor the LD file it writes
    proc = subprocess.Popen(haploview_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Predict total pairs (upper bound) from the INFO file located in run_dir
    predicted_total = predict_total_pairs_from_info(haploview_info)

    # Haploview will write its LD file using the ped filename in OUTPUT_DIR
    haploview_ld = Path(str(haploview_ped) + ".LD")
    if not haploview_ld.exists():
        haploview_ld = Path(str(haploview_ped) + ".ld")

    monitor_thread = threading.Thread(target=monitor_ld_file, args=(haploview_ld, predicted_total, proc), daemon=True)
    monitor_thread.start()

    stdout, stderr = proc.communicate()
    monitor_thread.join()

    if proc.returncode != 0:
        print("ERROR: Haploview failed:", file=sys.stderr)
        print(stderr, file=sys.stderr)
        if not keep_intermediate and "tempdir" in locals():
            try:
                tempdir.cleanup()
            finally:
                if copied_to_output and not keep_intermediate:
                    try:
                        if haploview_ped.exists():
                            haploview_ped.unlink()
                    except Exception:
                        pass
                    try:
                        if haploview_info.exists():
                            haploview_info.unlink()
                    except Exception:
                        pass
                # remove live LD file if created and we are not keeping intermediates
                try:
                    if 'haploview_ld' in locals() and haploview_ld.exists():
                        haploview_ld.unlink()
                except Exception:
                    pass
        sys.exit(1)

    # Verify Haploview produced the LD file at the expected location
    if not haploview_ld.exists():
        print(
            f"ERROR: Haploview did not produce expected LD file: {haploview_ld}",
            file=sys.stderr,
        )
        if not keep_intermediate and "tempdir" in locals():
            try:
                tempdir.cleanup()
            finally:
                if copied_to_output and not keep_intermediate:
                    try:
                        if haploview_ped.exists():
                            haploview_ped.unlink()
                    except Exception:
                        pass
                    try:
                        if haploview_info.exists():
                            haploview_info.unlink()
                    except Exception:
                        pass
        sys.exit(1)

    # Use the Haploview-written LD file directly as the recombinant_path
    recombinant_path = haploview_ld
    print(f"\tRecombinant file (live in output) → {recombinant_path}")

    # --- Step 3: IQ-TREE (optional) ------------------------------------------
    if not skip_iqtree:
        print("\n[Step 3] Running IQ-TREE...")

        # Use the alignment stem as the prefix for files inside the run folder
        run_prefix = alignment_file.stem
        iq_output_prefix = run_iqtree_dir / run_prefix

        extra_args = [
            "-m",
            "GTR",
            "-bb",
            "1000",
            "-alrt",
            "1000",
            "-nt",
            "AUTO",
        ]

        tree_utils.run_iqtree(
            alignment_file,
            iq_output_prefix,
            extra_args=extra_args,
        )

        contree_local = run_iqtree_dir / f"{run_prefix}.contree"
        time_pairs = tree_utils.parse_lengths(contree_local)
    else:
        print("\n[Step 3] Skipping IQ-TREE as requested (--skip-iqtree).")
        time_pairs = None

    # --- Step 4: SNP matching -------------------------------------------------
    print("\n[Step 4] Matching recombinant SNP pairs...")

    sample_data = matching.read_ped_file(ped_path)
    print("read ped file")
    pairs = matching.read_recombinant_file(recombinant_path)
    print("read recombinant file")

    sample_to_pairs = matching.find_matching_pairs(
        sample_data,
        pairs,
    )

    # --- Step 4: Summarize & export -------------------------------------------

    print("\n[Step 5] Summarizing and exporting results...")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    dataframe = matching.summarize_matches(
        sample_to_pairs,
        time_pairs,
    )

    n_groups = matching.getNumGroups(dataframe)

    colors_prefix = matching.getColors(
        n_groups,
        prefix="background-color",
        alpha=0.4,
    )

    # If the caller didn't set a custom output path, place the Excel into the run directory
    output_excel_path = Path(output_excel)
    if output_excel_path.parent == OUTPUT_DIR and output_excel_path.name == "matches_output.xlsx":
        output_excel_path = run_dir / f"{alignment_file.stem}_matches_output.xlsx"

    matching.save_to_excel(
        dataframe,
        output_excel_path,
        colors_prefix,
    )

    colors_no_prefix = matching.getColors(n_groups)

    # Color tree only if IQ-TREE was run and contree exists
    if not skip_iqtree:
        try:
            tree_utils.color_tree(
                dataframe,
                contree_local,
                run_dir / run_prefix,
                colors_no_prefix,
            )
        except Exception:
            print("Warning: Failed to color tree (IQ-TREE output may be missing).")

    # --- Cleanup temp files if needed ---
    if not keep_intermediate and "tempdir" in locals():
        try:
            tempdir.cleanup()
        finally:
            if 'copied_to_output' in locals() and copied_to_output and not keep_intermediate:
                try:
                    if haploview_ped.exists():
                        haploview_ped.unlink()
                except Exception:
                    pass
                try:
                    if haploview_info.exists():
                        haploview_info.unlink()
                except Exception:
                    pass
            # remove Haploview LD file if it lives in run_dir and intermediates not kept
            try:
                if 'recombinant_path' in locals() and recombinant_path.exists():
                    # only delete if it's the live output (in run_dir)
                    if recombinant_path.parent == run_dir:
                        recombinant_path.unlink()
            except Exception:
                pass
        print("Temporary intermediate files deleted.")

    print("\nDone.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the ViralRecombinant matching pipeline.", add_help=True
    )

    parser.add_argument(
        "-a",
        "--alignment",
        required=True,
        help="Path to alignment FASTA file.",
    )

    parser.add_argument(
        "--haploview-jar",
        required=False,
        help="Path to haploview.jar file (required unless using --pedinfo-only).",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=OUTPUT_DIR / "matches_output.xlsx",
        help="Output Excel filename.",
    )

    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep intermediate PED, INFO, and recombinant files (default: use temp files and delete after run)",
    )

    parser.add_argument(
        "--pedinfo-only",
        action="store_true",
        help="Only generate PED and INFO files from FASTA, then exit.",
    )

    parser.add_argument(
        "--skip-iqtree",
        action="store_true",
        help="Skip running IQ-TREE and tree coloring (useful if you only want Haploview outputs).",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    # If the user requests only PED/INFO generation, always keep intermediates
    # so the generated files are available for inspection.
    if args.pedinfo_only and not args.keep_intermediate:
        print("Note: --pedinfo-only implies --keep-intermediate; preserving PED/INFO files.")
        args.keep_intermediate = True
    if not args.pedinfo_only and not args.haploview_jar:
        print(
            "ERROR: --haploview-jar is required unless using --pedinfo-only.",
            file=sys.stderr,
        )
        sys.exit(2)
    main(
        alignment_file=args.alignment,
        haploview_jar=args.haploview_jar,
        output_excel=args.output,
        keep_intermediate=args.keep_intermediate,
        pedinfo_only=args.pedinfo_only,
        skip_iqtree=args.skip_iqtree,
    )
