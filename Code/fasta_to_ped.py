"""
fasta_to_ped.py

Converts a multiple sequence alignment (FASTA) into Haploview-compatible
PED and INFO files for use with the ViralRecombinant pipeline.

Pipeline steps
--------------
1. Parse FASTA and collect all alignment columns (positions).
2. Detect coverage boundaries — find the smallest range where every sample 
   has at least one valid (non-gap) character, trim to this region.
3. Keep only positions where at least one sequence differs from the rest
   (variable positions).
4. Drop positions that contain any ambiguous base — anything outside
   {A, T, C, G, -} (case-insensitive).
5. Drop positions where the total count of valid bases is less than the
   number of sequences (incomplete columns).
6. Drop singleton positions — positions where only one sequence has a
   non-consensus base.
7. Split positions that contain more than one non-consensus allele type into
   one virtual position per minority allele.  Each virtual position contrasts
   exactly one minority allele against the consensus (everything else).
   Example — Pos 900: (2A, 1T, 300G) → two virtual positions:
       Pos 900[A]: 2A vs 301G
       Pos 900[T]: 1T vs 302G
8. Determine the consensus base at each surviving (virtual) position.
9. Encode each sample's base at each (virtual) position:
       matches consensus  →  11
       differs from consensus →  12
10. Write PED file:  tab-separated, header row = "Position <pos1> <pos2> ..."
                    data rows   = "<sample_name> <val1> <val2> ..."
                    (Positions reported in original alignment coordinates)
11. Write INFO file: two columns — position (bp) and a marker ID
                    (used by Haploview for LD computation)
                    (Positions reported in original alignment coordinates)

Author: ViralRecombinant contributors
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

VALID_BASES: frozenset[str] = frozenset("ATCG-")

# Tie-breaking priority when two bases are equally frequent at a position
BASE_PRIORITY: dict[str, int] = {"A": 0, "T": 1, "C": 2, "G": 3, "-": 4}

# Haploview encoding
MATCH_CODE = 11  # base == consensus
MISMATCH_CODE = 12  # base != consensus


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------


def parse_fasta(fasta_path: str | Path) -> dict[str, str]:
    """
    Parses a FASTA file and returns an ordered dict of {sample_name: sequence}.

    All sequences are uppercased. Raises ValueError if the file is empty,
    no sequences are found, or sequences differ in length (not a valid MSA).

    Args:
        fasta_path: Path to the FASTA file.

    Returns:
        dict mapping sample name → aligned sequence string (all uppercase).
    """
    sequences: dict[str, str] = {}
    current_name: Optional[str] = None
    current_seq: list[str] = []

    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    with open(fasta_path, "r") as fh:
        for line_num, raw_line in enumerate(fh, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Save the previous record
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq).upper()
                current_name = line[1:].strip()  # strip the ">"
                current_seq = []
                if not current_name:
                    raise ValueError(f"Empty sequence header at line {line_num}")
            else:
                if current_name is None:
                    raise ValueError(
                        f"Sequence data before any header at line {line_num}"
                    )
                current_seq.append(line)

    # Save the last record
    if current_name is not None:
        sequences[current_name] = "".join(current_seq).upper()

    if not sequences:
        raise ValueError(f"No sequences found in {fasta_path}")

    # Validate that all sequences are the same length (required for MSA)
    lengths = {name: len(seq) for name, seq in sequences.items()}
    unique_lengths = set(lengths.values())
    if len(unique_lengths) > 1:
        bad = [(n, l) for n, l in lengths.items() if l != max(unique_lengths)]
        raise ValueError(
            f"Sequences differ in length — not a valid MSA. "
            f"Example mismatches: {bad[:3]}"
        )

    return sequences


# ---------------------------------------------------------------------------
# Coverage-based gap filtering
# ---------------------------------------------------------------------------


def find_coverage_boundaries(sequences: dict[str, str]) -> tuple[int, int]:
    """
    Finds the smallest range where every sample has at least one valid character.
    
    This function:
    1. Finds leading gaps for each sequence (gaps from start until first non-gap)
    2. Finds trailing gaps for each sequence (gaps from end until last non-gap)
    3. Finds the MAX leading gaps (furthest point any sequence starts)
    4. Finds the MIN trailing position (earliest point any sequence ends)
    5. Returns the range where ALL sequences have valid data
    
    Args:
        sequences: {sample_name: aligned_sequence}
    
    Returns:
        Tuple of (start_idx, end_idx) where start_idx is the position after the
        maximum leading gaps, and end_idx is the position where the earliest
        trailing gap starts (0-based, where end_idx is exclusive).
    """
    if not sequences:
        return 0, 0
    
    seq_length = len(next(iter(sequences.values())))
    
    # Calculate leading and trailing gaps for each sequence
    max_leading = 0
    min_trailing_start = seq_length
    
    gap_details = []
    
    for sample_name, seq in sequences.items():
        # Count leading gaps
        leading = 0
        for i in range(len(seq)):
            if seq[i] == '-':
                leading += 1
            else:
                break
        
        # Count trailing gaps
        trailing = 0
        for i in range(len(seq) - 1, -1, -1):
            if seq[i] == '-':
                trailing += 1
            else:
                break
        
        # Position where this sequence's trailing gaps start (end of valid region)
        trailing_start_pos = seq_length - trailing
        
        max_leading = max(max_leading, leading)
        min_trailing_start = min(min_trailing_start, trailing_start_pos)
        
        gap_details.append({
            'name': sample_name,
            'leading': leading,
            'trailing': trailing,
            'valid_end': trailing_start_pos
        })
    
    # Find sequences with the extremes for reporting
    max_leading_seq = max(gap_details, key=lambda x: x['leading'])
    min_end_seq = min(gap_details, key=lambda x: x['valid_end'])
    
    start_idx = max_leading
    end_idx = min_trailing_start
    
    print(
        f"Coverage boundaries detected (range where ALL samples have valid data):\n"
        f"  Sequence with most leading gaps: {max_leading_seq['name']}\n"
        f"    Leading gaps: {max_leading_seq['leading']}\n"
        f"  Sequence with most trailing gaps: {min_end_seq['name']}\n"
        f"    Trailing gaps: {seq_length - min_end_seq['valid_end']}\n"
        f"  Valid region: position {start_idx + 1} to {end_idx} (1-based, {end_idx - start_idx} bp)\n"
    )
    
    return start_idx, end_idx


# ---------------------------------------------------------------------------
# Position-level filtering
# ---------------------------------------------------------------------------


def _column_bases(sequences: dict[str, str], pos: int) -> list[str]:
    """Returns the list of bases at alignment position `pos` across all samples."""
    return [seq[pos] for seq in sequences.values()]


def find_variable_positions(sequences: dict[str, str]) -> list[int]:
    """
    Returns the 1-based alignment positions where at least one sequence
    has a different base from the rest.

    Args:
        sequences: {sample_name: aligned_sequence}

    Returns:
        Sorted list of 1-based variable positions.
    """
    if not sequences:
        return []

    seq_length = len(next(iter(sequences.values())))
    variable: list[int] = []

    for pos in range(seq_length):  # 0-based internally
        bases = _column_bases(sequences, pos)
        if len(set(bases)) > 1:
            variable.append(pos + 1)  # convert to 1-based for output

    return variable


def filter_ambiguous_positions(
    sequences: dict[str, str], positions: list[int]
) -> list[int]:
    """
    Removes positions that contain any base outside {A, T, C, G, -}.

    Args:
        sequences: {sample_name: aligned_sequence}
        positions: 1-based positions to evaluate.

    Returns:
        Filtered list of 1-based positions.
    """
    clean: list[int] = []
    for pos in positions:
        bases = set(_column_bases(sequences, pos - 1))  # 0-based access
        if bases.issubset(VALID_BASES):
            clean.append(pos)
    return clean


def filter_incomplete_positions(
    sequences: dict[str, str], positions: list[int]
) -> list[int]:
    """
    Removes positions where the total count of valid bases is less than
    the number of sequences (i.e. some bases are ambiguous/missing even
    after ambiguity filtering — belt-and-suspenders check).

    Args:
        sequences: {sample_name: aligned_sequence}
        positions: 1-based positions to evaluate.

    Returns:
        Filtered list of 1-based positions.
    """
    n_seqs = len(sequences)
    clean: list[int] = []
    for pos in positions:
        bases = _column_bases(sequences, pos - 1)
        valid_count = sum(1 for b in bases if b in VALID_BASES)
        if valid_count == n_seqs:
            clean.append(pos)
    return clean


def filter_singleton_positions(
    sequences: dict[str, str], positions: list[int]
) -> list[int]:
    """
    Removes positions where only one sequence differs from the consensus
    (singletons add noise and are excluded by Haploview convention).

    Args:
        sequences: {sample_name: aligned_sequence}
        positions: 1-based positions to evaluate.

    Returns:
        Filtered list of 1-based positions.
    """
    clean: list[int] = []
    for pos in positions:
        bases = _column_bases(sequences, pos - 1)
        consensus = _consensus_base(bases)
        n_different = sum(1 for b in bases if b != consensus)
        if n_different > 1:
            clean.append(pos)
    return clean


def trim_sequences_to_coverage(
    sequences: dict[str, str], start_idx: int, end_idx: int
) -> dict[str, str]:
    """
    Trims all sequences to the region defined by coverage boundaries.
    
    Args:
        sequences: {sample_name: aligned_sequence}
        start_idx: 0-based start position (inclusive)
        end_idx: 0-based end position (exclusive)
    
    Returns:
        Trimmed sequences dictionary.
    """
    trimmed = {}
    for sample_name, seq in sequences.items():
        trimmed[sample_name] = seq[start_idx:end_idx]
    return trimmed


# ---------------------------------------------------------------------------
# Multi-allelic splitting (Step 6)
# ---------------------------------------------------------------------------


def split_multiallelic_positions(
    sequences: dict[str, str], positions: list[int]
) -> list[tuple[int, str]]:
    """
    Expands positions with more than one non-consensus allele into one virtual
    position per minority allele.

    For a position where the consensus is G and there are both A's and T's as
    minority alleles, this yields two virtual positions:
        (pos, 'A')  — contrasts A against everything else (treated as G)
        (pos, 'T')  — contrasts T against everything else (treated as G)

    Positions with only a single non-consensus allele type pass through as-is,
    represented as (pos, minority_allele).

    Args:
        sequences: {sample_name: aligned_sequence}
        positions: 1-based positions surviving all prior filters.

    Returns:
        Ordered list of (1-based position, minority_allele) tuples.
        Positions with a single minority allele type appear once.
        Multi-allelic positions appear once per distinct minority allele.
    """
    result: list[tuple[int, str]] = []
    for pos in positions:
        bases = _column_bases(sequences, pos - 1)
        consensus = _consensus_base(bases)
        minority_alleles = sorted(
            {b for b in bases if b != consensus},
            key=lambda b: BASE_PRIORITY.get(b, 99),
        )
        for allele in minority_alleles:
            result.append((pos, allele))
    return result


# ---------------------------------------------------------------------------
# Consensus & encoding
# ---------------------------------------------------------------------------


def _consensus_base(bases: list[str]) -> str:
    """
    Returns the most-frequent base in `bases`.
    Ties are broken by BASE_PRIORITY (A > T > C > G > -).

    Args:
        bases: list of single-character bases (already uppercased).

    Returns:
        Single-character consensus base.
    """
    counts = Counter(b for b in bases if b in VALID_BASES)
    if not counts:
        return "N"  # should not happen after filtering, but defensive
    max_count = max(counts.values())
    candidates = [b for b, c in counts.items() if c == max_count]
    # Sort by priority to break ties deterministically
    candidates.sort(key=lambda b: BASE_PRIORITY.get(b, 99))
    return candidates[0]


def compute_consensus(
    sequences: dict[str, str], virtual_positions: list[tuple[int, str]]
) -> dict[tuple[int, str], str]:
    """
    Computes the consensus base for each virtual position.

    For a virtual position (pos, minority_allele) the consensus is the
    most-frequent base when all bases *other than* minority_allele are
    counted together — i.e. the majority allele at that position.

    Args:
        sequences: {sample_name: aligned_sequence}
        virtual_positions: list of (1-based position, minority_allele) tuples.

    Returns:
        dict mapping (1-based position, minority_allele) → consensus base.
    """
    consensus: dict[tuple[int, str], str] = {}
    for pos, minority_allele in virtual_positions:
        bases = _column_bases(sequences, pos - 1)
        # Consensus = most frequent base ignoring the minority allele
        majority_bases = [b for b in bases if b != minority_allele]
        consensus[(pos, minority_allele)] = _consensus_base(majority_bases)
    return consensus


def encode_sample(
    sequence: str,
    virtual_positions: list[tuple[int, str]],
    consensus: dict[tuple[int, str], str],
) -> list[int]:
    """
    Encodes a single sample's bases at each virtual position using the
    Haploview convention: 11 if the base matches the consensus, 12 if it
    differs.

    For a virtual position (pos, minority_allele), the sample is encoded as:
        12  if its base == minority_allele
        11  otherwise (base matches the majority / consensus)

    Args:
        sequence: The sample's aligned sequence string.
        virtual_positions: list of (1-based position, minority_allele) tuples.
        consensus: Mapping of (1-based position, minority_allele) → consensus base.

    Returns:
        List of integer codes (11 or 12) in virtual-position order.
    """
    codes: list[int] = []
    for pos, minority_allele in virtual_positions:
        base = sequence[pos - 1]  # 0-based access
        code = MISMATCH_CODE if base == minority_allele else MATCH_CODE
        codes.append(code)
    return codes


# ---------------------------------------------------------------------------
# File writers
# ---------------------------------------------------------------------------


def write_ped_file(
    sequences: dict[str, str],
    virtual_positions: list[tuple[int, str]],
    consensus: dict[tuple[int, str], str],
    output_path: str | Path,
) -> None:
    """
    Writes a Haploview-compatible Linkage format PED text file.
    No header row.
    6 preamble columns: FamID, IndID, PatID, MatID, Sex, Pheno.
    Genotypes are tab-separated allele pairs (1 1 or 1 2).

    Column order follows virtual_positions; multi-allelic positions appear
    as consecutive columns, one per minority allele.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as fh:
        # One row per sample (Header row intentionally omitted for Haploview)
        index = 1 # Start sample index at 1 for FamID/IndID
        for sample_name, sequence in sequences.items():
            codes = encode_sample(sequence, virtual_positions, consensus)

            # Format 11 -> "1 1" and 12 -> "1 2" (single cell per site: space-separated alleles)
            formatted_genotypes = ["1 1" if c == MATCH_CODE else "1 2" for c in codes]

            # Force 6 preamble columns. Using sample_name for IndID
            row = f"{index}\t{sample_name}\t0\t0\t0\t1\t" + "\t".join(
                formatted_genotypes
            )

            fh.write(row + "\n")
            index += 1

    print(
        f"PED file written → {output_path}  ({len(virtual_positions)} columns, {len(sequences)} samples)"
    )


def write_info_file(
    virtual_positions: list[tuple[int, str]],
    output_path: str | Path,
) -> None:
    """
    Writes a Haploview INFO file.
    Format (two tab-separated columns, no header):
        <position_bp>  <position_bp>

    Multi-allelic splits at the same genomic position each get their own
    row with the same bp coordinate, matching Haploview's column order.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as fh:
        for pos, _allele in virtual_positions:
            fh.write(f"{pos}\t{pos}\n")

    print(f"INFO file written → {output_path}  ({len(virtual_positions)} markers)")


# ---------------------------------------------------------------------------
# Summary / diagnostics
# ---------------------------------------------------------------------------


def _pluralize(n: int, word: str) -> str:
    return f"{n} {word}{'s' if n != 1 else ''}"


def print_filter_summary(
    n_sequences: int,
    original_length: int,
    trimmed_length: int,
    n_variable: int,
    n_after_ambiguous: int,
    n_after_incomplete: int,
    n_after_singleton: int,
    n_after_split: int,
) -> None:
    """Prints a human-readable summary of how many positions survived each filter."""
    dropped_by_coverage = original_length - trimmed_length
    dropped_nonvariable = trimmed_length - n_variable
    dropped_ambiguous = n_variable - n_after_ambiguous
    dropped_incomplete = n_after_ambiguous - n_after_incomplete
    dropped_singleton = n_after_incomplete - n_after_singleton
    added_by_split = n_after_split - n_after_singleton

    print(
        f"\n{'='*60}\n"
        f"  FASTA → PED/INFO  Filter Summary\n"
        f"{'='*60}\n"
        f"  Sequences in alignment         : {n_sequences:>7}\n"
        f"  Original alignment columns     : {original_length:>7}\n"
        f"  — Coverage trim (removed)      : {dropped_by_coverage:>7}\n"
        f"  After coverage filter          : {trimmed_length:>7}\n"
        f"  — Non-variable (removed)       : {dropped_nonvariable:>7}\n"
        f"  After variable filter          : {n_variable:>7}\n"
        f"  — Ambiguous base (removed)     : {dropped_ambiguous:>7}\n"
        f"  After ambiguity filter         : {n_after_ambiguous:>7}\n"
        f"  — Incomplete col (removed)     : {dropped_incomplete:>7}\n"
        f"  After completeness filter      : {n_after_incomplete:>7}\n"
        f"  — Singletons (removed)         : {dropped_singleton:>7}\n"
        f"  After singleton filter         : {n_after_singleton:>7}\n"
        f"  + Split multi-allelic          : {added_by_split:>7}\n"
        f"  Final virtual positions        : {n_after_split:>7}\n"
        f"{'='*60}\n"
    )


# ---------------------------------------------------------------------------
# Main pipeline function
# ---------------------------------------------------------------------------


def fasta_to_ped(
    fasta_path: str | Path,
    ped_output: str | Path,
    info_output: str | Path,
    verbose: bool = True,
) -> tuple[list[tuple[int, str]], dict[tuple[int, str], str]]:
    """
    Full pipeline: FASTA → PED + INFO files.

    Args:
        fasta_path:  Path to the input MSA FASTA file.
        ped_output:  Destination path for the PED file.
        info_output: Destination path for the INFO file.
        verbose:     If True, print a filter summary to stdout.

    Returns:
        (virtual_positions, consensus) — the final list of (1-based position,
        minority_allele) tuples and their consensus bases, in case the caller
        needs them for further processing.

    Raises:
        FileNotFoundError: If fasta_path does not exist.
        ValueError: If the FASTA is malformed or no positions survive filtering.
    """
    # Step 1 — Parse
    print(f"Parsing FASTA: {fasta_path}")
    sequences = parse_fasta(fasta_path)
    n_sequences = len(sequences)
    original_length = len(next(iter(sequences.values())))
    print(f"  {n_sequences} sequences, alignment length {original_length} bp")

    # Step 1.5 — Detect coverage boundaries and trim
    print("\nDetecting coverage boundaries...")
    start_idx, end_idx = find_coverage_boundaries(sequences)
    sequences = trim_sequences_to_coverage(sequences, start_idx, end_idx)
    trimmed_length = len(next(iter(sequences.values())))
    
    # Position offset for output files (to convert trimmed coords back to original)
    position_offset = start_idx

    # Step 2 — Variable positions
    positions = find_variable_positions(sequences)
    n_variable = len(positions)

    # Step 3 — Remove ambiguous-base positions
    positions = filter_ambiguous_positions(sequences, positions)
    n_after_ambiguous = len(positions)

    # Step 4 — Remove incomplete columns
    positions = filter_incomplete_positions(sequences, positions)
    n_after_incomplete = len(positions)

    # Step 5 — Remove singletons
    positions = filter_singleton_positions(sequences, positions)
    n_after_singleton = len(positions)

    # Step 6 — Split multi-allelic positions into one virtual position per
    #           minority allele (positions with a single minority allele type
    #           pass through unchanged as a single (pos, allele) entry).
    virtual_positions = split_multiallelic_positions(sequences, positions)
    n_after_split = len(virtual_positions)

    if verbose:
        print_filter_summary(
            n_sequences,
            original_length,
            trimmed_length,
            n_variable,
            n_after_ambiguous,
            n_after_incomplete,
            n_after_singleton,
            n_after_split,
        )

    if n_after_split == 0:
        raise ValueError(
            "No positions survived filtering. "
            "Check your alignment for ambiguous bases or low diversity."
        )

    # Step 7 — Consensus (per virtual position)
    consensus = compute_consensus(sequences, virtual_positions)

    # Adjust virtual positions back to original alignment coordinates
    virtual_positions_original = [
        (pos + position_offset, allele) for pos, allele in virtual_positions
    ]

    # Steps 8-10 — Encode and write
    write_ped_file(sequences, virtual_positions, consensus, ped_output)
    write_info_file(virtual_positions_original, info_output)

    return virtual_positions_original, consensus


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fasta_to_ped",
        description="Convert a multiple sequence alignment FASTA to Haploview PED + INFO files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("fasta", help="Input MSA FASTA file")
    p.add_argument(
        "--ped",
        default=None,
        help="Output PED file path (default: <fasta_stem>_ped.txt alongside input)",
    )
    p.add_argument(
        "--info",
        default=None,
        help="Output INFO file path (default: <fasta_stem>_info.txt alongside input)",
    )
    p.add_argument(
        "--quiet", action="store_true", help="Suppress filter summary output"
    )
    return p


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    fasta_path = Path(args.fasta)

    ped_output = (
        Path(args.ped)
        if args.ped
        else fasta_path.parent / (fasta_path.stem + "_ped.txt")
    )
    info_output = (
        Path(args.info)
        if args.info
        else fasta_path.parent / (fasta_path.stem + "_info.txt")
    )

    try:
        fasta_to_ped(fasta_path, ped_output, info_output, verbose=not args.quiet)
    except (FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()