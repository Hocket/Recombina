"""
This file contains utilities for for phylogenetic tree analysis and modification

Author: Alex Poyer
"""

import subprocess
import re
from Bio import Phylo


def run_iqtree(input_alignment, output_prefix, iqtree_path="iqtree3", extra_args=None):
    """
    Runs IQ-TREE on the given alignment file.

    Args:
        input_alignment (str): Path to the input alignment file (FASTA, PHYLIP, etc.).
        output_prefix (str): Prefix for output files.
        iqtree_path (str): Path to the IQ-TREE executable (default: "iqtree3").
        extra_args (list): Additional command-line arguments for IQ-TREE.

    Returns:
        int: The return code from the IQ-TREE process.
    """
    cmd = [str(iqtree_path), "-s", str(input_alignment), "-pre", str(output_prefix)]
    if extra_args:
        cmd.extend(map(str, extra_args))
    print("Running:", " ".join(cmd))
    result = subprocess.run(cmd)
    return result.returncode


def parse_lengths(tree_file):
    tree = Phylo.read(tree_file, "newick")
    time_pairs = {}
    for clade in tree.get_terminals():
        epi_isl = extract_epi_isl(clade.name)
        if epi_isl:
            time_pairs[epi_isl] = clade.branch_length
    return time_pairs


def extract_epi_isl(name):
    """Extracts the EPI_ISL identifier from a sample name string."""
    match = re.search(r"EPI_ISL_\d+", name)
    return match.group(0) if match else name


def color_tree(data, tree_file, output_file, colors):
    """writes comments on branches to state color

    Args:
        data (dataFrame): pandas dataframe representing a tree
        tree_file (string): A newick style tree file to be colored
        output_file (string): The output file for the colored newick tree
    """

    # pass through and add color to tree
    tree = Phylo.read(tree_file, "newick")
    nwk = str(output_file) + ".nwk"
    Phylo.write(tree, nwk, "newick")
    print("Tree saved to", nwk)

    df = data.set_index("epi_isl")

    groups = df["Group"].dropna().unique()
    groups_sorted = sorted(groups)
    group_to_color = {group: colors[i] for i, group in enumerate(groups_sorted)}

    for clade in tree.get_terminals():
        sample_name = clade.name
        epi_isl = extract_epi_isl(sample_name)
        if epi_isl is None:
            continue
        clade.name = epi_isl
        unique_recombinant = (
            df.at[epi_isl, "UniqueRecombinant"] if epi_isl in df.index else None
        )
        
        # COLORS USE HEX CODE
        # FIX: Ensure unique_recombinant is a string before calling .lower()
        if isinstance(unique_recombinant, str) and unique_recombinant.lower() == "unique":
            clade.comment = "&!color=#ff0000"
            continue

        group = df.at[epi_isl, "Group"] if epi_isl in df.index else None
        if group in group_to_color:
            clade.comment = f"&!color={group_to_color[group]}"

    nexus = str(output_file) + ".nexus"
    Phylo.write(tree, nexus, "nexus")
    print("Color tree saved to", nexus)