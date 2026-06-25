"""
This file contains functions for reading PED and LD files and processing recombinant SNP pairs

Author: Alex Poyer
"""

import os
import pandas as pd
import re
import tree_utils
import matplotlib.pyplot as plt


# This function may further be optimized by utilizing pandas dataframes. Works for small datasets
def read_ped_file(ped_file, ped_header="Position"):
    """
    Reads a PED file and extracts genotype data for each sample, merging duplicate positions.

    Args:
        ped_file (str): Path to the PED text file.
        ped_header (str, optional): Expected header for the first column. Defaults to "Position".

    Raises:
        ValueError: If the PED file is empty or has an invalid format.

    Returns:
        dict: sample_data
            A dictionary mapping each sample name (str) to its genotype data.
            Duplicate positions are merged using logical OR.
            Format:
                {
                    "sample_name1": {position1: value1, position2: value2, ...},
                    "sample_name2": {position1: value1, position2: value2, ...},
                    ...
                }
            Where positions are SNP positions (int) and values are genotype calls (int, e.g., 0 or 1).
    """
    sample_data = {}

    with open(ped_file, "r") as f:
        lines = f.readlines()
        if not lines:
            raise ValueError("ped.txt is empty")

        first_line = lines[0].strip()
        parts = re.split(r"\s+", first_line)
        parts = [p.strip() for p in parts]

        # =========================================================
        # 1. LEGACY FORMAT LOGIC (Has 'Position' Header)
        # =========================================================
        if parts[0].lower() == ped_header.lower():
            positions = [int(pos) for pos in parts[1:]]
            pos_indices = {}
            for idx, pos in enumerate(positions):
                pos_indices.setdefault(pos, []).append(idx)

            for line in lines[1:]:
                parts = re.split(r"\s+", line.strip())
                parts = [p.strip() for p in parts]
                if len(parts) != len(positions) + 1:
                    continue
                sample_name = parts[0]
                try:
                    values = [int(val) for val in parts[1:]]
                    merged = {}
                    for pos, idxs in pos_indices.items():
                        merged[pos] = int(any(values[i] for i in idxs))
                    sample_data[sample_name] = merged
                except ValueError:
                    continue
            return sample_data

        # =========================================================
        # 2. NEW LINKAGE FORMAT LOGIC (No Header, 6 Preamble Cols)
        # =========================================================
        else:
            # Infer the corresponding INFO file path to get the positions
            info_file = (
                str(ped_file)
                .replace("_ped.txt", "_info.txt")
                .replace("ped.txt", "info.txt")
            )

            if not os.path.exists(info_file):
                raise ValueError(
                    f"Linkage format detected, but matching INFO file not found at: {info_file}"
                )

            positions = []
            with open(info_file, "r") as f_info:
                for line in f_info:
                    info_parts = line.strip().split()
                    if len(info_parts) >= 2:
                        positions.append(
                            int(info_parts[1])
                        )  # Position is col 2 in Linkage INFO

            pos_indices = {}
            for idx, pos in enumerate(positions):
                pos_indices.setdefault(pos, []).append(idx)

            for line in lines:
                parts = re.split(r"\s+", line.strip())
                if len(parts) < 6:
                    continue

                sample_name = parts[1]  # In Linkage format, IndID is the 2nd column
                genotypes = parts[6:]  # Genotype pairs start at the 7th column

                values = []
                # Linkage genotypes are in pairs (e.g., 1 1 or 1 2)
                # We convert them back to 0 (match) and 1 (mismatch) for ViralRecombinant logic
                for i in range(len(positions)):
                    try:
                        a1 = genotypes[i * 2]
                        a2 = genotypes[i * 2 + 1]
                        if a1 == "1" and a2 == "1":
                            values.append(0)  # Consensus match
                        else:
                            values.append(1)  # Mismatch
                    except IndexError:
                        values.append(0)

                merged = {}
                for pos, idxs in pos_indices.items():
                    merged[pos] = int(any(values[i] for i in idxs))
                sample_data[sample_name] = merged

            return sample_data


# primative approach for parsing and storing data -> efficient for small/medium datasets
def read_recombinant_file(recombinant_file, CIHi_filter=0.9):
    """
    Parses a Haploview LD export file and returns all recombinant SNP position pairs.
    with CIHi < 0.9 by default.

    Args:
        recombinant_file (str): Path to the LD text file.
        CIHi_filter (float, optional): maximum CIHi for a pair to be considered (exclusive)
            default: 0.9

    Raises:
        ValueError: If the file does not contain columns L1, L2, and CIhi.

    Returns:
        List: pairs
            A list of recombinant SNP position pairs to check for in each sample.
            Format:
                [
                    (l1_1, l2_1),
                    (l1_2, l2_2),
                    ...
                ]
            Where each tuple contains two SNP positions (int).
    """
    try:
        df = pd.read_csv(
            recombinant_file,
            sep="\t",
            usecols=["L1", "L2", "CIhi"],
            dtype={"CIhi": float},  # Parse as float immediately
        )
    except ValueError:
        raise ValueError("recombinant_file must contain columns: L1, L2, CIhi")

    if CIHi_filter is not None:
        df = df[df["CIhi"] < CIHi_filter]

    df["L1"] = df["L1"].astype(str).str.replace("SNP_", "", regex=False)
    df["L2"] = df["L2"].astype(str).str.replace("SNP_", "", regex=False)

    # Safely convert to numeric (coercing any weird Haploview artifacts to NaN, then dropping them)
    df["L1"] = pd.to_numeric(df["L1"], errors="coerce")
    df["L2"] = pd.to_numeric(df["L2"], errors="coerce")
    df = df.dropna(subset=["L1", "L2"])

    pairs = list(zip(df["L1"].astype(int), df["L2"].astype(int)))

    return pairs


def find_matching_pairs(sample_data, pairs):
    """Finds which recombinant SNP pairs are present in each sample

    Args:
        sample_data (dict): A dictionary mapping each sample name (str) to its genotype data.
            See read_ped_file for format.
        pairs (list): A list of recombinant SNP position pairs to check for in each sample.
            See read_recombinant_file for format.

    Returns:
        dict: sample_to_pairs
            A dictionary mapping each sample name (str) to the set of recombinant pairs present in that sample.
            Format:
                {
                    "sample_name1": ((l1_1, l2_1), (l1_2, l2_2), ...),
                    "sample_name2": ((l1_3, l2_3), ...),
                    ...
                }
            Where each value is a tuple of recombinant pairs (each pair is a tuple of two SNP positions, int).
    """
    sample_to_pairs = {}
    for sample_name in sorted(sample_data.keys()):
        # print(f"matching {sample_name}") # Debug
        sample_dict = sample_data[sample_name]
        matching_pairs = []
        for l1, l2 in pairs:
            if sample_dict.get(l1, 0) == 1 and sample_dict.get(l2, 0) == 1:
                matching_pairs.append((l1, l2))
        sample_to_pairs[sample_name] = tuple(sorted(matching_pairs))
    # print("pairs:", sample_to_pairs) # Debug
    return sample_to_pairs


def summarize_matches(sample_to_pairs, phylo_times=None):
    """
    Summarizes which samples share the same set of recombinant pairs.

    Args:
        sample_to_pairs (dict): A dictionary mapping each sample name (str) to the set of recombinant pairs present in that sample.
            See find_matching_pairs for format.
        phylo_times (dict, optional): A dictionary mapping each sample name (str) to branch length / phylogenetic time (int)

    Returns:
        DataFrame: df
            A dataframe containing rows :
                [epi_isl identifier, sample name, number of pairs, uniqueness, pair identities, samples it shares identities with]
    """
    pairs_to_samples = {}
    # uses pairset as the key and samples with that pairset as the values
    for sample_name, matching_pairs in sample_to_pairs.items():
        pairs_to_samples.setdefault(matching_pairs, []).append(sample_name)

    rows = []
    for pair_set, samples in pairs_to_samples.items():
        num_pairs = len(pair_set)
        pair_identities = "; ".join([f"L1={l1}, L2={l2}" for l1, l2 in pair_set])

        # gathering phylogenetic times for this unique pair set
        sample_phylos = []
        for sample in samples:
            epi_isl = tree_utils.extract_epi_isl(sample)
            phylo_time = (
                phylo_times.get(epi_isl, float("inf")) if phylo_times else float("inf")
            )
            sample_phylos.append((sample, epi_isl, phylo_time))
        # sort by phylogenetic time (ascending by default)
        sample_phylos.sort(key=lambda x: x[2])

        seq_to_group = {}
        group_sets = []

        # Assign labels based on number of pairs and phylogenetic times.
        # - If num_pairs == 0 => mark as Non-recombinant
        # - Otherwise, find the minimal phylogenetic time; if multiple samples share
        #   the minimal time, mark those as Ambiguous; a single minimal gets Unique;
        #   all others are Progeny.
        if num_pairs == 0:
            for sample, epi_isl, phylo_time in sample_phylos:
                unique = "Non-recombinant"
                shared_with = ""
                rows.append(
                    [
                        epi_isl,
                        sample,
                        num_pairs,
                        unique,
                        pair_identities,
                        shared_with,
                        phylo_time if phylo_time != float("inf") else "",
                    ]
                )
        else:
            # determine minimal phylogenetic time among the group
            min_time = min(t for _, _, t in sample_phylos)
            min_count = sum(1 for _, _, t in sample_phylos if t == min_time)

            for idx, (sample, epi_isl, phylo_time) in enumerate(sample_phylos):
                if len(samples) == 1:
                    unique = "Unique"
                else:
                    if phylo_time == min_time and min_count > 1:
                        unique = "Ambiguous"
                    elif phylo_time == min_time:
                        unique = "Unique"
                    else:
                        unique = "Progeny"
                shared_with = ", ".join([s for _, s, _ in sample_phylos if s != epi_isl])
                rows.append(
                    [
                        epi_isl,
                        sample,
                        num_pairs,
                        unique,
                        pair_identities,
                        shared_with,
                        phylo_time if phylo_time != float("inf") else "",
                    ]
                )
    df = pd.DataFrame(
        rows,
        columns=[
            "epi_isl",
            "Sample",
            "Num Pairs",
            "UniqueRecombinant",
            "Pair Identities",
            "Shared With",
            "Phylogenetic Time",
        ],
    )
    # Apply grouping by shared recombinant pairs associated with group numbers
    group_by_relationship(df)
    return df


def group_by_relationship(df):
    """
        Groups samples by shared recombinant pairs and returns a
        dictionary mapping group numbers to lists of epi_isl identifiers.

    Args:
        df (DataFrame): A dataframe containing rows with 'epi_isl' and 'Shared With' columns

    Returns:
        dict: A dictionary mapping group numbers to lists of epi_isl identifiers
    """
    seq_to_group = {}
    group_sets = []

    for esl_id, shared_list, num_pairs in df[
        ["epi_isl", "Shared With", "Num Pairs"]
    ].values:

        # FIX: Non-recombinants (0 pairs) should not be placed into a shared group
        if num_pairs == 0:
            continue

        if shared_list:
            ids = [x.strip() for x in shared_list.split(",")]
            ids.append(esl_id)
            group_set = frozenset(ids)
        else:
            continue

        if group_set not in group_sets:
            group_sets.append(group_set)
        group_num = group_sets.index(group_set) + 1

        seq_to_group[esl_id] = group_num

    df["Group"] = df["epi_isl"].map(seq_to_group)


def save_to_excel(df, output_excel, colors):
    """
    Saves the DataFrame to an Excel file


    Args:
        excel_data (list): list with summarized sample data
        output_excel (str): output excel file name
    """
    export_cols = [
        "Sample",
        "epi_isl",
        "Num Pairs",
        "Phylogenetic Time",
        "UniqueRecombinant",
        "Group",
        "Shared With",
        "Pair Identities",
    ]
    styled = highlight_group_cells(df, colors)
    styled.to_excel(output_excel, index=False, columns=export_cols)
    print(f"Results saved to {output_excel}")


def highlight_group_cells(df, colors):
    """
    Applys background color to dataframe based on groups defined by shared recombinant pairs.
    Each group gets a different color, and samples in the same group share the same color.

    Args:
        df (DataFrame): A dataframe containing rows with 'group' column that
                defines group membership based on shared recombinant pairs.

    Returns:
        (DataFrame): A dataframe that has been styled with background colors based on group membership.
    """
    groups = df["Group"].dropna().unique()
    groups_sorted = sorted(groups, key=lambda x: (isinstance(x, str), x))

    group_to_color = {group: colors[i] for i, group in enumerate(groups_sorted)}

    def row_style(row):
        # FIX: Force sequences with 0 recombinant pairs to be styled light gray
        # Use a light gray background and dark text so cells remain readable
        if row["Num Pairs"] == 0:
            return ["background-color: #D3D3D3; color: #000000"] * len(row)

        color = group_to_color.get(row["Group"], "")
        return [color] * len(row)

    return df.style.apply(row_style, axis=1)


def getColors(num_groups, alpha=1.0, prefix=None):
    """
    Generates a list of distinct, visually lightened color strings for styling groups in a DataFrame.

    Args:
        num_groups (int): The number of unique groups/colors needed.
        alpha (float, optional): The alpha value for blending with white (0-1).
            defualts to 1.0 (no blending). Lower values produce lighter colors.
        prefix (str, optional): Optional prefix for each color string (e.g., "background-color").
            If provided, the prefix will be prepended to each color value, separated by ": ".

    Returns:
        list of str: List of color strings in hex format, optionally prefixed, suitable for use in pandas Styler or CSS.
    """
    cmap = plt.get_cmap("tab20" if num_groups <= 20 else "hsv")

    # Uses discrete indices instead of linspace to ensure distinct colors
    n_colors = cmap.N  # total colors in the colormap

    if prefix is not None:
        prefix = prefix.strip() + ": "
    else:
        prefix = ""

    colors = [
        prefix + blend_with_white(*cmap(i % n_colors)[:3], alpha=alpha)
        for i in range(num_groups)
    ]
    return colors


def blend_with_white(r, g, b, alpha=0.4):
    """
    Blends an RGB color with white to produce a lighter shade, and returns a CSS background-color string.

    Args:
        r (float): Red component of the color (0-1)
        g (float): Green component of the color (0-1)
        b (float): Blue component of the color (0-1)
        alpha (float, optional): Alpha component for blending (0-1). Defaults to 0.4.

    Returns:
        str: color-string in format '#RRGGBB' where RRGGBB is the blended color in hexadecimal.
    """
    r2 = int((r * alpha + 1 * (1 - alpha)) * 255)
    g2 = int((g * alpha + 1 * (1 - alpha)) * 255)
    b2 = int((b * alpha + 1 * (1 - alpha)) * 255)
    return f"#{r2:02X}{g2:02X}{b2:02X}"


def getNumGroups(df):
    """
    Extracts number of groups from a dataframe

    Args:
        df (DataFrame): A dataframe containing rows with 'group' column that
                defines group membership based on shared recombinant pairs.
    """
    groups = df["Group"].dropna().unique()
    groups_sorted = sorted(groups, key=lambda x: (isinstance(x, str), x))
    return len(groups_sorted)
