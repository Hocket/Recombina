# Recombina Toolkit — Complete User Guide

A comprehensive toolkit for analyzing viral sequences to detect recombination events, perform SNP matching, and build phylogenetic trees.

## Table of Contents

1. [What This Tool Does](#what-this-tool-does)
2. [What You Need to Install](#what-you-need-to-install)
3. [Installation Guide](#installation-guide)
4. [Project Structure](#project-structure)
5. [How to Use](#how-to-use)
6. [Understanding the Pipeline](#understanding-the-pipeline)
7. [Input and Output Files](#input-and-output-files)
8. [Example Workflow](#example-workflow)
9. [Troubleshooting](#troubleshooting)
10. [FAQ](#faq)

---

## What This Tool Does

Recombina helps you analyze viral DNA sequences to:

- **Convert sequence alignments to SNP data** — Takes a FASTA file with multiple aligned sequences and extracts variable (SNP) positions
- **Filter low-quality data** — Removes incomplete or ambiguous sequence regions
- **Compute linkage disequilibrium (LD)** — Uses Haploview to calculate how SNP variants are linked together
- **Build phylogenetic trees** — Creates evolutionary trees showing relationships between sequences
- **Generate reports** — Produces Excel summaries and colored tree visualizations

### Who Should Use This?

- Molecular biologists analyzing viral outbreaks or evolution
- Researchers studying viral recombination
- Anyone working with multiple aligned DNA/RNA sequences in FASTA format

### What You'll Get

- **PED/INFO files** — Standard format used by population genetics software (compatible with Haploview, PLINK, etc.)
- **LD heatmaps** — Visual representation of how SNPs are linked
- **Phylogenetic trees** — Evolutionary relationships between your sequences
- **Excel summaries** — Statistical summaries of findings

---

## What You Need to Install

Before you start, gather these tools. Don't worry — most are free and we'll guide you through each step.

### Required

- **Python 3.8 or newer** — The programming language this tool runs on
  - **Option A (Standard):** [Download from python.org](https://www.python.org/downloads/)
    - **On Windows:** Install with "Add Python to PATH" checked ✓
    - **On Mac/Linux:** Usually pre-installed (check: open terminal, type `python3 --version`)
  - **Option B (Recommended for scientists):** Use Miniconda/Anaconda instead
    - [Download Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) (lightweight, recommended)
    - Or [Download Anaconda](https://www.anaconda.com/download) (includes extra tools)
    - See "Installation Guide — Using Conda" section below for setup

- **Java** — Required to run Haploview (a sequence analysis tool)
  - [Download Java from oracle.com](https://www.oracle.com/java/technologies/downloads/)
  - Or use OpenJDK: `brew install openjdk` (Mac) or `apt install default-jdk` (Linux)
  - **Check if installed:** Open terminal/command prompt, type `java -version`

### Included in This Repository

- **Haploview 4.1** — Included, no separate download needed ✓
- **Python dependencies** — Installed automatically in step 3 below

### Optional (For Better Features)

- **IQ-TREE** — Required if you want to build phylogenetic trees
  - [Download from iqtree.org](http://www.iqtree.org/#download)
  - Without this: Pipeline still works, but skips tree building
  - **On Mac:** `brew install iqtree` (If homebrew is installed)
  - **On Linux:** `apt install iqtree`

- **tqdm** — Shows nice progress bars while running
  - Installed automatically with pip in step 3 ✓

---

## Installation Guide

### Step 1: Get the Code

**Option A: Download (Easiest for Beginners)**
1. Go to the project repository
2. Click the green "Code" button → "Download ZIP"
3. Unzip the folder to a location you'll remember (e.g., `Documents/Recombina/`)

**Option B: Clone with Git** (if you know Git)
```bash
git clone https://github.com/Hocket/Recombina.git
cd Recombina
```

---
## Installation Option 1: Using Standard Python + pip (Recommended for Beginners)

### Step 2: Create a Python Environment

A Python "environment" is like a separate workspace for this project. It keeps all the Python packages needed by this tool separate from your other Python projects.

**On Windows (Command Prompt or PowerShell):**

```bash
python -m venv .venv
.venv\Scripts\activate
```

**On Mac or Linux (Terminal):**

```bash
python3 -m venv .venv
source .venv/bin/activate
```

You'll know it worked when you see `(.venv)` at the start of your command prompt/terminal.

### Step 3: Install Python Packages

With the environment activated (you should see `(.venv)` in your prompt), run:

```bash
pip install -r Data/requirements/requirements.txt
pip install tqdm
```

This installs all the Python libraries the tool needs. It may take 1-2 minutes.

### Step 4: Verify Installation

To make sure everything is installed correctly:

```bash
java -version
python --version
```

Both commands should print version numbers (not "command not found").

**If Java doesn't work:**
- On Windows: Add Java to your PATH environment variable (search "Edit environment variables" in Windows)
- On Mac/Linux: Try `which java` to find where it's installed

### Step 5: (Optional) Install IQ-TREE for Tree Building

If you want to build phylogenetic trees:

**On Mac:**
```bash
brew install iqtree
```

**On Linux:**
```bash
apt install iqtree
```

**On Windows:**
- Download from [iqtree.org](http://www.iqtree.org/#download)
- Add the IQ-TREE folder to your PATH (search "Edit environment variables")

Check it's installed:
```bash
iqtree --version
```

---
## Installation Option 2: Using Conda (Recommended for Scientists & Advanced Users)

Conda is a package manager that handles both Python and non-Python packages (like Java, IQ-TREE) in one place. It's especially useful if you:
- Use multiple Python projects with different dependencies
- Work with bioinformatics tools regularly
- Want to easily share your environment setup with collaborators
- Are on a shared computing cluster

### Step 1: Install Miniconda or Anaconda

**Miniconda (Lightweight — Recommended):**
- [Download Miniconda](https://docs.conda.io/projects/miniconda/en/latest/)
- Select the version for your OS (Windows/Mac/Linux) and Python 3.10+
- Run the installer and follow prompts
- **On Mac/Linux:** You may need to run `source ~/miniconda3/bin/activate` after installation

**Anaconda (Full Suite — Includes Extra Tools):**
- [Download Anaconda](https://www.anaconda.com/download)
- Larger download but includes many pre-installed packages
- Same installation process as Miniconda

**Verify installation:**
```bash
conda --version
```

### Step 2: Get the Code

```bash
git clone https://github.com/Hocket/ViralRecombinant.git
cd ViralRecombinant
```

Or download the ZIP file as described in Option 1.

### Step 3: Create a Conda Environment from File

The easiest method — we provide an environment file:

```bash
conda env create -f Data/requirements/environment.yml
```

This creates an environment called `Recombina` with all dependencies.

**If the environment file doesn't exist, create one manually:**

```bash
conda create -n Recombina python=3.10 -y
conda activate Recombina
pip install -r Data/requirements/requirements.txt
pip install tqdm
```

### Step 4: Activate the Environment

Every time you want to use the tool, activate the conda environment:

```bash
conda activate Recombina
```

You'll know it worked when you see `(Recombina)` in your prompt instead of `(base)` or nothing.

### Step 5: Install Additional Tools via Conda (Optional)

**Install Java (if not already installed):**

```bash
conda install -c conda-forge openjdk
```

Verify:
```bash
java -version
```

**Install IQ-TREE (for phylogenetic trees):**

```bash
conda install -c bioconda iqtree
```

Verify:
```bash
iqtree --version
```

### Step 6: Verify Installation

```bash
java -version
python --version
iqtree --version  # only if you installed it
```

All should return version numbers.

### Step 7: Deactivate When Done

When you're finished using the tool, deactivate the environment:

```bash
conda deactivate
```

---

## Comparing Installation Methods

| Feature | pip (Standard) | Conda (Advanced) |
|---------|---|---|
| **Easiest for beginners** | ✓ | - |
| **Works on all systems** | ✓ | ✓ |
| **Handles non-Python packages** | - | ✓ (Java, IQ-TREE) |
| **Good for multiple projects** | - | ✓ |
| **Fast to set up** | ✓ | - (first time slower) |
| **Easy to share setup** | - | ✓ (environment.yml) |
| **Best for clusters/servers** | - | ✓ |

**Recommendation:**
- **Just want to run this tool once?** → Use **pip** (Option 1)
- **Do bioinformatics regularly?** → Use **Conda** (Option 2)
- **Not sure?** → Start with **pip**, switch to conda later if needed

---

## Managing Your Conda Environment

**View all your conda environments:**
```bash
conda env list
```

**Export your environment (to share with others):**
```bash
conda env export > my_environment.yml
```

**Remove an environment (if you don't need it):**
```bash
conda env remove -n Recombina
```

**Update all packages in the environment:**
```bash
conda activate Recombina
conda update --all
```

**List installed packages:**
```bash
conda activate Recombina
conda list
```

---

## Project Structure

Here's what's in your folder and what each part does:

```
Recombina/
│
├── Code/                          # The Code behind Recombina
│   ├── fasta_to_ped.py            # Converts FASTA → SNP data (PED/INFO files)
│   ├── run_matching.py            # Main workflow (runs everything)
│   ├── matching.py                # Helper functions for SNP analysis
│   ├── tree_utils.py              # Phylogenetic tree utilities
|   └── Haploview4.1.jar           # Linkage disequilibrium tool (included)
│
├── Data/
│   ├── InputFiles/                # PUT YOUR SEQUENCE FILES HERE
│   │   └── (empty — add your FASTA files)
│   │
│   ├── OutputFiles/               # Results appear here automatically
|   |   ├── IQTree_out/
│   │   └── (created when you run the tool)
│   │
│   └── requirements/
│       └── requirements.txt        # List of Python packages
│
└── README.md                       # Documentation
```

---

## How to Use

### The Basic Command

The simplest way to run the tool:

```bash
python Code/run_matching.py -a Data/InputFiles/your_alignment.fasta --haploview-jar Code/Haploview4.1.jar
```

**What this does:**
1. Reads your FASTA file from `Data/InputFiles/`
2. Extracts SNPs (variable positions)
3. Filters out low-quality data
4. Runs Haploview to compute linkage disequilibrium
5. Builds a phylogenetic tree (if IQ-TREE is installed)
6. Creates output files in `Data/OutputFiles/`

### Common Commands

**Generate only PED/INFO files (no Haploview/trees):**

```bash
python Code/run_matching.py -a Data/InputFiles/alignment.fasta --pedinfo-only
```

Use this if:
- You just want the SNP data in standard format
- You want to import into other tools (PLINK, etc.)
- You don't need the LD analysis or trees

**Skip IQ-TREE:**

```bash
python Code/run_matching.py -a Data/InputFiles/alignment.fasta \
  --haploview-jar Code/Haploview4.1.jar --skip-iqtree
```

Use this if:
- You just want the SNP data in standard format
- You want to import into other tools (PLINK, etc.)
- You don't need the LD analysis or trees

**Keep intermediate files for inspection:**

```bash
python Code/run_matching.py -a Data/InputFiles/alignment.fasta \
  --haploview-jar Code/Haploview4.1.jar --keep-intermediate
```

Use this if:
- You want to examine the PED/INFO files that Haploview uses
- You want to keep the LD file for later inspection
- You're debugging issues

**Full analysis with all options:**

```bash
python Code/run_matching.py -a Data/InputFiles/alignment.fasta \
  --haploview-jar Code/Haploview4.1.jar --keep-intermediate
```

### Important Notes

- **Filename matters:** Output folders are named after your input file. `alignment.fasta` → output in `Data/OutputFiles/alignment/`
- **Re-running:** If you run twice, the tool won't overwrite results. It creates `alignment_1/`, `alignment_2/`, etc.
- **Always activate the environment:** Before running any command, make sure you see `(.venv)` in your prompt. If not, run:
  - Windows: `.venv\Scripts\activate`
  - Mac/Linux: `source .venv/bin/activate`

---

## Understanding the Pipeline

### What Happens Step-by-Step

When you run the tool, it performs these operations on your FASTA file:

#### Step 1: Parse Input
- Reads your FASTA alignment file
- Checks that all sequences have the same length
- Counts total sequences and alignment length

#### Step 2: Detect Coverage Boundaries
- Looks at gaps (dashes `-`) at the start and end of each sequence
- Finds the region where **every sequence has at least one valid character**
- Trims off the incomplete parts
- **Example:** If one sequence starts with 100 gaps and another ends with 50 gaps, those regions are trimmed from all sequences

#### Step 3: Find Variable Positions
- Compares all sequences position by position
- Identifies positions where at least two different bases appear
- Ignores positions that are identical across all sequences (not informative)

#### Step 4: Filter Ambiguous Bases
- Removes positions containing invalid characters (anything outside `A, T, C, G, -`)
- Keeps analysis clean and interpretable

#### Step 5: Filter Incomplete Columns
- Removes positions where some sequences have `N` (unknown) or other non-standard bases
- Ensures high data quality

#### Step 6: Remove Singletons
- Removes positions where only ONE sequence differs from the majority
- These are usually sequencing errors, not real variants
- Reduces noise in the analysis

#### Step 7: Split Multi-Allelic Sites
- Some positions have 3+ different bases
- Each minority base gets its own column for statistical analysis

#### Step 8: Determine Consensus
- For each position, identifies the most common base (consensus)
- Used as the "reference" for encoding

#### Step 9: Encode Genotypes
- For each sequence at each position:
  - `1 1` = matches consensus base (same as most common)
  - `1 2` = differs from consensus (different base)
- This is standard population genetics encoding

#### Step 10: Write Output Files
- **PED file:** SNP data in standard linkage format
- **INFO file:** Position descriptions for Haploview

---

## Input and Output Files

### Input Format: FASTA Files

Your input should be a **multiple sequence alignment** in FASTA format.

**What FASTA looks like:**
```
>Sample_1
ATCGATCGATCGATCGATCGATCG
>Sample_2
ATCGATCGATCGATCGATCGATCG
>Sample_3
ATCGATCGTTCGATCGATCGATCG
```

**Rules:**
- Each sequence starts with `>` followed by a sample name
- The name can contain spaces, dates, ID numbers, etc.
- Sequences can be on one line or split across multiple lines
- All sequences MUST be the same length (aligned)
- Valid bases: `A`, `T`, `C`, `G`, `-` (gap), `N` (unknown)
- Case doesn't matter (ATCG or atcg both work)

**Example of a valid header:**
```
>hMpxV/DRC/HGRK-1L/2024|EPI_ISL_18886301|2024-01-15
```

### Output Files

After running, you'll find results in `Data/OutputFiles/alignment_name/`:

#### PED File (`*_ped.txt`)
- Standard linkage format used by genetics software
- Tab-separated columns
- First 6 columns: Family ID, Individual ID, Father ID, Mother ID, Sex, Phenotype
- Remaining columns: Genotype values for each SNP
- **Use this for:** Importing into PLINK, other population genetics tools

#### INFO File (`*_info.txt`)
- Two columns (both identical): SNP position
- One row per variable position
- **Use this for:** Reference with PED file, LD calculations

#### Haploview LD File (`*.ld`, if running Haploview)
- Linkage disequilibrium matrix
- Shows which SNPs are statistically linked
- **Use this for:** Understanding SNP associations

#### IQ-TREE Output (if IQ-TREE installed)
- `*.treefile` — The phylogenetic tree in Newick format
- `*.contree` — Consensus tree
- `*.log` — IQ-TREE analysis log
- **Use this for:** Visualizing evolutionary relationships

#### Excel Summary (if created)
- Human-readable statistics about your analysis
- Number of SNPs found, filtered, kept
- Tree statistics if available

### Understanding Filter Summary

When the tool runs, you'll see output like:

```
============================================================
  FASTA → PED/INFO  Filter Summary
============================================================
  Sequences in alignment         :     213
  Original alignment columns     :  198855
  — Coverage trim (removed)      :   19744
  After coverage filter          :  179111
  — Non-variable (removed)       :  165266
  After variable filter          :   13845
  — Ambiguous base (removed)     :   12819
  After ambiguity filter         :    1026
  — Incomplete col (removed)     :       0
  After completeness filter      :    1026
  — Singletons (removed)         :     559
  After singleton filter         :     467
  + Split multi-allelic          :       0
  Final virtual positions        :     467
============================================================
```

**What this tells you:**
- Started with 213 sequences, 198,855 positions
- Coverage trimming removed 19,744 incomplete columns
- Variable filter kept only positions with differences (removed 165,266 identical positions)
- Final result: 467 usable SNPs from 13,845 variable positions
- This is normal — most positions are either identical or unreliable

---

## Example Workflow

### Scenario: You Have 50 Monkeypox Sequences

1. **Prepare your file**
   - You have a FASTA file: `mpox_sequences.fasta`
   - Copy it to `Data/InputFiles/mpox_sequences.fasta`

2. **Activate the environment**
   ```bash
   # Windows
   .venv\Scripts\activate
   
   # Mac/Linux
   source .venv/bin/activate
   ```

3. **Run the pipeline**
   ```bash
   python Code/run_matching.py \
     -a Data/InputFiles/mpox_sequences.fasta \
     --haploview-jar Code/Haploview4.1.jar --keep-intermediate
   ```

4. **Wait for completion**
   - Takes seconds to minutes depending on file size
   - You'll see progress output on screen
   - Haploview window may pop up

5. **Check results**
   ```
   Data/OutputFiles/mpox_sequences/
   ├── mpox_sequences_ped.txt      # SNP genotypes
   ├── mpox_sequences_info.txt     # SNP positions
   ├── mpox_sequences.ld           # Linkage data
   ├── summary.xlsx                # Summary statistics
   └── IQTree_out/
       ├── mpox_sequences.treefile # Your tree
       └── (other IQ-TREE files)
   ```

6. **View results**
   - Open `summary.xlsx` in Excel or Google Sheets
   - Open `.treefile` in a tree viewer (FigTree, Dendroscope, etc.)
   - Use PED/INFO files with other tools

---

## Troubleshooting

### "Python command not found"

**Problem:** When you type `python`, you get "command not found" or similar error.

**Solution:**
- **Windows:** Python wasn't added to PATH during installation. Reinstall Python with "Add Python 3.x to PATH" checked ✓
- **Mac:** Use `python3` instead of `python` for all commands
- **Linux:** Use `python3` instead of `python`

### "Java not found" or Haploview won't run

**Problem:** Error mentions Java or Haploview when running the pipeline.

**Solution:**
1. Check Java is installed: `java -version`
2. If not installed, download from oracle.com or use:
   - Mac: `brew install openjdk`
   - Linux: `apt install default-jdk`
3. Try the command again

### "No module named 'pandas'" or similar package error

**Problem:** Error during execution about missing Python packages.

**Solution:**
1. Make sure the environment is activated (see `(.venv)` in prompt)
2. Reinstall packages:
   ```bash
   pip install -r Data/requirements/requirements.txt
   pip install tqdm
   ```
3. Try running again

### "FASTA validation error" or "Sequences differ in length"

**Problem:** Error reading your FASTA file.

**Solution:**
- Your sequences aren't aligned (different lengths)
- Check that all sequences have the same length
- Use a sequence alignment tool (Mafft, ClustalW, etc.) if needed
- Make sure there are no empty lines in the middle of sequences

### Output folder says "alignment_1", "alignment_2", etc.

**This is normal!** 
- The tool creates a new folder each run to avoid overwriting results
- If you run with the same file twice, results go to separate folders
- To clean up: delete old folders in `Data/OutputFiles/`

### Pipeline runs but output looks wrong

**Possible issues:**
1. **Too few SNPs found** — Your sequences might be very similar or from the same virus strain
2. **All SNPs filtered** — Your data might contain many ambiguous bases or gaps; check filter summary
3. **Tree not built** — IQ-TREE not installed; it's optional but needed for trees

**Quick fix:**
- Re-run with `--keep-intermediate` flag to inspect intermediate files
- Check that input FASTA is formatted correctly

### "Out of memory" error when Haploview runs

**Problem:** Haploview crashes with OutOfMemoryError.

**Solution:**
If you need to run Haploview manually with more memory:
```bash
java -Xmx4G -jar Code/Haploview4.1.jar
```

(Adjust `4G` to more/less memory as needed)

---

## FAQ

### Q: Can I use nucleotide sequences other than FASTA?

**A:** No, only FASTA format is supported. If you have a different format, use a converter tool like:
- SeqKit (online or command-line)
- Format-checking websites (NCBI, EBI)

### Q: What if some of my sequences are much shorter?

**A:** You need to align them first. Use alignment tools:
- **Mafft** — Fast and accurate
- **ClustalW** — Classic alignment tool
- **MUSCLE** — Good balance of speed and accuracy
- Web servers — NCBI Blast, EBI tools

### Q: Can I run this on Windows/Mac/Linux?

**A:** Yes, all three platforms are supported. The setup is slightly different but the tool works identically.

### Q: How large can my input file be?

**A:** Practical limits:
- **Good:** 50-1,000 sequences, 20,000-500,000 bp
- **May be slow:** 1,000-10,000 sequences or very long alignments (>1 million bp)
- **Very slow:** 10,000+ sequences or extreme alignments

For extremely large files, you might need:
- A computer with more RAM (8GB+)
- More time for processing (hours)
- Possibly splitting into smaller chunks

### Q: Do I need IQ-TREE?

**A:** No, it's optional. Without it:
- Pipeline still works
- You get PED/INFO files and LD analysis
- You don't get phylogenetic trees

### Q: Can I use RNA sequences?

**A:** Yes! The tool works with RNA (uses U instead of T), though it's designed for DNA. Just make sure:
- Sequences are aligned
- Format is valid FASTA
- All sequences use the same base encoding

### Q: What do the output genotype codes mean?

**A:** In the PED file:
- `1 1` = Homozygous reference (matches consensus base)
- `1 2` = Heterozygous (differs from consensus base)
- `0 0` = Missing data

This follows standard population genetics encoding used by PLINK and similar tools.

### Q: Can I import the output into PLINK or other tools?

**A:** Yes! The PED and INFO files are standard formats:
- **PLINK:** Can read PED files directly
- **Haploview:** Can read PED/INFO (already run by default)
- **Population genetics tools:** Most tools accept this format

### Q: How do I interpret the phylogenetic tree?

**A:** The tree shows evolutionary relationships:
- Branch length = evolutionary distance (more divergent = longer branch)
- Samples close together = similar sequences (small genetic distance)
- Samples far apart = different sequences (large genetic distance)
- Use FigTree or Dendroscope to visualize and customize trees

### Q: The tree looks wrong — what do I do?

**A:** Common issues:
1. **Very similar sequences** — Difficult to build accurate trees; normal
2. **Mixed quality** — Some good sequences, some poor quality
3. **Recombination** — Evolutionary trees may be misleading if recombination occurred
4. **Too few SNPs** — Trees need variation to build

Try:
- Inspecting the `.iqtree` log file for IQ-TREE statistics
- Using different evolution models (beyond default scope)

### Q: Can I run this from a USB drive or cloud storage?

**A:** Not recommended. Python environments work best on local disk. Instead:
1. Install on your local computer
2. Copy your input file to the project
3. Run normally
4. Copy results out when done

### Q: What if I want to modify the filtering parameters?

**A:** You'll need to edit the Python code:
1. Open `Code/fasta_to_ped.py` in a text editor
2. Look for filter thresholds (singleton removal, ambiguity checking, etc.)
3. Modify values and save
4. Re-run the pipeline

This requires some Python knowledge — consult documentation or ask for help.

---

## Getting Help

If you hit problems:

1. **Check this guide** — Read the Troubleshooting section above
2. **Examine filter summary** — Does it show expected SNP counts?
3. **Inspect intermediate files** — Run with `--keep-intermediate` to see raw PED/INFO
4. **Check log files** — Look in output folder for `.log` or `.iqtree` files
5. **Ask for help** — Open a GitHub issue or contact the developers

---

## Citation

If you use this tool in research, please cite:

```
Poyer, A. T, Feehley, M. C, & Feehley, P. J. (2026). Recombina: A toolkit for viral recombination and SNP matching analysis (Version 1.0.0) [Computer software]. GitHub. https://github.com/Hocket/Recombina
```

---

## Summary

You now know how to:
- ✓ Install the tool
- ✓ Prepare your input files
- ✓ Run the complete pipeline
- ✓ Understand what each step does
- ✓ Interpret the outputs
- ✓ Troubleshoot common issues

**Next steps:**
1. Install Python and Java (if not already done)
2. Clone/download this repository
3. Follow the Installation Guide (3 steps: setup environment, install packages, verify)
4. Prepare your FASTA file
5. Run your first analysis!

Good luck with your viral sequence analysis! 🧬