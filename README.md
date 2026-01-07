# Oropouche SNPs Landscape (L / M / S)

[![R](https://img.shields.io/badge/R-%E2%89%A5%204.1-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#license)

An **R** pipeline to catalogue SNPs (with coding impact) and generate a *publication-ready* **mutational landscape** figure for Oropouche segments **L, M, and S** — now provided as a **single script** with **run flags** (Part 1 / Part 2) and **automatic caching** (skip Part 1 if outputs already exist).

---

## Table of contents
- [Repository layout](#repository-layout)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [How to run](#how-to-run)
  - [Flags: run Part 1 / Part 2](#flags-run-part-1--part-2)
  - [Automatic skip (cache) for Part 1](#automatic-skip-cache-for-part-1)
- [Pipeline overview](#pipeline-overview)
- [Plot visual customizations](#plot-visual-customizations)
  - [Colours](#colours)
  - [Export (PNG/SVG/PDF)](#export-pngsvgpdf)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [How to cite](#how-to-cite)

---

## Repository layout

Recommended structure:

```
.
├── scripts/
│   ├── Mapping-and-plot_SNPs_landscape.R         
│   └── Mapping-and-plot_SNPs_landscape_pt-br.R   
├── example_data/                    
│   ├── L_aln.fasta
│   ├── M_aln.fasta
│   └── S_aln.fasta
└── outputs/                         
    ├── L_snps_catalogue.csv
    ├── M_snps_catalogue.csv
    ├── S_snps_catalogue.csv
    ├── LMS_snps_catalogue.csv
    ├── LMS_snps_catalogue_summary.csv
    └── figures/
        ├── mutational_landscape_L_M_S.png
        ├── mutational_landscape_L_M_S.svg
        └── mutational_landscape_L_M_S.pdf
```

---

## Inputs

### Required
- Nucleotide **aligned FASTA (MSA)**, one per segment (L/M/S).
  - Must contain **one reference** and **N samples** (e.g., `Seq_01`, `Seq_02`, …).
  - All sequences must have the **same length** (valid alignment).
- These FASTA files are read from **`PARAM$input_dir`** (Part 1).

---

## Outputs

### Catalogue (Part 1)
- `L_snps_catalogue.csv`, `M_snps_catalogue.csv`, `S_snps_catalogue.csv`
- `LMS_snps_catalogue.csv` (merged L+M+S)
- `LMS_snps_catalogue_summary.csv` (summary by sample/segment/effect)

### Figure (Part 2)
- `mutational_landscape_L_M_S.png`
- `mutational_landscape_L_M_S.svg`
- `mutational_landscape_L_M_S.pdf`

---

## How to run

This pipeline is intended for **interactive execution in RStudio**.

### Step 0) Open the script
In the RStudio **Files** panel, open:

- `scripts/Mapping-and-plot_SNPs_landscape.R`

### Step 1) Edit parameters (BLOCO 1)
In the script, locate and edit the parameter blocks:

- `PARAM` (Part 1: FASTA inputs, reference headers, output folder)
- `PARAM_PLOT` (Part 2: plot inputs/outputs and figure export options)

At minimum, confirm:
- `PARAM$input_dir` (where your aligned FASTAs are located)
- `PARAM$output_dir` (where CSVs and Figures will be written; also the default input for Part 2)
- `PARAM$segments$L/M/S$fasta_aln` and `PARAM$segments$L/M/S$ref_id`

### Step 2) Run in RStudio
Choose one:
- Select the **entire script** and click **Run** (or `Ctrl+A` → `Ctrl+Enter`);
- Run **block-by-block** (recommended): select each block and press `Ctrl+Enter`.

---


### Note on Part 2 input
By default, Part 2 reads `LMS_snps_catalogue.csv` from **`PARAM$output_dir`** (because it is produced by Part 1).  
If you want to plot a catalogue generated elsewhere, set `PARAM_PLOT$workdir` and/or `PARAM_PLOT$arquivo_entrada` accordingly.

## Flags: run Part 1 / Part 2

At the top of the script (BLOCO 0), set:

```r
RUN_PART1   <- TRUE   # generate SNP catalogue (CSVs)
RUN_PART2   <- TRUE   # generate the mutational landscape plot
FORCE_PART1 <- FALSE  # recalc Part 1 even if LMS_snps_catalogue.csv already exists
```

Common use cases:

- **Run everything (Part 1 + Part 2)**  
  `RUN_PART1 <- TRUE` and `RUN_PART2 <- TRUE`

- **Run only Part 1 (generate/refresh CSVs)**  
  `RUN_PART1 <- TRUE` and `RUN_PART2 <- FALSE`

- **Run only Part 2 (plot), without recalculating**  
  `RUN_PART1 <- FALSE` and `RUN_PART2 <- TRUE`  
  (requires `LMS_snps_catalogue.csv` to exist)

---

## Automatic skip (cache) for Part 1

Part 1 checks whether the merged catalogue exists in **`PARAM$output_dir`**:

- `file.path(PARAM$output_dir, PARAM$merged_out_csv)`  
  (default: `.../LMS_snps_catalogue.csv`)

If it exists and `FORCE_PART1 <- FALSE`, **Part 1 is skipped automatically** and the pipeline proceeds to Part 2 (if enabled).

To force recalculation (e.g., after changing FASTAs/ref_id/filters), set:

- `FORCE_PART1 <- TRUE`

---

## Pipeline overview

1. Load the aligned FASTA and identify the reference sequence.
2. Scan alignment columns and call SNPs vs the reference (ignoring indels/gaps if configured).
3. Convert alignment coordinates → reference coordinates (gap-free reference positions).
4. Translate and classify mutation effects:
   - `synonymous` → **Syn**
   - `nonsynonymous`, `stop_gained`, `stop_lost` → **Nonsyn**
   - otherwise → **Unk**
5. Merge L/M/S and produce the final landscape plot faceted by sample/Node.

![Mutational landscape (L/M/S)](outputs/figures/mutational_landscape_L_M_S.png)

---

# Plot visual customizations

The sections below describe what to edit **inside the Part 2 block** of the combined script.

## Colours

### Mutation type (fixed)
- Syn: `#1f78b4` (blue)
- Nonsyn: `#e31a1c` (red)
- Unk: `grey50`

```r
scale_color_manual(
  name   = "SNPs landscape of L, M and S segments",
  breaks = c("Syn", "Nonsyn", "Unk"),
  values = c(Syn = "#1f78b4", Nonsyn = "#e31a1c", Unk = "grey50"),
  labels = c(
    Syn    = "Synonymous",
    Nonsyn = "Non-synonymous (incl. stop)",
    Unk    = "Unknown/ambiguous"
  )
)
```

---

## Export (PNG/SVG/PDF)

Recommendation:
- **PNG**: 300–600 dpi
- **SVG/PDF**: vector formats (ideal for Inkscape/Illustrator)
- Overlapping labels may happen; open the SVG in Inkscape and adjust if needed.

The script exports by default to:

- `file.path(PARAM_PLOT$out_dir_fig, "mutational_landscape_L_M_S.*")`

---

## Troubleshooting

### Reference not found
Ensure `ref_id` **exactly matches** the reference header in the aligned FASTA located in `PARAM$input_dir`.

### FASTA is not an alignment
All sequences must have the **same length**.

### Part 2: input CSV not found
If you run Part 2 only (`RUN_PART1 <- FALSE`), ensure the file exists at:

- `file.path(PARAM_PLOT$workdir, PARAM_PLOT$arquivo_entrada)`

---

## License

Recommended: **MIT** (permissive and commonly used in research repositories).  
Add a `LICENSE` file at the repository root.

---

## How to cite

Suggestion:

> Thiago Sousa. Oropouche SNP plotting pipeline (L/M/S) (vX.Y.Z) [Computer software]. Orcid: orcid.org/0000-0001-9809-8883. Accessed: YYYY-MM-DD.  
> Thiagojsousa. Oropouche SNP plotting pipeline (L/M/S) [Computer software]. GitHub: Thiagojsousa/Oropouche_Mutational_Landscape_SNPs.git. Accessed: YYYY-MM-DD. Website: thiagojsousa.com.br.
