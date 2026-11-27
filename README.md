# Allelic diversity

This repository contains a set of Python scripts to **explore allelic diversity in large collections of bacterial isolates**, based on a core-genome or whole-genome MLST (cg/wgMLST) allelic matrix and a metadata table.

The scripts allow users to investigate the distribution of alleles across different groups of isolates (e.g. lineages), identify group-specific and shared alleles, and visualize allelic diversity in the context of a species phylogenetic tree.

## Scripts

### 1. `allelic_heatmap.py` (v2.1.0)

This script generates an **interactive heatmap of allelic diversity associated with a species phylogenetic tree**. Alleles are colored according to their distribution across user-defined groups (based on a selected metadata variable), allowing users to explore patterns of allelic conservation, sharing and group-specificity. The visualization can help identify unusual patterns of allele distribution across phylogenetic groups and highlight potential recombination events or horizontal gene transfer. In addition to the interactive heatmap, the script generates several complementary tables describing allele distribution, exclusivity and dominance across groups.

#### Input

- **Allelic matrix** (`*.tsv`): cg/wgMLST allelic matrix, with samples as rows and loci as columns.
- **Metadata table** (`*.tsv`): metadata associated with the samples, with sample identifiers in the first column.
- **Phylogenetic tree** (`*.nwk`): Newick tree containing the sample identifiers used to order the heatmap.

Optionally, users can also provide:

- **Color configuration file** (`*.tsv`): a table defining the colors assigned to each metadata category.
- **Samples of interest** (`*.tsv` or `*.txt`): a list of sample identifiers to investigate separately from the remaining collection.

#### Main features

**a. Allele classification and color assignment**

For each locus, alleles are classified according to their distribution across the categories defined by the selected metadata variable (e.g. lineage).

- **Exclusive alleles:** alleles found in only one category. These are represented using a lighter version of the corresponding category color.
- **Shared alleles:** alleles found in two or more categories. These are represented using the corresponding category color in a darker shade. The assigned color reflects the category in which the allele is dominant.
- **Conserved alleles:** loci with only one allele among the regular samples are represented in black. Additionally, users can define thresholds to highlight alleles that are widespread across multiple categories or dominant in several categories.
- **Missing alleles:** represented in white.

Allele dominance can be determined using either:

- **Frequency mode (default):** the proportion of isolates within each category carrying the allele.
- **Count mode:** the absolute number of isolates within each category carrying the allele.

**b. Tree visualization**

Two additional options are available to improve tree visualization:

- Display sample identifiers at the tree tips, using a user-defined font size.
- Display the horizontal branch-length scale, with a user-defined tick-label font size.
  
**c. Samples of interest**

Users can optionally provide a list of samples of interest using `--sample_of_interest`. These samples are displayed alongside the remaining collection in the heatmap but are **excluded from all category-based calculations**, including allele frequency, allele counts and dominant-category assignment. This allows users to investigate specific isolates independently of the background population. Alleles found exclusively among the selected samples of interest are highlighted in purple by default:

- Alleles found in multiple samples of interest are classified as `exclusive (group sample_of_interest)`.
- Alleles found in a single sample of interest are classified as `exclusive (this sample_of_interest)`.

The default purple color can be customized through the color configuration file.

**d. Customizable colors**

Colors can either be assigned automatically or specified by the user through a TSV configuration file. When colors are automatically generated, the script saves the corresponding configuration to `<output_prefix>_colors.tsv`, allowing users to reuse the same colors in subsequent analyses.

#### Output

The script generates the following files:

- **`<output_prefix>.html`** – Interactive HTML visualization containing the phylogenetic tree and allelic heatmap. Users can explore individual cells by hovering over them to retrieve sample, locus, allele and category information. The visualization also includes options to download the heatmap as PNG or SVG.

- **`<output_prefix>_allelic_detailed_matrix.tsv`** – Detailed allelic matrix containing the original allele information, annotated with the classification of each allele:
  - `exclusive`: allele found in only one category.
  - `shared`: allele found in multiple categories.
  - `conserved`: allele associated with a conserved locus.
  - `exclusive (group sample_of_interest)` or `exclusive (this sample_of_interest)`: alleles exclusive to the selected samples of interest.
  - `missing`: missing allele information.

- **`<output_prefix>_allele_summary.tsv`** – Summary table listing the alleles classified as exclusive or shared across the different categories, including their corresponding loci. Alleles exclusive to samples of interest are also reported.

- **`<output_prefix>_allele_summary_counts.tsv`** – Summary of the number of distinct loci and alleles classified as exclusive or shared across the different categories, including samples of interest when applicable.

- **`<output_prefix>_allele_dominance.tsv`** – Table describing the dominance of each allele across the different categories. For each locus and allele, it reports the dominant category, allele classification, absolute counts and relative frequencies in each category.

- **`<output_prefix>_colors.tsv`** – Color configuration file generated automatically when no custom color file is provided. It includes the colors assigned to the metadata categories and the special `sample_of_interest` category.

- **`<output_prefix>.log`** – Log file containing information about the execution, including the parameters used, execution timestamps and generated output files.

#### Usage

```bash
python allelic_heatmap.py -h
```

Available arguments:

```text
  -h, --help            show this help message and exit
  -a ALLELIC_MATRIX, --allelic_matrix ALLELIC_MATRIX
                        Allelic matrix in TSV format.
  -m METADATA, --metadata METADATA
                        Metadata table in TSV format.
  -c COLORS, --colors COLORS
                        (Optional) TSV file indicating the colors to use in the heatmap. This file should contain the following columns in the header:
                        'category' (indicating the different categories observed in the specified metadata column) and color (indicating the color code to
                        use for the category). If '--sample_of_interest' is provided, this file should also include a row where 'category' corresponds to
                        the expression 'sample_of_interst' and 'color' should indicate the color code for alleles exclusively observed in samples of
                        interest. If this file is not provided colors are auto-assigned and saved to <prefix>_colors.tsv.
  -hcol HEADER, --header HEADER
                        Metadata column to use for categories.
  -t TREE, --tree TREE  Newick tree file.
  -o OUTPUT, --output_prefix OUTPUT
                        Output prefix.
  -mode {frequency,count}, --mode {frequency,count}
                        Dominant category mode: 'frequency' (default) or 'count'.
  -ncat NCAT, --ncat NCAT
                        Optional: if an allele is present in at least N categories, color it grey (nearly conserved).
  -ncat_dominant NCAT_DOMINANT, --ncat_dominant NCAT_DOMINANT
                        Optional: if an allele is dominant in at least N categories, color it grey (nearly conserved).
  --sample_of_interest SAMPLE_OF_INTEREST
                        [OPTIONAL] TSV or TXT file with the list of samples of interest. Sample IDs must be in the first column and no headers should be
                        provided.
  --tree-labels SIZE    [OPTIONAL] Show sample names at the tree tips. SIZE is the font size in points (e.g. 6).
  --tree-scale SIZE     [OPTIONAL] Show the horizontal tree branch-length scale at the bottom. SIZE is the tick-label font size in points (e.g. 8).
```

**Important:** When providing a custom color file together with `--sample_of_interest`, the color file must also contain a row with `category = sample_of_interest` for the selected metadata column.

---

### 2. `allelic_distribution.py`

This script provides a complementary approach to exploring allelic diversity by quantifying the proportion of isolates within each group that carry a given allele. Unlike `allelic_heatmap.py`, this script focuses on the **quantitative distribution of alleles across user-defined groups**, allowing users to identify differences in allele prevalence between lineages or other metadata categories.

#### Input

- **Allelic matrix** (`*.tsv`): cg/wgMLST allelic matrix, with samples as rows and loci as columns.
- **Metadata table** (`*.tsv`): metadata associated with the samples, with sample identifiers in the first column.

#### Output

- **`<output_prefix>_summary.tsv`** – Table in which each row corresponds to an allele and each column to a group of samples (defined by the selected metadata variable). Values indicate the proportion of isolates within each group carrying the corresponding allele.
- **`<output_prefix>_distribution.tsv`** – Table in which each row corresponds to a locus and each column to a group, describing the allelic distribution within each group.

#### Usage

```bash
python allelic_distribution.py -h
```

Available arguments:

```text
-h, --help
    Show this help message and exit.

-m METADATA, --metadata METADATA
    [MANDATORY] Metadata file in TSV format.

-a ALLELES, --allele ALLELES
    [MANDATORY] Allelic matrix in TSV format.

-c GROUP_COLUMN, --group-column GROUP_COLUMN
    [MANDATORY] Name of the metadata column containing
    the groups of interest.

-g GROUP_INTEREST, --group-interest GROUP_INTEREST
    [OPTIONAL] Comma-separated list of groups to select.
    If omitted, all groups are considered.

-o OUTPUT, --output OUTPUT
    [OPTIONAL] Prefix for output files.
    Default: Allele_distribution.
```

---

## Input file requirements

For both scripts, the following requirements should be met:
- Sample identifiers must match between the allelic matrix and metadata table.
- Sample identifiers in the phylogenetic tree must match those in the allelic matrix.
- The first column of both the allelic matrix and metadata table must contain sample identifiers.
- Missing alleles should be represented by `0`.
- The samples-of-interest file must contain sample identifiers in its first column, without a header.

---

## Command-line examples

Example input files are available in the [Zenodo repository](https://doi.org/10.5281/zenodo.17177579). We recommend downloading and unzipping the corresponding folders before running the examples below.

The following examples use the *Chlamydia trachomatis* dataset provided in that repository.

### 1. Generate an allelic heatmap ordered according to a species phylogenetic tree

The following command generates an interactive heatmap in which alleles are colored according to the lineage in which they are dominant, determined using allele frequency. Alleles present in at least three lineages are highlighted in black.

```bash
python allelic_heatmap.py \
    -a 02_input_to_setup_cgMLST_1230_background/Ct_alleles_missing_code_0.tsv \
    -m 02_input_to_setup_cgMLST_1230_background/Ct_metadata.tsv \
    -t 03_example_ReporTree_output_1230_genomes/HC/Ct_HC_single_HC.nwk \
    -hcol lineage \
    -o TEST_allelic_heatmap \
    -mode frequency \
    -ncat 3
```

### 2. Generate an allelic heatmap highlighting dominant alleles

Alternatively, users can highlight alleles that are dominant in at least three different lineages.

```bash
python allelic_heatmap.py \
    -a 02_input_to_setup_cgMLST_1230_background/Ct_alleles_missing_code_0.tsv \
    -m 02_input_to_setup_cgMLST_1230_background/Ct_metadata.tsv \
    -t 03_example_ReporTree_output_1230_genomes/HC/Ct_HC_single_HC.nwk \
    -hcol lineage \
    -o TEST_allelic_heatmap \
    -mode frequency \
    -ncat_dominant 3
```

### 3. Investigate samples of interest

To investigate specific samples independently of the remaining collection, users can provide a text file containing the sample identifiers of interest.

For example, create a file named `samples_of_interest.txt`:

```text
Sample_001
Sample_002
Sample_003
```

The following command generates a heatmap in which these samples are highlighted and excluded from category-based calculations. Alleles found exclusively in these samples are represented in purple.

```bash
python allelic_heatmap.py \
    -a 02_input_to_setup_cgMLST_1230_background/Ct_alleles_missing_code_0.tsv \
    -m 02_input_to_setup_cgMLST_1230_background/Ct_metadata.tsv \
    -t 03_example_ReporTree_output_1230_genomes/HC/Ct_HC_single_HC.nwk \
    -hcol lineage \
    -o TEST_allelic_heatmap_SOI \
    -mode frequency \
    -ncat_dominant 3 \
    --sample_of_interest samples_of_interest.txt \
    --tree-labels 6 \
    --tree-scale 8
```

### 4. Obtain the allelic distribution per lineage

The following command generates tables describing the proportion of isolates carrying each allele across the different lineages.

```bash
python allelic_distribution.py \
    -a 02_input_to_setup_cgMLST_1230_background/Ct_alleles_missing_code_0.tsv \
    -m 02_input_to_setup_cgMLST_1230_background/Ct_metadata.tsv \
    -c lineage \
    -o TEST_allelic_distribution
```

## Installation

These scripts are written in Python 3.10. To facilitate installation and ensure that the required dependencies are available, a Conda environment configuration file (`allelic_diversity.yml`) is provided.

To install and activate the environment:

```bash
cd /PATH/TO/allelic_diversity/

conda env create -f allelic_diversity.yml

conda activate allelic_diversity
```

To check that the scripts are correctly installed:

```bash
python allelic_heatmap.py -h

python allelic_distribution.py -h
```

**Note:** The interactive HTML visualization uses Plotly's CDN to load the JavaScript library. Therefore, an internet connection is required when opening the HTML file.

---

## Citation

If you use any of these scripts, please cite:

[Lodhia, Z., Mixão, V., Isidro, J. et al. A new core-genome MLST approach for Chlamydia trachomatis genomic surveillance and research. Nat Commun 17, 9228 (2026). https://doi.org/10.1038/s41467-026-77098-8](https://www.nature.com/articles/s41467-026-77098-8)

---

## Funding

This work was funded by national funds through FCT – Foundation for Science and Technology, I.P., in the frame of Individual CEEC 2022.00851.CEECIND/CP1748/CT0001 ([doi: 10.54499/2022.00851.CEECIND/CP1748/CT0001](https://doi.org/10.54499/2022.00851.CEECIND/CP1748/CT0001)) and of the doctoral fellowship SFRH/BD/147446/2019 ([doi: 10.54499/SFRH/BD/147446/2019](https://doi.org/10.54499/SFRH/BD/147446/2019)), and by the European Union project "Sustainable use and integration of enhanced infrastructure into routine genome-based surveillance and outbreak investigation activities in Portugal" – GENEO [101113460], on behalf of the EU4H programme [EU4H-2022-DGA-MS-IBA-01-02].
