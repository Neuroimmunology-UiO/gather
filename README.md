
[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green)](LICENSE)
[![Bioconda](https://img.shields.io/badge/bioconda-available-brightgreen)](https://bioconda.github.io/)
[![Singularity Image](https://img.shields.io/badge/singularity-image-ff69b4)](https://github.com/Neuroimmunology-UiO/gather/releases)



# 🕷️ GATHeR: Graph-based Accurate Tool for immunoglobulin HEavy- and light-chain Reconstruction

This repository provides the environment setup for `gather`, including all necessary dependencies.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Processing 10x Genomics Chromium Single-Cell RNA-Seq Data for Assembly](#processing-10x-genomics-chromium-single-cell-rna-seq-data-for-assembly)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Quick Setup (Linux, via Singularity – Recommended)

To simplify installation and ensure reproducibility, we provide a prebuilt Singularity image: `GATHER.sif`.

**Before you begin:**  
Ensure [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) is installed on your system. You can check with:

```bash
singularity --version
```
### Download the Image

The `GATHER.sif` image is available for direct download from the [Releases section](https://github.com/Neuroimmunology-UiO/gather/releases) of our GitHub repository.

Once downloaded, you can immediately run any GATHeR tool from within the container. For example:

```bash
singularity exec GATHER.sif sc_asm.py --help
```

### Quick Setup (macOS & Linux, via Docker – Recommended for Mac)

To simplify installation and ensure reproducibility, we also provide a **Docker image**.  

**Before you begin:**  
Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is installed on your system. You can check with:

```bash
docker --version
```

Now, you can immediately run any GATHeR tool using the following Docker command as an example:

```bash
docker run --rm -it \
  -u "$(id -u):$(id -g)" \
  -v /path/to/your/scRNAseq_folder:/data \
  --entrypoint bash \
  gather:1.0 \
  -lc "sc_asm.py \
        --seq_1 /data/<read1.fastq> \
        --seq_2 /data/<read2.fastq> \
        --output_dir /data/output \
        --num_jobs 4"
```
### Conda-Based Installation (Alternative to Singularity)

If you prefer a local development environment or wish to integrate GATHeR into an existing Conda workflow, you can install it via Conda using [Bioconda](https://bioconda.github.io/). This method is platform-independent and works well for users running Linux or macOS with Conda installed.

We recommend creating a dedicated environment for GATHeR:

```bash
conda create --name gather-env python=3.8
conda activate gather-env
```

Once the environment is active, install the `gather` package and its core dependencies from Bioconda and Conda-Forge:

```bash
conda install -c bioconda gather
```

### Dependencies

In addition to the Python-based modules installed with Conda, GATHeR relies on several **external tools** and **R packages** that are required for assembly, annotation, and phylogenetic analysis.  

#### External Tools
- **NCBI BLAST+** — required for sequence similarity searches  
- **SPAdes** — de novo RNA assembly
- **IgPhyML** — phylogenetic inference of B-cell receptor lineages  
- **Immcantation helper scripts** — utilities for configuring IgBLAST/IMGT references  

#### R Packages
- **BiocManager** — package manager for Bioconductor  
- **treeio**, **ggtree**, **dowser** — Bioconductor packages for lineage tree analysis  
- **dplyr**, **ggrepel** — CRAN packages for data wrangling and visualization  

To install these dependencies in one step, we provide an automated setup script:  

📄 [`setup_gather_dependencies.sh`](https://github.com/Neuroimmunology-UiO/gather/blob/main/scripts/setup_dependencies.sh)

```bash
# Make the script executable
chmod +x setup_gather_dependencies.sh

# Run the installer
./setup_gather_dependencies.sh
```


### Prepare IgBLAST and Reference Databases

To assign V(D)J genes and annotate the junctional regions, we rely on [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/) and reference data from the IMGT database. Setting up the IgBLAST environment requires a few one-time setup steps.

For convenience, we provide an automated setup script based on the [Change-O IgBLAST setup guide](https://changeo.readthedocs.io/en/stable/examples/igblast.html). The script handles download, extraction, and configuration of the necessary files and tools.

You can download the setup script from our GitHub repository:

📄 [`setup_igblast_env.sh`](https://github.com/Neuroimmunology-UiO/gather/blob/main/scripts/setup_igblast_env.sh)

After downloading, make the script executable and run it:

```bash
chmod +x setup_igblast_env.sh
./setup_igblast_env.sh
```

> **Note:** You can adjust the IgBLAST version or installation directories inside the script if needed. After running it, restart your shell or run `source ~/.bashrc` to ensure the `IGDATA` environment variable is correctly set.

For more details and troubleshooting, refer to the official [Change-O IgBLAST setup guide](https://changeo.readthedocs.io/en/stable/examples/igblast.html).

## Usage

GATHeR supports single-cell RNA sequencing data from technologies like Smart-seq2/3 and 10x Genomics Chromium.

### Merging Paired-End Reads

If your reads are paired-end:

```bash
sc_asm.py --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz --output_dir .
```

### Single-End Data

```bash
sc_asm.py --seq_1 {cell_file_name}_R1.fastq.gz --output_dir .
```
### Output Files

- `{cell_file_name}_merged.unitigs.fa`: Unitigs from the cDBG graph  
- `{cell_file_name}_merged.contigs.fa`: Assembled contigs (transcriptome, GATHeR algorithm)  
- `{cell_file_name}_merged.BCR.fa`: Annotated BCR sequences (SPAdes and GATHeR algorithms denoted as algo_1 and algo_2 in the headers)
  
### K-mer Size

- Default: `k = 25`
- Rule of thumb: `k < read length - 20`
- For short reads (e.g., 30 bp), the feasible k-mer size range is constrained; values between 21 and 25 are generally appropriate, with 27 being a possible upper limit depending on the context.
- 
## Processing 10x Genomics Chromium Single-Cell RNA-Seq Data for Assembly

To assemble data generated by 10x Genomics Chromium single-cell RNA-seq, the raw sequencing reads must first be processed using Cell Ranger, which organizes the output files in the following structure:

```
project_folder/
├── fastq/
│   ├── sample_name_S1_L001_R1_001.fastq.gz  # Read 1: contains cell barcodes and UMIs
│   ├── sample_name_S1_L001_R2_001.fastq.gz  # Read 2: contains the cDNA read
│   └── sample_name_S1_L001_I1_001.fastq.gz  # Index read (optional; used for sample demultiplexing)
```

The list of detected cell barcodes is typically located in:

```
outs/filtered_feature_bc_matrix/barcodes.tsv.gz
```

### Prepare the Data for Assembly

1. Concatenate all R1 reads:

    ```bash
    cat *R1*.fastq.gz > Merged_R1.fastq.gz
    ```

2. Concatenate all R2 reads:

    ```bash
    cat *R2*.fastq.gz > Merged_R2.fastq.gz
    ```

3. Decompress the barcode file:

    ```bash
    gunzip barcodes.tsv.gz
    ```

4. Run the assembly script:

    ```bash
    10x_asm.py --barcode_file barcodes.tsv \
               --r1_fastq Merged_R1.fastq.gz \
               --r2_fastq Merged_R2.fastq.gz \
               --output_dir_name cells
    ```

**Note**: The `--num_jobs` parameter defaults to 8, but for optimal performance, set it to the maximum number of available CPU cores.

## Clonality analysis and constant region analysis

After assembling B-cell receptor (BCR) sequences for each cell, the reconstructed heavy and light chains—saved in separate FASTA files with distinct headers—can be collected across all cells for clonality and constant region analysis.
To prepare the input, copy all `{cell_file_name}_merged.BCR.fa` files from the individual cell directories (if they are not already in a shared location) into a single directory, e.g. `BCRs_DIR`.Then, you can run the `postproc.py` script to perform V(D)J gene assignment, CDR3 parsing, productivity assessment, and, optionally, clonal grouping and lineage reconstruction.

### Step 1: Analyze Light Chains

For light chains, use the following command:

```bash
postproc.py --bcrs_dir <PATH>/BCRs_DIR --chain light --output_dir .
```

This will produce:

- `light_chains_db-pass.tsv`: a **Change-O–like tab-delimited file** containing detailed annotations of the light chain sequences, including gene calls, junction regions, productivity, and alignment features.

### Step 2: Analyze Heavy Chains with Clonal Inference and Visualization

To get only `heavy_chains_db-pass.tsv`: a **Change-O–like tab-delimited file** including  **constant region polymorphism analysis**:

```bash
postproc.py --bcrs_dir <PATH>/BCRs_DIR --chain heavy --output_dir .
```
#### For heavy chains, you can include additional flags to enable clonal grouping and downstream visualization:

```bash
postproc.py --bcrs_dir <PATH>/BCRs_DIR \
                      --chain heavy \
                      --clonality \
                      --lineage_tree \
                      --output_dir .
```
This command will generate:
- `heavy_chains_db-best.tsv` /  `light_chains_db-best.tsv`: a tab-delimited file with full BCR annotations using only the best assembled sequence for each cell.
- `heavy_chains.fmt7` / `light_chains.fmt7`: Raw IgBLAST output in format 7 (tabular/structured text).
- `heavy_chains_db-best_clone-pass.tsv`: a tab-delimited file with full BCR annotations **plus clonal assignments**.
- A BCR **clonal network plot**, where sequences are visualized as nodes connected by lines if they belong to the same clone (`--clonality`).
- Phylogenetic **lineage trees** of individual clones reconstructed from inferred germline sequences (`--lineage_tree`). This step can also be performed with [Dowser](https://github.com/immcantation/dowser?tab=readme-ov-file) using the GATHER output tables.

## License

Licensed under the Apache License. See the [LICENSE](LICENSE) file for details.

