
[![Python 3.12+](https://img.shields.io/badge/python-3.12%2B-blue)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-Apache%202.0-green)](LICENSE)
[![Bioconda](https://img.shields.io/badge/bioconda-available-brightgreen)](https://bioconda.github.io/)
[![Singularity Image](https://img.shields.io/badge/singularity-image-ff69b4)](https://github.com/Neuroimmunology-UiO/gather/releases)

# 🕷️ GATHeR: Graph-based Accurate Tool for immunoglobulin HEavy- and light-chain Reconstruction

This repository provides the environment setup for `gather`, including all necessary dependencies. Follow the instructions below to create and activate the environment using `conda`.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Processing 10x Genomics Chromium Single-Cell RNA-Seq Data for Assembly](#processing-10x-genomics-chromium-single-cell-rna-seq-data-for-assembly)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Quick Setup with Singularity (Recommended)

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

In addition to the Python-based tools installed above, GATHeR relies on several external tools and libraries. These include:

<!-- You'll fill in this section next -->

---

Let me know what you’d like to include in the `Dependencies` subsection — for example, if it should list tools like `bcalm`, `spades`, `blast`, or anything else that’s critical for the pipeline.


### Step 1: Prepare IgBLAST and Reference Databases

To assign V(D)J genes and annotate the junctional regions, we rely on [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/) and reference data from the IMGT database. Setting up the IgBLAST environment requires a few one-time setup steps.

The setup process is supported by a set of helper scripts available from the [Immcantation project](https://immcantation.org/), which can be found in their [GitHub repository](https://github.com/immcantation/immcantation/tree/master/scripts). These scripts simplify the configuration of IgBLAST and IMGT references. You will need the following setup scripts:

- `fetch_igblastdb.sh` — downloads core IgBLAST databases
- `fetch_imgtdb.sh` — downloads IMGT reference sequences
- `clean_imgtdb.py` — cleans downloaded IMGT sequences (optional for some datasets)
- `imgt2igblast.sh` — converts IMGT references to IgBLAST-compatible format

To use these scripts, copy all the tools from the `/scripts` folder in the [Immcantation repository](https://github.com/immcantation/immcantation/tree/master/scripts) into a directory in your `PATH`.

### Step 1: Prepare IgBLAST and Reference Databases

To perform V(D)J annotation and clonality analysis, IgBLAST must be configured with IMGT reference sequences. This setup includes downloading IgBLAST binaries, fetching reference databases, and converting IMGT germline files into a format compatible with Change-O.

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
cat {cell_file_name}_R1.fastq.gz {cell_file_name}_R2.fastq.gz > {cell_file_name}_merged.fastq.gz
```

Then run:

```bash
sc_asm.py --seq_merged {cell_file_name}_merged.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz
```

### Single-End Data

```bash
sc_asm.py --seq_merged {cell_file_name}_R1.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz
```
### Output Files

- `{cell_file_name}_merged.unitigs.fa`: Unitigs from the cDBG graph  
- `{cell_file_name}_merged.contigs.fa`: Assembled contigs (transcriptome, GATHeR algorithm)  
- `{cell_file_name}_merged.BCR_algo1.fa `: Annotated BCR sequences (GATHeR algorithm)
- `{cell_file_name}_merged.BCR_algo2.fa `: Annotated BCR sequences (SPAdes algorithm)
- `{cell_file_name}_merged.BCR_contiguous.fa`: Filtered BCR sequences based on the weight and contiguity using the contigs from the both algos.

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

## Clonality analysis

After successfully assembling B-cell receptor (BCR) sequences for each cell, the reconstructed heavy and light chains—each saved in separate FASTA files with distinct headers—can be collected across all cells (`Collected_heavy_chains.fasta` and `Collected_light_chains.fasta`) to perform clonality analysis.

Once merged, you can run the `clonality_analysis.py` script to perform V(D)J gene assignment, CDR3 parsing, productivity assessment, and, optionally, clonal grouping and lineage reconstruction.

### Step 1: Analyze Light Chains

For light chains, use the following command:

```bash
clonality_analysis.py --fasta_file Collected_light_chains.fasta --chain light
```

This will produce:

- `light_chains_db-pass.tsv`: a **Change-O–like tab-delimited file** containing detailed annotations of the light chain sequences, including gene calls, junction regions, productivity, and alignment features.

### Step 2: Analyze Heavy Chains with Clonal Inference and Visualization

For heavy chains, you can include additional flags to enable clonal grouping and downstream visualization:

```bash
clonality_analysis.py --fasta_file Collected_heavy_chains.fasta \
                      --chain heavy \
                      --clonal_plot \
                      --run_combined \
                      --lineage_tree
```

This command will generate:

- `heavy_chains_db-pass_clone-pass.tsv`: a tab-delimited file with full BCR annotations **plus clonal assignments** and **constant region polymorphism analysis**.
- A BCR **clonal network plot**, where sequences are visualized as nodes connected by lines if they belong to the same clone (`--clonal_plot`).
- Phylogenetic **lineage trees** of individual clones reconstructed from inferred germline sequences (`--lineage_tree`).

These outputs follow the standard formats used in the [Change-O](https://changeo.readthedocs.io/) framework and are compatible with downstream tools for repertoire diversity, lineage evolution, and visualization.


## License

Licensed under the Apache License. See the [LICENSE](LICENSE) file for details.

