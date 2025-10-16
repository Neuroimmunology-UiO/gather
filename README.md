[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![Bioconda](https://img.shields.io/badge/bioconda-available-brightgreen)](https://bioconda.github.io/)
[![Singularity Image](https://img.shields.io/badge/singularity-image-ff69b4)](https://github.com/Neuroimmunology-UiO/gather/releases)




# 🕷️ GATHeR: Graph-based Accurate Tool for Immunoglobulin HEavy- and Light-chain Reconstruction

GATHeR assembles BCR heavy and light chains from single-cell RNA-seq using de Bruijn graphs. It reconstructs paired BCR sequences and extends them into constant regions for accurate isoform, subclass, and allele assignment. Supports Smart-seq2/3 and 10x Genomics data.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Processing 10x Genomics Chromium 5′ GEX Data for Assembly](#processing-10x-genomics-chromium-5-gex-data-for-assembly)
- [Citation](#-citation)

## Installation

### Quick Setup (Linux, via Singularity — Recommended)

To simplify installation and ensure reproducibility, we provide a prebuilt Singularity image: `GATHER.sif`, which has been tested on Ubuntu 24.04.1.

**Before you begin**  
Ensure [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) is installed on your system:
```bash
singularity --version
```
### Download the Image

The `GATHER.sif` image is available for direct download from the [Releases section](https://github.com/Neuroimmunology-UiO/gather/releases) of our GitHub repository.

Once downloaded, you can immediately run any GATHeR tool from within the container. For example:

```bash
singularity exec GATHER.sif sc_asm.py --help
```
---
### Quick Setup (macOS & Linux, via Docker – Recommended for Mac)

We also provide a Docker image, which has been tested on MAC OS XX and Ubuntu 24.04.1.  

**Before you begin:**  
Ensure [Docker Desktop](https://www.docker.com/products/docker-desktop/) is installed on your system:

```bash
docker --version
```
Pull the Docker image and run:

```bash
docker pull seyedmos/gather:1.0
```

```bash
docker run --rm -it \
-u "$(id -u):$(id -g)" \
-v /path/to/your/scRNAseq_folder:/data \
--entrypoint bash \
seyedmos/gather:1.0 \
-lc "sc_asm.py \
--seq_1 /data/<read1.fastq> \
--seq_2 /data/<read2.fastq> \
--output_dir /data/output \
--num_jobs 4"
```
---

### Conda-based Installation (Alternative to Singularity)

If you prefer a local setup or want to integrate GATHeR into an existing Conda workflow, you can install it via Conda. Installation (tested on Ubuntu 24.04.1 with an Intel Core i5-1345U, 12 cores, and 16 GB RAM) took about 30 minutes. We recommend creating a dedicated environment to avoid dependency conflicts:

```bash
conda create -n gather-env -c bioconda gather
conda activate gather-env
```
#### Dependencies (only for Conda-based installation)
In addition to the Python modules installed through Conda, **GATHeR** requires the following external tools and R packages to support assembly, annotation, and phylogenetic analysis.

##### External Tools
- **NCBI BLAST+ v2.16.0** — sequence similarity searches
- **SPAdes v4.2.0** — de novo RNA assembly
- **IgPhyML** — phylogenetic inference of B-cell receptor lineages
  
##### R Packages
- **BiocManager v1.30.25** — package manager for Bioconductor
- **treeio v1.26.0**, **ggtree v3.14.0**, **dowser v1.10.0** — Bioconductor packages for lineage tree analysis
- **dplyr v1.1.4**, **ggrepel v0.9.5** — CRAN packages for data wrangling and visualization

To simplify installation, we provide an automated setup script:

📄 [`setup_gather_dependencies.sh`](https://github.com/Neuroimmunology-UiO/gather/blob/main/scripts/setup_dependencies.sh)

```bash
chmod +x setup_gather_dependencies.sh
./setup_gather_dependencies.sh
```

#### Prepare IgBLAST and Reference Databases (only for Conda-based installation)

To assign V(D)J genes and annotate junctional regions, GATHeR relies on [IgBLAST](https://www.ncbi.nlm.nih.gov/igblast/) and reference data from the IMGT database. Setting up the IgBLAST environment requires a few one-time configuration steps.

For convenience, we provide a setup script (based on the [Change-O IgBLAST guide](https://changeo.readthedocs.io/en/stable/examples/igblast.html)) that automates downloading, extracting, and configuring IgBLAST and its reference files.

📄 [`setup_igblast_env.sh`](https://github.com/Neuroimmunology-UiO/gather/blob/main/scripts/setup_igblast_env.sh)

```bash
chmod +x setup_igblast_env.sh
./setup_igblast_env.sh
```

> **Note:** You can adjust the IgBLAST version or installation directory inside the script if needed. After installation, restart your shell or run `source ~/.bashrc` to ensure the `IGDATA` environment variable is correctly set.

For troubleshooting and further details, see the [Change-O IgBLAST setup guide](https://changeo.readthedocs.io/en/stable/examples/igblast.html).

---

## Usage

GATHeR supports single-cell RNA sequencing data from technologies like Smart-seq2/3 and 10x Genomics Chromium.

### Paired-End Reads
```bash
sc_asm.py --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz --output_dir .
```

### Single-End Reads
```bash
sc_asm.py --seq_1 {cell_file_name}_R1.fastq.gz --output_dir .
```

### Output Files
- `{cell_file_name}_merged.unitigs.fa`: Unitigs from the cDBG (de Bruijn) graph
- `{cell_file_name}_merged.contigs.fa`: Assembled contigs (transcriptome, GATHeR algorithm)  
- `{cell_file_name}_merged.BCR.fa`: Annotated BCR sequences; headers indicate the algorithm (algo_1=SPAdes, algo_2=GATHeR)
  
### K-mer Size
- Default: `k = 25`
- Rule of thumb: `k < read length - 20`
- For short reads (e.g., 30 bp), the feasible k-mer size range is constrained; values between 21 and 25 are generally appropriate, with 27 being a possible upper limit depending on the context.

## Processing 10x Genomics Chromium 5' GEX Data for Assembly

To assemble data generated by 10x Genomics Chromium 5' gene-expression (GEX) single-cell RNA-seq, the raw sequencing reads must first be processed using Cell Ranger, which organizes the output files in the following structure:

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

## Clonality Analysis and Constant-Region Analysis

After assembling B-cell receptor (BCR) sequences for each cell, the reconstructed heavy and light chains—saved in separate FASTA files with distinct headers—can be collected across all cells for clonality and constant region analysis.
To prepare the input, copy all `{cell_file_name}_merged.BCR.fa` files from the individual cell directories (if they are not already in a shared location) into a single directory, e.g. `BCRs_DIR`.Then, you can run the `postproc.py` script to perform V(D)J gene assignment, CDR3 parsing, productivity assessment, and, optionally, clonal grouping and lineage reconstruction.

### Step 1: Analyze Light Chains
For light chains, use the following command:

```bash
postproc.py --bcrs_dir <PATH>/BCRs_DIR --chain light --output_dir .
```

This will produce:
- `light_chains_db-pass.tsv`: a **Change-O–like tab-delimited file** containing detailed annotations of the light chain sequences, including gene calls, junction regions, productivity, and alignment features.

### Step 2: Analyze Heavy Chains (Optional Clonal Inference and Visualization)

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
- Phylogenetic **lineage trees** of individual clones reconstructed from inferred germline sequences (`--lineage_tree`). This step can also be performed with [Dowser](https://github.com/immcantation/dowser?tab=readme-ov-file) using the above tables.

## License

Licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## 📖 Citation

If you use this software in your research, please cite:

> **GATHeR: Graph-based Accurate Tool for immunoglobulin HEavy- and light-chain Reconstruction**  
> Seyedraoufi, S., Gornitzka, M. B., Lossius, A.  
> *bioRxiv, 2025*  
> [https://doi.org/10.1101/2025.09.21.677530](https://doi.org/10.1101/2025.09.21.677530)
