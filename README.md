#![Logo](gatherr.jpg "Our Logo")

# GATHeR: Graph-based Accurate Tool for immunoglobulin HEavy- and light-chain Reconstruction

This repository provides the environment setup for `gather`, including all necessary dependencies. Follow the instructions below to create and activate the environment using `conda`.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites

Ensure you have Conda installed. If you don't have Conda, download and install it from the [official Anaconda website](https://www.anaconda.com/products/individual) or [Miniconda site](https://docs.conda.io/en/latest/miniconda.html).

## Installation

### Using the Environment File

To create and activate the `gather-env` environment using the provided `environment.yml` file, follow these steps:

1. **Clone the Repository**:

    ```bash
    git clone https://github.com/Neuroimmunology-UiO/gather.git
    cd gather
    ```

2. **Create the Environment**:

    Use the `environment.yml` file to create the environment. This will install all necessary dependencies.

    ```bash
    conda env create -f environment.yml
    ```

3. **Activate the Environment**:

    Once the environment is created, activate it using the following command:

    ```bash
    conda activate gather-env
    ```

### Installing Dependencies Separately

If you prefer to install dependencies manually or using `requirements.txt`, you can follow these instructions:

1. **Create and activate a conda environment** (optional but recommended):

    ```bash
    conda create --name gather-env python=3.8
    conda activate gather-env
    ```

2. **Install the dependencies using pip**:

    ```bash
    pip install -r requirements.txt
    ```

3. **Install `bcalm` (if not available on PyPI)**:

    ```bash
    conda install -c bioconda bcalm
    ```

## Usage

**GATHeR** is compatible with popular single-cell RNA sequencing (scRNA-seq) technologies, including plate-based full-length transcript methods (Smart-seq2 / Smart-seq3), and 10x Genomics Chromium.

If you have scRNA-seq data (separated for each cell) in `fastq` or `fastq.gz` format, then for paired-end data, you first need to merge the reads. In Linux, this can be done with:

```bash
cat {cell_file_name}_R1.fastq.gz {cell_file_name}_R2.fastq.gz > {cell_file_name}_merged.fastq.gz
```

Then run:

```bash
sc_asm.py --seq_merged {cell_file_name}_merged.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz
```

### Output Files

After a successful run, the following files will be generated:

- `{cell_file_name}_merged.fastq.unitigs.fa` – Unitigs from the compact de Bruijn graph  
- `{cell_file_name}_merged.fastq.contigs.fa` – Assembled contigs (transcriptome) using the GATHeR algorithm  
- `{cell_file_name}_merged.fastq.BCR.fa` – Annotated BCR sequences using the GATHeR algorithm  

### Single-End Data

For single-end reads, you can run:

```bash
sc_asm.py --seq_merged {cell_file_name}_R1.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz
```

### SPAdes Mode for Low Read Depth / Naive B Cells

If dealing with low read depth or naive B-cells, we recommend using the `--use_spades` option and reducing the `--min_freq` parameter (recommended range: 1–5):

```bash
sc_asm.py --seq_merged {cell_file_name}_merged.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz --min_freq 3 --use_spades
```

This produces additional outputs:

- `{cell_file_name}_merged.fastq.BCR_contigs.fa` – Annotated BCRs from SPAdes  
- `transcripts.fasta` – Assembled contigs (transcriptome) from SPAdes  

### Notes on K-mer Size

- Default k-mer size: `25`
- Recommended rule: `k < read length - 20`
- For short reads (e.g., 30 bp), valid k-mer sizes are around 21–25, maybe up to 27.

## Contributing

If you wish to contribute to this project, please follow these steps:

1. **Fork the repository**.
2. **Create a new branch** for your feature or bug fix.
3. **Make your changes** and commit them with clear, descriptive messages.
4. **Push your changes** to your forked repository.
5. **Create a pull request** to have your changes reviewed and merged into the main repository.

## License

This project is licensed under the Apache License. See the [LICENSE](LICENSE) file for more details.
