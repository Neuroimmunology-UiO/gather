#![Logo](gatherr.jpg "Our Logo")

# GATHeR: Graph-based Accurate Tool for immunoglobulin HEavy- and light-chain Reconstruction

This repository provides the environment setup for `gather`, including all necessary dependencies. Follow the instructions below to create and activate the environment using `conda`.

## Table of Contents

- [Prerequisites](#-prerequisites)
- [Installation](#-installation)
- [Usage](#-usage)
- [Contributing](#-contributing)
- [License](#-license)

## ✅ Prerequisites

Ensure you have Conda installed. If you don't have Conda, download and install it from the [official Anaconda website](https://www.anaconda.com/products/individual) or [Miniconda site](https://docs.conda.io/en/latest/miniconda.html).

## Installation

### Using the Environment File

To create and activate the `gather-env` environment using the provided `environment.yml` file:

1. **Clone the Repository**:

    ```bash
    git clone https://github.com/Neuroimmunology-UiO/gather.git
    cd gather
    ```

2. **Create the Environment**:

    ```bash
    conda env create -f environment.yml
    ```

3. **Activate the Environment**:

    ```bash
    conda activate gather-env
    ```

### Installing Dependencies Separately

1. **Create and activate a conda environment** (optional but recommended):

    ```bash
    conda create --name gather-env python=3.8
    conda activate gather-env
    ```

2. **Install Python dependencies**:

    ```bash
    pip install -r requirements.txt
    ```

3. **Install `bcalm` via Bioconda**:

    ```bash
    conda install -c bioconda bcalm
    ```

## Usage

**GATHeR** supports single-cell RNA sequencing data from technologies like Smart-seq2/3 and 10x Genomics Chromium.

### Merging Paired-End Reads

If your reads are paired-end:

```bash
cat {cell_file_name}_R1.fastq.gz {cell_file_name}_R2.fastq.gz > {cell_file_name}_merged.fastq.gz
```

Then run:

```bash
sc_asm.py --seq_merged {cell_file_name}_merged.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz
```

### Output Files

- `{cell_file_name}_merged.fastq.unitigs.fa`: Unitigs from the cDBG graph  
- `{cell_file_name}_merged.fastq.contigs.fa`: Assembled contigs (transcriptome, GATHeR)  
- `{cell_file_name}_merged.fastq.BCR.fa`: Annotated BCR sequences (GATHeR)

### Single-End Data

```bash
sc_asm.py --seq_merged {cell_file_name}_R1.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz
```

### SPAdes Mode for Low Coverage or Naive B Cells

Recommended for low read depth or naive B-cells:

```bash
sc_asm.py --seq_merged {cell_file_name}_merged.fastq.gz --seq_1 {cell_file_name}_R1.fastq.gz --seq_2 {cell_file_name}_R2.fastq.gz --min_freq 3 --use_spades
```

Additional outputs:

- `{cell_file_name}_merged.fastq.BCR_contigs.fa`: Annotated BCRs from SPAdes  
- `transcripts.fasta`: Assembled transcriptome (SPAdes)

### K-mer Size

- Default: `k = 25`
- ⚠️ Rule of thumb: `k < read length - 20`
- For read length = 30 bp ➝ choose k ≈ 21–25 (maybe 27)

## 🤝 Contributing

We welcome contributions!

1. Fork the repo
2. Create a feature branch
3. Commit changes clearly
4. Push to your fork
5. Open a pull request

## License

Licensed under the Apache License. See the [LICENSE](LICENSE) file for details.
