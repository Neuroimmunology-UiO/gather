# Single-Cell BCR Analysis Example (Smart-seq Data)

This section provides an example workflow using **Smart-seq single-cell RNA-seq data** from:  
- 6 human plasmablasts  
- 10 memory B cells  

## Features

By following the provided steps, you can:  
- Assemble heavy and light chains  
- Perform constant-region novel allele analysis  
- Conduct clonality analysis  
- Construct lineage trees  

---

## Example: Memory B Cells

In the folder `memory_B_cells/` you will find FASTQ files for **10 single cells**.  
To process them, please use our Singularity image.

### 1. Download the Singularity image
Download the required container (`GATHER.sif`) from [Releases section](https://github.com/Neuroimmunology-UiO/gather/releases).

### 2. Run the assembly loop
From inside the `memory_B_cells/` folder, run the following command:

```bash
for r1 in *_R1.fastq.gz
do
    # Find matching R2 file
    r2=${r1/_R1.fastq.gz/_R2.fastq.gz}

    # Extract sample name (remove suffix)
    sample=${r1%_R1.fastq.gz}

    echo "Running assembly for: $sample"

    singularity exec GATHER.sif sc_asm.py \
        --seq_1 "$r1" \
        --seq_2 "$r2" \
        --output_dir . \
        --num_jobs 4
done
```
### 2. Collect BCRs and perform 


