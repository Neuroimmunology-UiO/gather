# Single-Cell BCR Analysis (Smart-seq) — Code-Oriented Guide

Minimal, code-first walkthrough to assemble BCRs, post-process constant regions, run clonality, and generate lineage trees with **Singularity**.

## Example Data

Example input files can be downloaded from Figshare at the following link:


➡️ [Download example data from Figshare](https://doi.org/10.6084/m9.figshare.30152848)

---

## 0) Prerequisites

```bash
# Download container from:
# https://github.com/Neuroimmunology-UiO/gather/releases

# Verify the image is present (either local or absolute path)
ls -lh GATHER.sif  # or: ls -lh ~/Sys_admin/GATHER/GATHER.sif
```

(Optional) Set an absolute path once and reuse:
```bash
SIF="GATHER.sif"  # or: SIF="~/Sys_admin/GATHER/GATHER.sif"
```

---

## 1) Data Layout

```bash
# Work inside the dataset with FASTQs
cd Memory_B_cells/

# Expected naming:
#   <sample>_R1.fastq.gz
#   <sample>_R2.fastq.gz
ls *_R{1,2}.fastq.gz
```

---

## 2) Assemble BCRs (per cell)

```bash
# Outputs one FASTA per cell: <sample>__merged.BCR.fa
for r1 in *_R1.fastq.gz; do
  r2=${r1/_R1.fastq.gz/_R2.fastq.gz}
  sample=${r1%_R1.fastq.gz}

  echo "[assemble] ${sample}"

  singularity exec "${SIF:-GATHER.sif}" sc_asm.py \
    --seq_1 "$r1" \
    --seq_2 "$r2" \
    --output_dir . \
    --num_jobs 4
done
```

Expected files:
```text
<sample>_merged.BCR.fa
```

---

## 3) Post-Processing / Constant-Region Analysis

### Light chains
```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir . \
  --chain light \
  --output_dir "output"
```

### Heavy chains
```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir . \
  --chain heavy \
  --output_dir "output"
```

---

## 4) Clonality (Heavy Chain)

```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir . \
  --chain heavy \
  --output_dir "output" \
  --clonality
```
---

## 5) Lineage Trees (per clone)

```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir . \
  --chain heavy \
  --output_dir "output" \
  --clonality \
  --lineage_tree
```
---

## 6) Tips / Troubleshooting

```bash
# Ensure file pairing exists
for r1 in *_R1.fastq.gz; do
  test -f "${r1/_R1.fastq.gz/_R2.fastq.gz}" || echo "Missing R2 for $r1"
done

# Use absolute path to the image if needed
SIF="<PATH>/GATHER.sif"

# Adjust CPU usage during assembly (default shown: 4)
#   --num_jobs <N>
