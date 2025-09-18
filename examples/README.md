# Single-Cell BCR Analysis (Smart-seq) — Code-Oriented Guide

Minimal, code-first walkthrough to assemble BCRs, post-process constant regions, run clonality, and generate lineage trees with **Singularity**.

---

## 0) Prerequisites

```bash
# Download container from:
# https://github.com/Neuroimmunology-UiO/gather/releases

# Verify the image is present (either local or absolute path)
ls -lh GATHER.sif  # or: ls -lh ~/Sys_admin/GATHER/GATHER.sif
```

(Optional) set an absolute path once and reuse:
```bash
SIF="GATHER.sif"  # or: SIF="~/Sys_admin/GATHER/GATHER.sif"
```

---

## 1) Data Layout

```bash
# Work inside the dataset with FASTQs
cd memory_B_cells/

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
<sample>__merged.BCR.fa
```

---

## 3) Collect Assembled BCRs

```bash
mkdir -p BCRs
cp *__merged.BCR.fa BCRs/
bcr_dir="BCRs"

# Optional: preview a few headers
grep -m 5 -E "^>" BCRs/* | head -n 10
```

---

## 4) Post-Processing / Constant-Region Analysis

### Light chains
```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir "$bcr_dir" \
  --chain light \
  --output_dir "output"
```

### Heavy chains
```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir "$bcr_dir" \
  --chain heavy \
  --output_dir "output"
```

---

## 5) Clonality (Heavy Chain)

```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir "$bcr_dir" \
  --chain heavy \
  --output_dir "output" \
  --clonality
```
---

## 6) Lineage Trees (per clone)

```bash
singularity exec "${SIF:-GATHER.sif}" postproc.py \
  --bcrs_dir "$bcr_dir" \
  --chain heavy \
  --output_dir "output" \
  --clonality \
  --lineage_tree
```
---

## 7) Tips / Troubleshooting

```bash
# Ensure file pairing exists
for r1 in *_R1.fastq.gz; do
  test -f "${r1/_R1.fastq.gz/_R2.fastq.gz}" || echo "Missing R2 for $r1"
done

# Use absolute path to the image if needed
SIF="~/Sys_admin/GATHER/GATHER.sif"

# Adjust CPU usage during assembly (default shown: 4)
#   --num_jobs <N>

# Run from anywhere by using absolute paths
bcr_dir="$(pwd)/BCRs"
outdir="$(pwd)/output_heavy"
singularity exec "${SIF}" postproc.py --bcrs_dir "$bcr_dir" --chain heavy --output_dir "$outdir"
