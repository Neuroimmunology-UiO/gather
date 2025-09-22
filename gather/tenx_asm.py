#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 10:15:29 2025

@author: mojtaba
"""

import os
import pandas as pd
import argparse
import sys
import gzip
from itertools import product, combinations
from tqdm import tqdm
from utils import CommandExecutor, IO
import glob
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed

com = CommandExecutor()
io = IO()



def parse_arguments():
    parser = argparse.ArgumentParser(description="Split paired FASTQ files by barcodes and run assembler.")

    parser.add_argument(
        '--barcode_file',
        type=str,
        required=True,
        help='Path to the TSV file containing barcodes. This typically contains cell barcodes extracted from 10x Genomics data.'
    )

    parser.add_argument(
        '--r1_fastq',
        type=str,
        required=True,
        help='Path to the R1 FASTQ.gz file. This file contains the read sequences along with the cell barcodes and UMIs (Unique Molecular Identifiers) in the case of 10x Genomics data.'
    )

    parser.add_argument(
        '--r2_fastq',
        type=str,
        required=True,
        help='Path to the R2 FASTQ.gz file. This file contains the complementary read sequences for the paired-end reads in 10x Genomics protocol.'
    )

    parser.add_argument(
        '--scratch_dir',
        type=str,
        default=os.getcwd(),
        help='Scratch directory for processing outputs. Defaults to the current working directory if not provided.'
    )

    parser.add_argument(
        '--output_dir_name',
        type=str,
        default='cells',
        help='Name of the output directory where cell directories will be created (default: cells).'
    )

    parser.add_argument(
        '--max_errors',
        type=int,
        default=1,
        help='Maximum allowed mismatches for barcode matching (default: 1).'
    )

    parser.add_argument(
        '--num_jobs',
        type=int,
        default=8,
        help='Number of jobs to use for the assembler (default: 8).'
    )
    
    parser.add_argument(
        '--clean_cells_data',
        action='store_true',
        help='If set, clean all files in each subdir except the final BCR contigs file.'
    )


    return parser.parse_args()

def read_barcodes(barcode_file):
    bc_df = pd.read_csv(barcode_file, sep="\t", header=None)
    barcodes_list = bc_df[0].tolist()
    return [i.split('-')[0] for i in barcodes_list]

def fastq_paired_reader(r1_handle, r2_handle):
    """
    A generator that yields paired-end FASTQ reads in raw line form:
    (r1_header, r1_seq, r1_plus, r1_qual,
     r2_header, r2_seq, r2_plus, r2_qual)
    """
    while True:
        r1_header = r1_handle.readline()
        if not r1_header:
            break  # End of file
        r1_seq = r1_handle.readline()
        r1_plus = r1_handle.readline()
        r1_qual = r1_handle.readline()

        r2_header = r2_handle.readline()
        r2_seq = r2_handle.readline()
        r2_plus = r2_handle.readline()
        r2_qual = r2_handle.readline()

        if not r2_qual:
            break  # Incomplete pair at end of file

        yield (
            r1_header.strip(),
            r1_seq.strip(),
            r1_plus.strip(),
            r1_qual.strip(),
            r2_header.strip(),
            r2_seq.strip(),
            r2_plus.strip(),
            r2_qual.strip()
        )


def generate_variants(barcode, max_errors=1):
    """
    Generate all variants of `barcode` with up to `max_errors` mismatches.
    Returns a set of variant strings.
    """
    bases = ['A', 'C', 'G', 'T']
    results = {barcode}
    length = len(barcode)
    
    for n_err in range(1, max_errors + 1):
        for positions in combinations(range(length), n_err):
            for subs in product(bases, repeat=n_err):
                variant_list = list(barcode)
                for i, pos in enumerate(positions):
                    variant_list[pos] = subs[i]
                results.add(''.join(variant_list))
    return results


def split_paired_by_barcode(
    r1_fastq,
    r2_fastq,
    output_dir,
    valid_barcodes,
    max_errors=0,
    barcode_length=16,
    chunk_size=1000
):
    """
    Process paired-end FASTQ files (gzipped) with chunked writes per barcode.
    This is significantly faster than opening and closing files on every read.
    
    :param r1_fastq: Path to R1.fastq.gz
    :param r2_fastq: Path to R2.fastq.gz
    :param output_dir: Directory for output
    :param valid_barcodes: An iterable of valid (canonical) barcodes
    :param max_errors: If > 0, allow up to N mismatches in barcodes
    :param barcode_length: Number of bases (from R1) to use as the barcode
    :param chunk_size: Flush chunk_size reads per barcode to disk at a time
    """
    os.makedirs(output_dir, exist_ok=True)

    # --- Build the lookup dictionary (barcode_map) ---
    #    For exact matches, you can directly set bc_lookup[bc] = bc
    #    For up to max_errors mismatches, precompute neighbors -> canonical BC
    bc_lookup = {}
    valid_barcodes = list(valid_barcodes)
    valid_barcodes_set = set(valid_barcodes)

    # For speed, map each valid barcode to itself first (exact match).
    for bc in valid_barcodes:
        bc_lookup[bc] = bc

    # If mismatch tolerance is used, build variants -> canonical
    if max_errors > 0:
        print("Precomputing barcode variants...", file=sys.stderr)
        for valid_bc in tqdm(valid_barcodes, desc="Preprocessing Barcodes", unit="barcode"):
            for variant in generate_variants(valid_bc, max_errors=max_errors):
                # If collisions matter, decide how to handle them
                # (Below we overwrite; you might want to skip or track collisions)
                bc_lookup[variant] = valid_bc

    # --- Prepare buffers for chunked writing ---
    # We'll store lines for R1 and R2 in memory, keyed by the "canonical" barcode.
    r1_buffers = {}
    r2_buffers = {}
    # We also track the total number of reads in each buffer so we know when to flush
    buffer_counts = {}

    # Instead of making directories and opening files for each barcode right away,
    # we'll do it lazily upon the first flush. We'll store them here:
    barcode_paths = {}

    def flush_barcode(bc):
        """
        Flush the buffered reads for barcode `bc` to disk, then clear its buffers.
        """
        if bc not in r1_buffers or not r1_buffers[bc]:
            return
        bc_dir = barcode_paths[bc]["dir"]
        out_r1_path = barcode_paths[bc]["r1_path"]
        out_r2_path = barcode_paths[bc]["r2_path"]

        # Create the directory if it doesn't exist
        os.makedirs(bc_dir, exist_ok=True)

        # Append the buffered lines using gzip in text append mode
        with gzip.open(out_r1_path, 'at') as f1:
            f1.writelines(r1_buffers[bc])
        with gzip.open(out_r2_path, 'at') as f2:
            f2.writelines(r2_buffers[bc])

        # Clear the buffers
        r1_buffers[bc] = []
        r2_buffers[bc] = []
        buffer_counts[bc] = 0

    def lazy_init_barcode_paths(bc):
        """
        Setup the directory/file paths for a given barcode if not already done.
        """
        if bc in barcode_paths:
            return
        bc_dir = os.path.join(output_dir, bc)
        barcode_paths[bc] = {
            "dir": bc_dir,
            "r1_path": os.path.join(bc_dir, f"{bc}_R1.fastq.gz"),
            "r2_path": os.path.join(bc_dir, f"{bc}_R2.fastq.gz")
        }


    # Initialize counters
    total_records = 0
    accepted_records = 0

    # --- Process the gzipped FASTQ files in parallel ---
    with gzip.open(r1_fastq, 'rt') as fh1, gzip.open(r2_fastq, 'rt') as fh2:
        reader = fastq_paired_reader(fh1, fh2)

        for (
            r1_header, r1_seq, r1_plus, r1_qual,
            r2_header, r2_seq, r2_plus, r2_qual
        ) in tqdm(reader, desc="Processing Records", unit="read"):

            total_records += 1

            # Extract the barcode from R1
            bc_candidate = r1_seq[:barcode_length]
            resolved_bc = None

            # Check exact or mismatch-based
            if bc_candidate in bc_lookup:
                resolved_bc = bc_lookup[bc_candidate]

            if resolved_bc:
                accepted_records += 1
                # Lazy init paths/buffers
                if resolved_bc not in r1_buffers:
                    r1_buffers[resolved_bc] = []
                    r2_buffers[resolved_bc] = []
                    buffer_counts[resolved_bc] = 0
                    lazy_init_barcode_paths(resolved_bc)

                # Store the lines in memory
                r1_buffers[resolved_bc].append(
                    f"{r1_header}\n{r1_seq}\n{r1_plus}\n{r1_qual}\n"
                )
                r2_buffers[resolved_bc].append(
                    f"{r2_header}\n{r2_seq}\n{r2_plus}\n{r2_qual}\n"
                )

                buffer_counts[resolved_bc] += 1

                # Flush if buffer is large
                if buffer_counts[resolved_bc] >= chunk_size:
                    flush_barcode(resolved_bc)

    # --- Final flush of all barcodes ---
    for bc in r1_buffers:
        flush_barcode(bc)

    # --- Summary ---
    print(
        f"Processed {total_records} read pairs. "
        f"Accepted {accepted_records} ({(accepted_records / total_records * 100):.2f}%).",
        file=sys.stderr
    )
    
def process_bcr_output(path_to_bcr):
    """Processes a BCR file to get sorted B-cell receptor sequences."""
    headers_BCR, seqs_BCR = IO.read_fasta(path_to_bcr, imgt_fasta=False)

    filtered_H = [
        (headers_BCR[i], seqs_BCR[i])
        for i in range(len(headers_BCR))
        if 'IGH' in headers_BCR[i][0].split(',')[0] and len(headers_BCR[i][0].split(',')) > 2
    ]

    weights_H = [
        float(headers_BCR[i][-1])
        for i in range(len(headers_BCR))
        if 'IGH' in headers_BCR[i][0].split(',')[0] and len(headers_BCR[i][0].split(',')) > 2
    ]
    
    if not filtered_H:
        filtered_H = [
            (headers_BCR[i], seqs_BCR[i])
            for i in range(len(headers_BCR))
            if 'IGH' in headers_BCR[i][0].split(',')[0]
        ]
        weights_H = [
            float(headers_BCR[i][-1])
            for i in range(len(headers_BCR))
            if 'IGH' in headers_BCR[i][0].split(',')[0]
        ]

    filtered_L = [
        (headers_BCR[i], seqs_BCR[i])
        for i in range(len(headers_BCR))
        if (
            ('IGK' in headers_BCR[i][0].split(',')[0] or 'IGL' in headers_BCR[i][0].split(',')[0])
            and len(headers_BCR[i][0].split(',')) > 2
        )
    ]

    weights_L = [
        float(headers_BCR[i][-1])
        for i in range(len(headers_BCR))
        if (
            ('IGK' in headers_BCR[i][0].split(',')[0] or 'IGL' in headers_BCR[i][0].split(',')[0])
            and len(headers_BCR[i][0].split(',')) > 2
        )
    ]

    if not filtered_L:
        filtered_L = [
            (headers_BCR[i], seqs_BCR[i])
            for i in range(len(headers_BCR))
            if 'IGK' in headers_BCR[i][0].split(',')[0] or 'IGL' in headers_BCR[i][0].split(',')[0]
        ]
        weights_L = [
            float(headers_BCR[i][-1])
            for i in range(len(headers_BCR))
            if 'IGK' in headers_BCR[i][0].split(',')[0] or 'IGL' in headers_BCR[i][0].split(',')[0]
        ]

    combined_H = sorted(zip(weights_H, filtered_H), key=lambda x: x[0], reverse=True)
    combined_L = sorted(zip(weights_L, filtered_L), key=lambda x: x[0], reverse=True)

    top_combined = combined_H[:] + combined_L[:]

    sorted_sequences = [seq for _, (_, seq) in top_combined]
    sorted_headers = [' '.join(header) for _, (header, _) in top_combined]

    return sorted_headers, sorted_sequences

def combine_bcrs(headers_1, headers_2, seqs_1, seqs_2):
    """
    Combine BCRs from the two assmbly algorithms to get the most of boths. 
    In the case of simialr BCRs, we append the one whith longest lenght.
    
    """
    
    headers_1_genes = [entry[0] for entry in headers_1]
    headers_2_genes = [entry[0] for entry in headers_2]
    
    unique_indices_list1 = [i for i, elem in enumerate(headers_1_genes) if elem not in headers_2_genes]
    unique_indices_list2 = [i for i, elem in enumerate(headers_2_genes) if elem not in headers_1_genes]
    
    headers_1_uq = [headers_1[i] for i in unique_indices_list1]
    headers_2_uq = [headers_2[i] for i in unique_indices_list2]
    
    seqs_1_uq = [seqs_1[i] for i in unique_indices_list1]
    seqs_2_uq = [seqs_2[i] for i in unique_indices_list2]
    
    shared_indices_list2 = [i for i, elem in enumerate(headers_2_genes) if elem in headers_1_genes]
    
    headers_2_sh = [headers_2[i] for i in shared_indices_list2]
    seqs_2_sh = [seqs_2[i] for i in shared_indices_list2]
    
    headers_f = headers_1_uq + headers_2_uq + headers_2_sh
    seqs_f = seqs_1_uq + seqs_2_uq + seqs_2_sh
    

    return headers_f, seqs_f




def clean_directory(subdir, path_to_keep):
    """Remove all files in subdir except for the specified path_to_keep."""
    for filename in os.listdir(subdir):
        file_path = os.path.join(subdir, filename)
        if file_path != path_to_keep and os.path.isfile(file_path):
            os.remove(file_path)


def run_single_cell_assembler_wrapper(
    subdir, 
    num_jobs, 
    scratch_dir, 
    clean_cells_data
):
    """
    A wrapper function that processes a single subdirectory (cell)
    with run_single_cell_assembler logic. This function is meant
    to be run in parallel.
    """
    # Path to the final aggregated BCR FASTA
    all_bcrs_fasta = os.path.join(scratch_dir, "ALL_BCRs.fasta")
    all_unitigs_fasta = os.path.join(scratch_dir, "ALL_unitigs.fasta")

    r1_fastq = glob.glob(os.path.join(subdir, '*_R1.fastq.gz'))
    r2_fastq = glob.glob(os.path.join(subdir, '*_R2.fastq.gz'))

    # Ensure both files are found
    if len(r1_fastq) == 1 and len(r2_fastq) == 1:
        try:
            # Determine the script's path (you can hard-code or find dynamically)
            script_name = "sc_asm.py"
            script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), script_name)

            # Check if the script exists
            if not os.path.exists(script_path):
                raise FileNotFoundError(f"Assembly script not found: {script_path}")

            success, error_message = com.run_single_cell_assembler(
                script_path=script_path,
                seq_1=r2_fastq[0],
                spades_path=shutil.which("spades.py"),
                blast_dir=os.path.dirname(shutil.which("blastn")),
                use_spades=True,
                min_freq=1,
                num_jobs=num_jobs
            )

            if success:
                print(f"Assembler ran successfully in {subdir}.")

                # Process BCR output
                basename = os.path.basename(subdir)
                path_to_bcr = os.path.join(subdir, f"{basename}_R2.BCR.fa")
                path_to_unitigs = os.path.join(subdir, f"{basename}_R2.unitigs.fa")

                # Check if neither of the key BCR files exists
                if not os.path.exists(path_to_bcr):
                    print(f"No BCR files found in {subdir}. Removing directory.")
                    shutil.rmtree(subdir)
                    return

                # Read existing FASTA files safely
                headers_bcr, sequences_bcr = [], []
             
                if os.path.exists(path_to_bcr):
                    headers_bcr, sequences_bcr = IO.read_fasta(path_to_bcr, imgt_fasta=False)

                headers_unitigs, sequences_unitigs = IO.read_fasta(path_to_unitigs, imgt_fasta=False)

                # Combine the available sequences
                headers = headers_bcr 
                sequences = sequences_bcr
                headers_refined = [' '.join(header) if isinstance(header, (list, tuple)) else str(header) for header in headers]

                # Add sample name as prefix
                headers = [f"{basename}, {i}" for i in headers_refined]

                # Save to main fasta
                io.save_to_fasta(
                    headers=headers, 
                    n_sequences=sequences, 
                    file_path=all_bcrs_fasta, 
                    append=True, 
                    include_header=True
                )
                
                
                headers_u = [f"{os.path.basename(subdir)}_{i}" for i,j in enumerate(headers_unitigs)]

                # Save to main fasta
                io.save_to_fasta(
                    headers=headers_u, 
                    n_sequences=sequences_unitigs, 
                    file_path=all_unitigs_fasta, 
                    append=True, 
                    include_header=True
                )

                # Clean directory if flag is set
                if clean_cells_data:
                    clean_directory(subdir, path_to_bcr)
   

            else:
                print(f"Error in {subdir}: {error_message}")
                print(f"Removing directory: {subdir}")
                shutil.rmtree(subdir)

        except Exception as e:
            print(f"An unexpected error occurred in {subdir}: {str(e)}")
            print(f"Removing directory: {subdir}")
            shutil.rmtree(subdir)


def process_batch(
    batch_id, 
    barcode_batch, 
    r1_fastq, 
    r2_fastq, 
    output_dir, 
    max_errors
):
    """
    Only does the split now. (Removed assembler calls here.)
    """
    print(f"Processing batch {batch_id}: {len(barcode_batch)} barcodes")

    # Split fastq by barcodes (this can be parallelized via multiple batches)
    split_paired_by_barcode(
        r1_fastq, 
        r2_fastq, 
        output_dir, 
        barcode_batch, 
        max_errors
    )
    # Return immediately after the split
    return


def main():
    args = parse_arguments()

    # Read the barcodes
    barcodes_list = read_barcodes(args.barcode_file)

    r1_fastq = args.r1_fastq
    r2_fastq = args.r2_fastq

    # Directory that will receive all split cells
    output_dir = os.path.join(args.scratch_dir, args.output_dir_name)
    os.makedirs(output_dir, exist_ok=True)

    num_jobs = args.num_jobs         # total CPU threads available
    max_workers = min(num_jobs, os.cpu_count())  # limit concurrency

    # Decide how to chunk the barcodes among the workers
    batch_size = max(1, len(barcodes_list) // max_workers)

    # 1) SPLIT PHASE: do splitting in batches (in parallel)
    # -----------------------------------------------------
    print("=== Splitting FASTQ by barcode (parallel) ===")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for i in range(0, len(barcodes_list), batch_size):
            batch_id = (i // batch_size) + 1
            barcode_batch = barcodes_list[i : i + batch_size]

            futures.append(
                executor.submit(
                    process_batch,
                    batch_id,
                    barcode_batch,
                    r1_fastq,
                    r2_fastq,
                    output_dir,
                    args.max_errors
                )
            )

        # Wait for all splits to finish
        for f in as_completed(futures):
            try:
                f.result()
            except Exception as e:
                print(f"Exception in one of the splitting batches: {e}")

    print("=== All barcode splitting is done. ===")

    # 2) ASSEMBLY PHASE: gather all cell dirs, then run the assembler
    # ---------------------------------------------------------------
    # By now, output_dir should have subdirectories for each cell/barcode
    cell_dirs = []
    for root, dirs, _ in os.walk(output_dir):
        for d in dirs:
            cell_dirs.append(os.path.join(root, d))    
    
    assembler_jobs = 4

    # How many parallel processes can we run if each uses assembler_jobs threads?
    # e.g. if max_workers=16 and assembler_jobs=4, concurrency = 4
    concurrency = max_workers // assembler_jobs
    concurrency = max(1, concurrency)  # avoid zero or negative

    print(f"=== Running assembler on {len(cell_dirs)} cells. ===")
    print(f"Max CPU threads: {max_workers}, each assembler job uses: {assembler_jobs}, concurrency: {concurrency}")

    # We use ProcessPoolExecutor to parallelize
    with ProcessPoolExecutor(max_workers=concurrency) as executor:
        futures = []
        for subdir in cell_dirs:
            futures.append(
                executor.submit(
                    run_single_cell_assembler_wrapper,
                    subdir,
                    assembler_jobs,        # each job uses 4 threads
                    args.scratch_dir,
                    args.clean_cells_data
                )
            )

        # Track progress with as_completed
        for f in tqdm(as_completed(futures), total=len(cell_dirs), desc="Assembling", unit="cell", colour='cyan'):
            try:
                f.result()
            except Exception as e:
                print(f"Exception in assembly: {e}")

    print("=== All assembly processes completed. ===")


if __name__ == '__main__':
    main()
