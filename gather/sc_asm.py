#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 11:23:14 2025

@author: mojtaba
"""

from rich.console import Console
from gather import __version__, __author__, __description__, __url__
from utils import DBG, ProcessSeq, CommandExecutor, IO
import numpy as np
import networkx as nx
from tqdm import tqdm
import os
import argparse
import glob
import shutil

graph = DBG()
prcess = ProcessSeq()
com = CommandExecutor()
io = IO()


def display_custom_banner():
    console = Console()
    
    console.print('GATHER', style="bold cyan")
    console.print(f"Version: {__version__}", style="bold yellow")
    console.print(__description__, style="italic")
    console.print(f"Developed by {__author__}", style="green")
    console.print(f"GitHub: {__url__}", style="underline blue")
    console.print("\n")



def parse_arguments(package_path):
    parser = argparse.ArgumentParser(description='Process RNA sequences.')
    parser.add_argument('--seq_merged', required=True, help='Path to the merged seq in the paired-end sequences')
    parser.add_argument('--seq_1', required=True, help='Path to the paired-end sequence file 1')
    parser.add_argument('--seq_2', required=False, help='Path to the paired-end sequence file 2 (optional for single-end input)')
    parser.add_argument('--k', type=int, default=30, help='kmers size parameter, default is 25')
    parser.add_argument('--min_freq', type=int, default=5, help='Minimum frequency, default is 5')
    parser.add_argument('--output_dir', help='Output directory, default is the directory of the input_fastq file')
    parser.add_argument('--num_jobs', type=int, default=8, help='Number of jobs, default is -1')
    parser.add_argument('--merge_bcrs', action='store_true', help="If set, merge BCR sequences")
    parser.add_argument('--TCR', action='store_true', help="If set, assemble TCR sequences")
    parser.add_argument('--mouse', action='store_true', help="If set, use mouse database, otherwise use human")
    parser.add_argument('--bioalign', action='store_true', help="If set, use Biopython pairwise aligner")
    parser.add_argument('--use_spades', action='store_true', help="If set, use Spades assembler")
    
    parser.add_argument('--blast_dir',
                        default=os.path.dirname(shutil.which("blastn")) if shutil.which("blastn") else "~/Sys_admin/ncbi-blast-2.16.0+-x64-linux/ncbi-blast-2.16.0+/bin",
                        help='BLAST binary directory. Assumes blastn command is available in PATH.')
    parser.add_argument('--spades_path',
                        default=shutil.which("spades.py") if shutil.which("spades.py") else "~/Sys_admin/SPAdes-4.0.0-Linux/bin/spades.py",
                        help='SPAdes path. Assumes spades.py is available in PATH.')

    args = parser.parse_args()
    
    # Determine the correct database path
    if args.mouse:
        db_species = 'mouse'
    else:
        db_species = 'human'

    if args.TCR:
        db_type = 'TCR'
    else:
        db_type = 'BCR'

    args.db = os.path.join(package_path, 'database', db_species, db_type)

    # Set the default output_dir to be the directory of input_fastq if not specified
    if not args.output_dir:
        args.output_dir = os.path.dirname(os.path.abspath(args.seq_1))

    return args

def main():
    package_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    args = parse_arguments(package_path)
    display_custom_banner()

    # Using a single `scratch_dir` in this example as args.output_dir
    scratch_dir = args.output_dir
    concat_bcalm_file =  args.seq_merged
    input_rna_fastq_1 = args.seq_1
    input_rna_fastq_2 = args.seq_2 if args.seq_2 else None

    path_hc = os.path.join(args.db, 'IGHC.fasta')
    path_hv = os.path.join(args.db, 'IGHV.fasta')
    path_lc = os.path.join(args.db, 'IGLC.fasta')
    path_lv = os.path.join(args.db, 'IGLV.fasta')
    
    path_imgt_db =  os.path.join(args.db, 'IMGT.fasta') 
    
    k = args.k
    min_freq = args.min_freq
    
    os.chdir(scratch_dir)
    initial_items = set(os.listdir(scratch_dir))
    
    com.run_bcalm(concat_bcalm_file, k + 1, min_freq)    
    
    # input_unitigs_name = os.path.basename(concat_bcalm_file).replace('.gz', '.unitigs.fa')
    pattern = os.path.join(scratch_dir, '*.unitigs.fa')
    input_unitigs_name = glob.glob(pattern)[0]
    
    output_contigs_name = input_unitigs_name.replace('unitigs', 'contigs')
    
    
    
    input_spades_1 = os.path.join(scratch_dir, input_rna_fastq_1)
    input_spades_2 = os.path.join(scratch_dir, input_rna_fastq_2) if input_rna_fastq_2 else None
    input_path = os.path.join(scratch_dir, input_unitigs_name)
    output_path = os.path.join(scratch_dir, output_contigs_name)
    
    weights, overlap_indices, unitigs = io.parse_unitigs(input_path)
    
    unitigs_rc = [io.reverse_complement(i) for i in unitigs]
    
    begin, end = graph.convert_kmers_to_array(unitigs, k)
    begin_rc, end_rc = graph.convert_kmers_to_array(unitigs_rc, k)
    
    overlap_inds_list_0 = graph.find_matching_indices(begin, begin, begin_rc, weights,
                                                args.num_jobs)
    overlap_inds_list_1 = graph.find_matching_indices(end, end, end_rc, weights,
                                                args.num_jobs)
    
    rep_overlap_ind = overlap_inds_list_0 + overlap_inds_list_1
    
    overlap_indices_filtered = graph.remove_low_w_tuples(overlap_indices, rep_overlap_ind)
    overlap_indices_filtered = np.array(overlap_indices_filtered, dtype=object)
    
    G = nx.DiGraph()
    
    for i, unitig in enumerate(unitigs):
        G.add_node(i, sequence=unitig, weight=weights[i])
    
    for i, j, rc_1, rc_2 in overlap_indices_filtered:
        weight = (weights[i] + weights[j]) / 2.0
        tag = f"{rc_1}_{rc_2}"
        G.add_edge(i, j, weight=weight, tag=tag)
    
    wcc = list(nx.weakly_connected_components(G))
    wcc_s = [list(single_set)[0] for single_set in wcc if len(single_set) == 1]
    wcc_l = [l_set for l_set in wcc if len(l_set) > 1]
    
    optimum_paths, optimum_paths_W = graph.find_optimum_path(G, wcc_l, unitigs, weights)
    
    Contigs = [
        graph.merge_overlaps(path, overlap_indices_filtered, unitigs, unitigs_rc, k)
        for path in tqdm(
    optimum_paths,
    desc="Making contigs from optimum paths",
    bar_format="{l_bar}{bar:20}{r_bar}",
    colour="#808080",  
    )
    ]
    Contigs += [unitigs[i] for i in wcc_s if len(unitigs[i]) > 300]
    
    optimum_paths_W += [weights[i] for i in wcc_s if len(unitigs[i]) > 300]
    
    Contigs_headers = [f'{i}_{header}' for i, header in enumerate(optimum_paths_W)]
    
    io.save_to_fasta(headers=Contigs_headers, n_sequences=Contigs, file_path=output_path, include_header=True)
    
    header_hc, seqs_hc = io.read_fasta(path_hc, imgt_fasta=True)
    header_lc, seqs_lc = io.read_fasta(path_lc, imgt_fasta=True)
    header_hv, seqs_hv = io.read_fasta(path_hv, imgt_fasta=True)
    header_lv, seqs_lv = io.read_fasta(path_lv, imgt_fasta=True)
    
    seqs_hv = [seq.upper() for seq in seqs_hv]
    seqs_lv = [seq.upper() for seq in seqs_lv]
    seqs_hc = [seq.upper() for seq in seqs_hc]
    seqs_lc = [seq.upper() for seq in seqs_lc]
    
    header_c, seqs_contigs = io.read_fasta(output_path)
    header_c_float = [float(header[0]) for header in header_c]
    
    Contigs_headers = [f'{i}_{header}' for i, header in enumerate(header_c_float)]    
    io.save_to_fasta(headers=Contigs_headers, n_sequences=seqs_contigs, file_path=output_path, include_header=True)
    
    output_contigs_name = input_unitigs_name.replace('unitigs', 'BCR' if not args.TCR else 'TCR')
    output_contigs_name_spades = input_unitigs_name.replace('unitigs', 'BCR_contgus' if not args.TCR else 'TCR_contgus')
    output_bcr_path = os.path.join(scratch_dir, output_contigs_name)
    output_bcr_path_spades = os.path.join(scratch_dir, output_contigs_name_spades)
    
    output_spades_path = None
    seqs_contigs_spades = None
    if args.use_spades:
        spades_path = os.path.expanduser(args.spades_path)

        spades_success, spades_error = com.run_spades(spades_path, scratch_dir, input_spades_1, input_spades_2, is_rna=True)
        
        if not spades_success:
            print(f"Error running SPAdes: {spades_error}")
            print("We continue with the main algorithm!")
        else:
            print("SPAdes ran successfully.")
            output_spades_path = os.path.join(scratch_dir, 'transcripts.fasta')
            headers_spades, seqs_contigs_spades = io.read_fasta(output_spades_path)
    
    if args.bioalign:
        prcess.process_alignment(seqs_hc, header_hc, seqs_contigs, header_c, output_bcr_path, variable=False, num_jobs=args.num_jobs)
        prcess.process_alignment(seqs_hv, header_hv, seqs_contigs, header_c, output_bcr_path, variable=True, num_jobs=args.num_jobs)
        prcess.process_alignment(seqs_lc, header_lc, seqs_contigs, header_c, output_bcr_path, variable=False, num_jobs=args.num_jobs)
        prcess.process_alignment(seqs_lv, header_lv, seqs_contigs, header_c, output_bcr_path, variable=True, num_jobs=args.num_jobs)
       
    else:
        blast_dir = os.path.expanduser(args.blast_dir)
        blast_output_file = os.path.join(scratch_dir, "blast.out")
        prcess.process_alignment_blast(seqs_contigs, output_bcr_path, output_bcr_path_spades,
                                    blast_output_file, blast_dir, path_imgt_db, scratch_dir,
                                    output_path, output_spades_path, spades=args.use_spades,
                                    seqs_contigs_spades=seqs_contigs_spades, num_jobs=args.num_jobs)
    
    
        
    BCR_file_to_process = output_bcr_path_spades if args.use_spades else output_bcr_path   
        
    output_highest_expressed_BCR = input_unitigs_name.replace('unitigs', 'BCR_HE' if not args.TCR else 'TCR_HE')
    headers, sequences =  prcess.process_bcr_file(BCR_file_to_process)
    io.save_to_fasta(headers=headers, n_sequences=sequences, file_path=output_highest_expressed_BCR, include_header=True, append=False)
    
    current_items = set(os.listdir(scratch_dir))

    new_items = current_items - initial_items
    
    prcess.clean_scratch_dir(scratch_dir, new_items)
if __name__ == '__main__':
    main()
    