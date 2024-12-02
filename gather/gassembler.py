#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:23:05 2024

@author: NI group
"""

import numpy as np
import re
import networkx as nx
from joblib import Parallel, delayed
from tqdm import tqdm
from itertools import product
from Bio import Align, SeqIO
import pandas as pd
import os
from Bio.Seq import Seq
import subprocess
from difflib import SequenceMatcher
import argparse


def reverse_complement(dna):
    """
    Compute the reverse complement of a DNA sequence.

    Args:
        dna (str): The DNA sequence.

    Returns:
        str: The reverse complement of the DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    try:
        reverse_comp = ''.join(complement[base] for base in reversed(dna))
    except KeyError as e:
        raise ValueError(f"Invalid DNA base found: {e}")
        
    return reverse_comp

def parse_unitigs(fasta_file):
    """
    Parse unitigs from a FASTA file.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        tuple: A tuple containing:
            - A list of km values (list of floats).
            - A list of tuples containing read index, integer between L patterns, and signs.
            - A list of sequences (list of strings).
    """
    
    km_values = []
    tuple_list = []
    sequences = []
    current_seq = ""
    
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):  # It's a header line
                if current_seq:
                    sequences.append(current_seq)
                    current_seq = ""
                
                header = line.strip()
                
                # Extract the km value
                km_match = re.search(r'km:f:(\d+\.\d+)', header)
                if km_match:
                    km_values.append(float(km_match.group(1)))
                
                # Extract the read index
                read_index_match = re.match(r'>(\d+)', header)
                if read_index_match:
                    read_index = int(read_index_match.group(1))
                
                # Find all L patterns and process them
                L_patterns = re.findall(r'L:([+|-]):(\d+):([+|-])', header)
                for L in L_patterns:
                    middle_int = int(L[1])
                    first_sign = 'f' if L[0] == '+' else 'r'
                    second_sign = 'f' if L[2] == '+' else 'r'
                    tuple_list.append((read_index, middle_int, first_sign, second_sign))

            else:
                current_seq += line.strip()
        
        # Append the last sequence
        if current_seq:
            sequences.append(current_seq)

    return km_values, tuple_list, sequences

def save_to_fasta(headers=None, n_sequences=None, file_path=None, include_header=False, append=False):
    """
    Save a list of DNA sequences to a file in FASTA format.

    :param headers: List of headers for each sequence if include_header is True or None
    :param n_sequences: List of DNA sequences as strings
    :param file_path: Path to save the output FASTA file
    :param include_header: Boolean indicating whether to use provided headers
    :param append: Boolean indicating whether to append to the file (True) or overwrite it (False)
    """
    if n_sequences is None:
        raise ValueError("n_sequences should not be None")

    if file_path is None:
        raise ValueError("file_path should not be None")

    mode = 'a' if append else 'w'

    with open(file_path, mode) as fasta_file:
        for i, sequence in enumerate(n_sequences):
            if include_header and headers:
                fasta_file.write(f">{headers[i]}\n")
            else:
                seq_id = f"seq{i+1}"  # Generate a unique identifier
                fasta_file.write(">{seq_id}\n")
            
            fasta_file.write(sequence + "\n")


def read_fasta(file_path, imgt_fasta=False):
    """
    Read a FASTA file and extract headers and sequences.

    :param file_path: The path to the FASTA file to read.
    :param imgt_fasta: Boolean indicating whether the FASTA file follows IMGT format.
    :return: Tuple containing:
             - A list of headers, each header is a list of parts (depending on the format).
             - A list of sequences as strings.
    """
    headers = []
    sequences = []
    with open(file_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            if imgt_fasta:
                # Complex header handling for IMGT FASTA files
                header_parts = record.description.split('|')
                header_parts = [part.strip() if part.strip() else 'None' for part in header_parts]
            else:
                # Simple header handling
                header_parts = record.description.split()
            headers.append(header_parts)
            sequences.append(str(record.seq))
    return headers, sequences


def run_bcalm(input_file, kmer_size, abundance_min):
    """
    Run the bcalm command with specified parameters.

    :param input_file: Path to the input file.
    :param kmer_size: Kmer size to use.
    :param abundance_min: Minimum abundance threshold.
    :raises RuntimeError: If bcalm command fails.
    """
    command = [
        "bcalm",  # Command to run
        "-in", input_file,  # Input file
        "-kmer-size", str(kmer_size),  # Kmer size
        "-abundance-min", str(abundance_min), # Minimum abundance
    ]

    # Run the command and suppress output
    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    # Check if the command was successful
    if result.returncode == 0:
        print("BCALM ran successfully!")
    else:
        print("Error running bcalm:")
        raise RuntimeError("BCALM command failed with return code {}".format(result.returncode))



# Define a mapping from DNA characters to float values
DNA_TO_FLOAT = {'A': 1.0, 'T': 2.0, 'C': 3.0, 'G': 4.0}

def convert_kmers_to_array(list_of_strings, k):
    """
    Convert k-mers in each string of a list to 2D numpy arrays of float values.

    :param list_of_strings: List of DNA sequences as strings.
    :param k: Length of k-mers to extract from the beginning and end of each sequence.
    :return: Tuple of two 2D numpy arrays: beginnings and ends.
    """
    
    num_strings = len(list_of_strings)
    
    # Initialize two 2D numpy arrays for beginnings and ends
    beginnings = np.zeros((num_strings, k), dtype=float)
    ends = np.zeros((num_strings, k), dtype=float)
    
    for i, dna_str in enumerate(list_of_strings):
        # Extract and convert the beginning and end substrings
        beginning_str = dna_str[:k]
        end_str = dna_str[-k:]
        
        # Convert to float using the mapping
        beginnings[i] = [DNA_TO_FLOAT[ch] for ch in beginning_str]
        ends[i] = [DNA_TO_FLOAT[ch] for ch in end_str]
    
    return beginnings, ends


def find_matching_indices(array1, array2, array2_rc, weights, num_jobs):
    """
    Finds the indices of rows that are equal in two given 2D NumPy arrays
    The function also considers the reverse complement of array2.
    
    Parameters:
    array1 (np.ndarray): First 2D NumPy array.
    array2 (np.ndarray): Second 2D NumPy array.
    weights (np.ndarray): 1D NumPy array of weights corresponding to rows of arrays.
    
    Returns:
    list of tuples: list of indices that have lower weight comparing to their equal pairs.
    """
    
    def compare_arrays(array1, array2, num_jobs):
        rows, cols = array1.shape[0], array2.shape[0]

        # Split work to chunks to be processed by workers
        def compare_chunk(start, end):
            result_indices = []
            for i in range(start, end):
                comparison_result = np.all(array1[i, np.newaxis] == array2, axis=1)
                matched_indices = np.argwhere(comparison_result).flatten()
                for m_idx in matched_indices:
                    result_indices.append((i, m_idx))
            return result_indices

        # Split indices according to the number of jobs
        chunk_size = (rows + num_jobs - 1) // num_jobs
        chunks = [(i, min(i + chunk_size, rows)) for i in range(0, rows, chunk_size)]

        # Run the parallel computation in chunks
        parallel_results = Parallel(n_jobs=num_jobs)(
            delayed(compare_chunk)(start, end) for start, end in tqdm(chunks, desc='Processing unitigs', colour="#800080")
        )
        
        # Flatten the list of results
        all_indices = [indice for sublist in parallel_results for indice in sublist]

        return all_indices

    equal_indices = compare_arrays(array1, array2, num_jobs)
    equal_indices_rc = compare_arrays(array1, array2_rc, num_jobs)
    

    equal_indices_low_w_norc =[i if weights[j] > weights[i] else j for i, j in equal_indices if i != j]
    equal_indices_low_w_rc =[i if weights[j] > weights[i] else j for i, j in equal_indices_rc if i != j]
    
    equal_indices_low_w = equal_indices_low_w_norc + equal_indices_low_w_rc
    
    return equal_indices_low_w

def find_overlap_indices(array1, array2, array1_rc, array2_rc, num_jobs):
    """
    Finds the indices of rows that are equal in two given 2D NumPy arrays
    and calculates the average weight for each matching row pair. 
    The function also considers the reverse complement of array2.
    
    Parameters:
    array1 (np.ndarray): First 2D NumPy array.
    array2 (np.ndarray): Second 2D NumPy array.
    weights (np.ndarray): 1D NumPy array of weights corresponding to rows of arrays.
    
    Returns:
    list of tuples: Each tuple contains the indices of matching rows (i, j),
                    the average weight (w_avg), and a flag indicating if the row
                    in array2 is direct or reverse complement.
    """
   
    
    def compare_arrays(array1, array2, num_jobs):
        rows, cols = array1.shape[0], array2.shape[0]

        # Split work to chunks to be processed by workers
        def compare_chunk(start, end):
            result_indices = []
            for i in range(start, end):
                comparison_result = np.all(array1[i, np.newaxis] == array2, axis=1)
                matched_indices = np.argwhere(comparison_result).flatten()
                for m_idx in matched_indices:
                    result_indices.append((i, m_idx))
            return result_indices

        # Split indices according to the number of jobs
        chunk_size = (rows + num_jobs - 1) // num_jobs
        chunks = [(i, min(i + chunk_size, rows)) for i in range(0, rows, chunk_size)]

        # Run the parallel computation in chunks
        parallel_results = Parallel(n_jobs=num_jobs)(
            delayed(compare_chunk)(start, end) for start, end in tqdm(chunks, desc='Processing unitigs', colour="#800080", ncols=100)
        )
        
        all_indices = [indice for sublist in parallel_results for indice in sublist]

        return all_indices
    
    equal_indices = compare_arrays(array1, array2, num_jobs)
    equal_indices_rc = compare_arrays(array1, array2_rc, num_jobs)
    equal_indices_rc_1 = compare_arrays(array1_rc, array2, num_jobs)
    equal_indices_rc_2 = compare_arrays(array1_rc, array2_rc, num_jobs)

    # Create the result list with (i, j, direction_i, direction_j)
    overlap_1 =[ (i, j, 'f', 'f') for i, j in equal_indices if i != j]
    overlap_2 =[ (i, j, 'f', 'r') for i, j in equal_indices_rc if i != j]
    overlap_3 =[ (i, j, 'r', 'f') for i, j in equal_indices_rc_1 if i != j]
    overlap_4 =[ (i, j, 'r', 'r') for i, j in equal_indices_rc_2 if i != j]
    
    overlap_indices_total = overlap_1 + overlap_2 + overlap_3 + overlap_4
    
    return overlap_indices_total


def remove_low_w_tuples(tuples_list, indices_array):
    """
    Remove tuples from a list where the first or second element is present in a provided array.

    :param tuples_list: List of tuples, each containing (int, int, str, str).
    :param indices_array: A numpy array of indices.
    :return: List of filtered tuples.
    """
    # Convert the tuples list to a NumPy structured array to handle mixed data types
    dtype = [('first', int), ('second', int), ('third', 'U1'), ('fourth', 'U1')]
    structured_array = np.array(tuples_list, dtype=dtype)

    # Create a mask, using np.isin to check for the presence of the values in indices_array
    mask_first = np.isin(structured_array['first'], indices_array)
    mask_second = np.isin(structured_array['second'], indices_array)

    # Combine the masks and invert to get the rows we want to keep
    mask = ~mask_first & ~mask_second

    # Apply the mask to filter out the unwanted rows
    filtered_structured_array = structured_array[mask]

    # Convert filtered structured array back to a list of tuples
    filtered_tuples_list = [
        (row['first'], row['second'], row['third'], row['fourth']) 
        for row in filtered_structured_array
    ]
    
    return filtered_tuples_list


def dfs_all_paths(graph, start_node):
    """
    Perform a depth-first search to find all paths in the graph from the start node.

    :param graph: A networkx graph.
    :param start_node: The starting node for DFS.
    :return: A list of all paths from the start_node.
    """
    stack = [(start_node, [], set(), None)]
    paths = []

    while stack:
        current_node, path, visited, prev_edge_data = stack.pop()

        path.append(current_node)
        visited.add(current_node)

        is_leaf = True

        for neighbor in graph.neighbors(current_node):
            edge_data = graph.get_edge_data(current_node, neighbor)

            if neighbor not in visited:
                if prev_edge_data is None or edge_data['tag'][0] == prev_edge_data['tag'][-1]:
                    is_leaf = False
                    stack.append((neighbor, path.copy(), visited.copy(), edge_data))

        if is_leaf:
            paths.append(path)

    return paths

def find_optimum_path(G, wcc, unitigs, weights):
    """
    Find the optimum path in a subgraph of G defined by the weakly connected components (WCC).

    :param G: The input graph (networkx graph).
    :param wcc: A lsit conatining WCC sets.
    :param unitigs: an array of unitigs.
    :param weights: an array of weigths corresponding to unitigs.
    :return: List of nodes representing the optimum path and their weights. 
    """
    longest_paths = []
    long_w = []

    for cluster in tqdm(wcc, desc="Finding the optimum path in WCCs", colour="#800080", ncols=100):
        longest_path_size = len(cluster)
        max_path_avg_weight = 0
        longest_cluster_path = None
        
        for start_node in cluster:
            paths = dfs_all_paths(G, start_node)
            for path in paths:
                path_avg_weight = np.sum(np.fromiter((weights[node] for node in path), dtype=float)) / longest_path_size
                if path_avg_weight > max_path_avg_weight:
                    longest_cluster_path = path
                    max_path_avg_weight = path_avg_weight
        
        if longest_cluster_path is not None:
            longest_paths.append(longest_cluster_path)
            long_w.append(max_path_avg_weight)
    
    return longest_paths, long_w

def merging_unitigs(overlap_info, unitigs_max, unitigs, unitigs_rc, k, first_index=False):
    """
    Merge two unitigs based on overlap information and their orientations.

    :param overlap_info: Tuple containing (ind_1, ind_2, rc_1, rc_2) where
                         - ind_1: Index of the first unitig.
                         - ind_2: Index of the second unitig.
                         - rc_1: Orientation of the first unitig ('f' for forward, 'r' for reverse complement).
                         - rc_2: Orientation of the second unitig ('f' for forward, 'r' for reverse complement).
    :param unitigs_max: The current maximum unitig sequence being assembled.
    :param unitigs: List of unitig sequences.
    :param unitigs_rc: List of reverse complement unitig sequences.
    :param k: The length of the k-mer overlap.
    :param first_index: Boolean indicating if this is the first index to be merged.
    :return: The updated maximum unitig sequence after merging.
    """
    ind_1, ind_2, rc_1, rc_2 = overlap_info
    if first_index:
        if rc_1 == 'f' and rc_2 == 'f':
            unitigs_max = unitigs[ind_1]
            adding_unitig = unitigs[ind_2]
            unitigs_max = unitigs_max + adding_unitig[k:]
        elif rc_1 == 'f' and rc_2 == 'r':
            unitigs_max = unitigs[ind_1]
            adding_unitig = unitigs_rc[ind_2]
            unitigs_max = unitigs_max + adding_unitig[k:]        
        elif rc_1 == 'r' and rc_2 == 'f':
            unitigs_max = unitigs_rc[ind_1]
            adding_unitig = unitigs[ind_2]
            unitigs_max = unitigs_max + adding_unitig[k:]  
        else:
            unitigs_max = unitigs_rc[ind_1]
            adding_unitig = unitigs_rc[ind_2]
            unitigs_max = unitigs_max + adding_unitig[k:]
    else:
        if rc_2 == 'f':
            adding_unitig = unitigs[ind_2]
            unitigs_max = unitigs_max + adding_unitig[k:]
        else:
            adding_unitig = unitigs_rc[ind_2]
            unitigs_max = unitigs_max + adding_unitig[k:]
    
    return unitigs_max


def merge_overlaps(indices, overlap_inds_list, unitigs, unitigs_rc, k):
    """
    Merge multiple unitigs based on overlap information and given indices.

    :param indices: List of indices to merge.
    :param overlap_inds_list: List of tuples containing overlap information (ind_1, ind_2, rc_1, rc_2).
    :param unitigs: List of unitig sequences.
    :param unitigs_rc: List of reverse complement unitig sequences.
    :param k: The length of the k-mer overlap.
    :return: The merged unitig sequence.
    """
    indices_array = np.array(indices, dtype=np.int64)
    
    mask = np.isin(overlap_inds_list[:, 0], indices_array) | np.isin(overlap_inds_list[:, 1], indices_array)
    
    filtered_oil = overlap_inds_list[mask]
    unitigs_max = None
    for i in range(1, len(indices)):
        current_index = indices[i - 1]
        next_index = indices[i]
        # overlap_info = next((item for item in overlap_inds_list if item[0] == current_index and item[1] == next_index), None)
        mask = (filtered_oil[:, 0] == np.int64(current_index)) & (filtered_oil[:, 1] == np.int64(next_index))
        overlap_info = filtered_oil[mask][0]
        # Check and initialize unitigs_max when current_index is 0
        if (i-1) == 0 and unitigs_max is None:
            
            unitigs_max = merging_unitigs(overlap_info, unitigs_max, unitigs, unitigs_rc, k, first_index=True)
        else:
            unitigs_max = merging_unitigs(overlap_info, unitigs_max, unitigs, unitigs_rc, k, first_index=False)

        
    return unitigs_max

class AlignmentError(Exception):
    pass

def compute_aligned_fraction(identities, gaps, mismatches):
    """
    Compute the fraction of aligned characters.

    :param identities: Number of matching bases.
    :param gaps: Number of gaps in the alignment.
    :param mismatches: Number of mismatching bases.
    :return: Aligned fraction.
    """
    return identities / (identities + gaps + mismatches)

def local_alignment_best(seq1, seq2, match_score=3, mismatch_penalty=-3, gap_open_penalty=-10, gap_extend_penalty=-1):
    """
    Perform local alignment between two sequences and return the best alignment.

    :param seq1: First sequence.
    :param seq2: Second sequence.
    :param match_score: Score for a match.
    :param mismatch_penalty: Penalty for a mismatch.
    :param gap_open_penalty: Penalty for opening a gap.
    :param gap_extend_penalty: Penalty for extending a gap.
    """
    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_penalty
    aligner.open_gap_score = gap_open_penalty
    aligner.extend_gap_score = gap_extend_penalty

    # Perform the alignment
    alignments = aligner.align(seq1, seq2)

    # Check if alignments are found
    if not alignments:
        raise AlignmentError("No alignments found between sequences.")

    # Extract the alignment with the highest score
    best_alignment = alignments[0]
    aligned_counts = list(best_alignment.counts())
    gaps, identities, mismatches = aligned_counts
    
    return best_alignment.score, gaps, identities, mismatches, best_alignment


def process_combination(i, j, seq1, seq2):
    """
    Process a combination of sequence indices and perform alignments.
    
    :param i: Index of the first sequence.
    :param j: Index of the second sequence.
    :param seq1: First sequence.
    :param seq2: Second sequence.
    """
    try:
        score_forward, gaps_forward, identities_forward, mismatches_forward, _ = local_alignment_best(seq1, seq2)
        score_reverse, gaps_reverse, identities_reverse, mismatches_reverse, _ = local_alignment_best(seq1, reverse_complement(seq2))

        if score_forward > score_reverse:
            A_forward = compute_aligned_fraction(identities_forward, gaps_forward, mismatches_forward)
            return score_forward, A_forward, (i, j, 'forward')
        else:
            A_reverse = compute_aligned_fraction(identities_reverse, gaps_reverse, mismatches_reverse)
            return score_reverse, A_reverse, (i, j, 'reverse')
    except AlignmentError as e:
        # Log or handle the situation where no alignment was found
        print(f"AlignmentError for pair (seq1 index {i}, seq2 index {j}): {e}")
        # Return a score and A indicating no alignment could be made
        return float('-inf'), 0.0, (i, j, None)

def find_best_alignments(seqs1, seqs2, num_jobs):
    """
    Find the best alignments between two sets of sequences in parallel.

    :param seqs1: List of sequences from the first set.
    :param seqs2: List of sequences from the second set.
    :param num_jobs: Number of parallel jobs.
    :return: Sorted list of results based on alignment scores.
    """
    
    filtered_seqs2 = [s for s in seqs2 if len(s) > 0]
    
    # Use itertools.product to create the Cartesian product of the indices
    all_combinations = list(product(range(len(seqs1)), range(len(filtered_seqs2))))
    total_combinations = len(seqs1) * len(filtered_seqs2)
    
    results = Parallel(n_jobs=num_jobs)(
        delayed(process_combination)(i, j, seqs1[i], filtered_seqs2[j]) 
        for i, j in tqdm(all_combinations, total=total_combinations, desc='Processing pairwise alignments', colour="#800080", ncols=100)
    )
    
    # Sort the results based on A in descending order
    sorted_results = sorted(results, key=lambda x: x[0], reverse=True)

    return sorted_results

# def find_most_similar(data, search_list):
#     """
#     Finds the most similar string in search_list that matches a substring in the data.

#     :param data: The data string containing substrings.
#     :param search_list: List of strings to match against.
#     :return: The closest match found or None if no match is found.
#     """
#     substrings = data.split()
#     close_matches = []

#     # Loop through the substrings and compare with items in the search list
#     for substring in substrings:
#         if substring.startswith("IG"):
#             # Find the closest match to the current substring in the search list
#             match = difflib.get_close_matches(substring, search_list, n=1, cutoff=0.1)
#             if match:
#                 close_matches.append(match[0])

#     # Return the closest match if any match is found
#     if close_matches:
#         return close_matches[0]
#     else:
#         return None

# def process_alignment(seqs, headers, seqs_contigs, header_c, output_bcr_path, num_jobs):
#     """
#     Processes sequence alignments and saves the processed data.

#     :param seqs: List of sequences to be aligned.
#     :param headers: List of headers corresponding to seqs.
#     :param seqs_contigs: List of contig sequences.
#     :param header_c: List of header sequences for contigs.
#     :param output_bcr_path: Path to save the output BCR data.
#     :param num_jobs: Number of parallel jobs to run.
#     """
#     best_alignments = find_best_alignments(seqs, seqs_contigs, num_jobs)
    
#     flattened_data = [
#         (A0, A1, A2[0], A2[1], A2[2]) for A0, A1, A2 in best_alignments
#     ]
    
#     df = pd.DataFrame(flattened_data, columns=["SCORE", "ALIGN_FRACTION", "IMGT_i", "CONTIG_i", "RC"])
#     filtered_df = df[(df['SCORE'] >= 200.0) & (df['ALIGN_FRACTION'] >= 0.95)]
    
#     unique_unitig_indices = filtered_df['CONTIG_i'].unique().tolist()
    
#     output_data = []
    
#     for unitig_index in unique_unitig_indices:
#         relevant_imgt_indices = filtered_df[filtered_df['CONTIG_i'] == unitig_index]['IMGT_i'].tolist()
        
#         collected_headers = [headers[imgt_index] for imgt_index in relevant_imgt_indices]
#         last_elements = [header[-1] for header in collected_headers]
#         most_have_none = last_elements.count('None') > len(last_elements) // 2
        
#         rc_status = filtered_df[filtered_df['CONTIG_i'] == unitig_index]['RC'].mode().iloc[0]
        
#         combined_header = ' '.join([string for sublist in collected_headers for string in sublist])
#         search_list = ["IgG", "IgA", "IgM", "IgE", "IgD", "IgK", "IgL"]
#         ig_name = find_most_similar(combined_header, search_list)
#         seq = seqs_contigs[unitig_index]
#         header_seq = header_c[unitig_index][0]
#         ig_name += ',  W: ' + header_seq
        
#         output_seq = str(Seq(seq).reverse_complement()) if most_have_none and rc_status == 'forward' else seq
#         output_data.append((ig_name, combined_header, output_seq))
    
#     # Sort based on the weight
#     output_data.sort(key=lambda x: float(x[0].split("W: ")[-1]), reverse=True)
    
#     # Save the sequences in the sorted order
#     for ig_name, combined_header, output_seq in output_data:
#         save_to_fasta([ig_name], [output_seq], output_bcr_path, include_header=True, append=True)

def process_alignment(seqs, headers, seqs_contigs, header_c, output_bcr_path, tag=False, variable=False, num_jobs=1):
    """
    Processes sequence alignments and saves the processed data.

    :param seqs: List of sequences to be aligned.
    :param headers: List of headers corresponding to seqs.
    :param seqs_contigs: List of contig sequences.
    :param header_c: List of header sequences for contigs.
    :param output_bcr_path: Path to save the output BCR data.
    :param tag: Boolean flag to determine if tags need to be added to headers.
    :param variable: Boolean flag to determine if the sequence is variable or heavy chain.
    :param num_jobs: Number of parallel jobs to run.
    """
    best_alignments = find_best_alignments(seqs, seqs_contigs, num_jobs)
    
    flattened_data = [
        (A0, A1, A2[0], A2[1], A2[2]) for A0, A1, A2 in best_alignments
    ]
    
    df = pd.DataFrame(flattened_data, columns=["SCORE", "ALIGN_FRACTION", "IMGT_i", "CONTIG_i", "RC"])

    filtered_df = df[(df['SCORE'] >= 200.0) & (df['ALIGN_FRACTION'] >= 0.90)]
    
    unique_unitig_indices = filtered_df['CONTIG_i'].unique().tolist()
    
    output_data = []
    
    for unitig_index in unique_unitig_indices:
        imgt_index = filtered_df[filtered_df['CONTIG_i'] == unitig_index]['IMGT_i'].tolist()[0]
        
        collected_headers = headers[imgt_index]
        most_have_none = collected_headers[-1] == 'None'
        
        rc_status = filtered_df[filtered_df['CONTIG_i'] == unitig_index]['RC'].mode().iloc[0]
        
        search_list = ["IgG", "IgA", "IgM", "IgE", "IgD", "IgK", "IgL"]

        matches = [SequenceMatcher(None, substring, ig.upper()).ratio() 
                   for substring in collected_headers if substring.startswith("IG") 
                   for ig in search_list]

        ind_max = np.argmax(matches)
        ig_name = search_list[ind_max]
        seq = seqs_contigs[unitig_index]
        header_seq = header_c[unitig_index][0]
        
        if tag:
            # Add either variable_chain or heavy_chain based on the variable parameter
            if variable:
                ig_name += ', W: {}, variable_region'.format(header_seq)
                output_seq = str(Seq(seq).reverse_complement()) if most_have_none and rc_status == 'reverse' else seq
            else:
                ig_name += ', W: {}, constant_region'.format(header_seq)
                output_seq = str(Seq(seq).reverse_complement()) if most_have_none and rc_status == 'forward' else seq     
        else:
            ig_name += ', W: {}'.format(header_seq)
            output_seq = str(Seq(seq).reverse_complement()) if most_have_none and rc_status == 'forward' else seq
        
        output_data.append((ig_name, output_seq))
    
    output_data.sort(key=lambda x: float(re.search(r'\d+\.\d+', x[0]).group()), reverse=True)
    
    
    # Save the sequences in the sorted order
    for ig_name, output_seq in output_data:
        save_to_fasta([ig_name], [output_seq], output_bcr_path, include_header=True, append=True)
        

def classify_reads(headers, sequences):
    bcr_heavy_chain = ["IgG", "IgA", "IgM", "IgE", "IgD"]
    bcr_light_chain = ["IgK", "IgL"]
    
    classified_reads = {
        'BCR_heavy_chain': {},
        'BCR_light_chain': {}
    }
    
    for header, sequence in zip(headers, sequences):
        # Example: header = ['>IgG,', 'W:', '4947.541521739131']
        read_type = header[0].split(",")[0]
        
        weight = float(header[2].split(",")[0])
        
        if read_type in bcr_heavy_chain:
            if read_type not in classified_reads['BCR_heavy_chain']:
                classified_reads['BCR_heavy_chain'][read_type] = []
            classified_reads['BCR_heavy_chain'][read_type].append((read_type, weight, sequence))
            
        elif read_type in bcr_light_chain:
            if read_type not in classified_reads['BCR_light_chain']:
                classified_reads['BCR_light_chain'][read_type] = []
            classified_reads['BCR_light_chain'][read_type].append((read_type, weight, sequence))
    
    return classified_reads


def local_alignment_BCR(seq1, seq2, match_score=3, mismatch_penalty=-10, gap_open_penalty=-10, gap_extend_penalty=-10):
    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_penalty
    aligner.open_gap_score = gap_open_penalty
    aligner.extend_gap_score = gap_extend_penalty

    # Perform the alignment
    alignments = aligner.align(seq1, seq2)
    
    # Extract the alignment with the highest score
    best_alignment = alignments[0]
    
    return best_alignment

def merge_sequences(seq1, seq2, weight_1, weight_2, K):
    """
    Merges two sequences based on local alignment of the overlapping regions.

    Parameters:
    - seq1: First sequence (string).
    - seq2: Second sequence (string).
    - K: Length of the overlapping region to consider.
    - weight_func: A function that calculates the weight of a sequence.
    
    Returns:
    - The merged contig sequence (string).
    """
    
    
    # def glue_contigs(seq1, seq2, aligned_pairs, K):
    #     seq1_inds, seq2_inds = aligned_pairs[0]
    #     seq_merged = seq1[:-(K - seq1_inds[1])] + seq2[seq2_inds[1]:]
    #     return seq_merged
    
    def glue_contigs(seq1, seq2, aligned_pairs, K):
        
        seq1_inds1, seq1_inds2 = aligned_pairs[0]
        seq2_inds1, seq2_inds2 = aligned_pairs[1]
        seq_merged = seq1[:(len(seq1) + (K - seq1_inds2))] + seq2[seq2_inds2:]
        return seq_merged    
    
    def calculate_length(t):
        array1, array2 = t
        length1 = array1[1] - array1[0]
        length2 = array2[1] - array2[0]
        return max(length1, length2)

    begin_1 = seq1[:K]
    end_1 = seq1[-K:]
    begin_2 = seq2[:K]
    end_2 = seq2[-K:]
    
    avg_w =  (weight_1 + weight_2)/2.
    
    # Perform best local alignments
    best_alignment_1 = local_alignment_BCR(end_1, begin_2)
    aligned_pairs_1 = list(zip(best_alignment_1.aligned[0], best_alignment_1.aligned[1]))
    
    best_tuple_1 = max(aligned_pairs_1, key=calculate_length)
    len_align_1 = calculate_length(best_tuple_1)
    
    

    best_alignment_2 = local_alignment_BCR(end_2, begin_1)
    aligned_pairs_2 = list(zip(best_alignment_2.aligned[0], best_alignment_2.aligned[1]))
    
    best_tuple_2 = max(aligned_pairs_2, key=calculate_length)
    len_align_2 = calculate_length(best_tuple_2)

    # Check alignment lengths and merge sequences
    if len_align_1 >= 15 or len_align_2 >= 15:
        if len_align_1 >= len_align_2:
            
            merged_contig = glue_contigs(seq1, seq2, best_tuple_1,  K)
            return len_align_1, avg_w, merged_contig
        else:
            merged_contig = glue_contigs(seq2, seq1, best_tuple_2, K)
            return len_align_2, avg_w, merged_contig
    else:
        return None, None, None


def process_bcr_sequences(output_bcr_path):
    headers, sequences = read_fasta(output_bcr_path)
    classified_reads = classify_reads(headers, sequences)
    
    # Lists to store the resulting headers and sequences
    result_headers = []
    result_sequences = []

    # Process both 'BCR_heavy_chain' and 'BCR_light_chain'
    for chain_type in ['BCR_heavy_chain', 'BCR_light_chain']:
        weights_for_sorting = []
        for immunoglobulin_type in classified_reads[chain_type].keys():
            immunoglobins = classified_reads[chain_type][immunoglobulin_type]

            if not immunoglobins:
                continue

            heavy_seqs = [t[2] for t in immunoglobins]
            heavy_ws = [t[1] for t in immunoglobins]

            new_merged = heavy_seqs[0]
            new_merged_w = heavy_ws[0]
            remove_nodes = [0]

            while True:
                merged_seqs = []
                for i in range(len(heavy_seqs)):
                    if i not in remove_nodes:
                        score, avg_weight, merged_seq = merge_sequences(new_merged, heavy_seqs[i], new_merged_w, heavy_ws[i], 50)
                        merged_seqs.append((i, score, avg_weight, merged_seq))

                filtered_merged_seqs = [(i, j, k, l) for i, j, k, l in merged_seqs if j is not None]
                if not filtered_merged_seqs:
                    break

                filtered_merged_seqs_sort = sorted(filtered_merged_seqs, key=lambda x: x[2], reverse=True)
                new_merged = filtered_merged_seqs_sort[0][3]
                new_merged_w = filtered_merged_seqs_sort[0][2]
                remove_nodes.append(filtered_merged_seqs_sort[0][0])

            # Append the final merged sequence with the header
            header = f"{chain_type} {immunoglobulin_type} W: {new_merged_w}"
            result_headers.append(header)
            result_sequences.append(new_merged)
            weights_for_sorting.append(new_merged_w)
            

    return result_headers, result_sequences

def display_custom_banner():
    banner = r"""
  ____     _     _____  _   _  _____  ____  
 / ___|   / \   |_   _|| | | || ____||  _ \ 
| |  _   / _ \    | |  | |_| ||  _|  | |_) |
| |_| | / ___ \   | |  |  _  || |___ |  _ < 
 \____|/_/   \_\  |_|  |_| |_||_____||_| \_\
                                             
    """
    version = "Version: 0.1.0"
    print(banner)
    print(version)
    print("\n")


def parse_arguments(package_path):
    parser = argparse.ArgumentParser(description='Process RNA sequences.')
    parser.add_argument('--seq', help='Path to the input fastq.gz file')
    parser.add_argument('--k', type=int, default=25, help='kmers size parameter, default is 25')
    parser.add_argument('--min_freq', type=int, default=5, help='Minimum frequency, default is 5')
    parser.add_argument('--output_dir', help='Output directory, default is the directory of the input_fastq file')
    parser.add_argument('--num_jobs', type=int, default=-1, help='Number of jobs, default is -1')
    parser.add_argument('--merge_bcrs', action='store_true', help="If set, merge BCR sequences")
    parser.add_argument('--TCR', action='store_true', help="If set, assemble TCR sequences")
    parser.add_argument('--mouse', action='store_true', help="If set, use mouse database, otherwise use human")
    parser.add_argument('--variable_chain', action='store_true', help="If set, output both variable and heavy chains if different")

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
        args.output_dir = os.path.dirname(os.path.abspath(args.seq))

    return args

def main():
    package_path = '/home/mojtaba/Sys_admin/gather'  # Update this to the correct package path
    args = parse_arguments(package_path)
    display_custom_banner()
    
    scratch_dir = args.output_dir
    input_rna_fastq = args.seq
    input_unitigs_name = os.path.basename(input_rna_fastq).replace('.gz', '.unitigs.fa')
    
    path_hc = os.path.join(args.db, 'hc')
    path_hv = os.path.join(args.db, 'hv')
    path_lc = os.path.join(args.db, 'lc')
    path_lv = os.path.join(args.db, 'lv')
    
    output_contigs_name = input_unitigs_name.replace('unitigs', 'contigs')
    
    input_bcalm = os.path.join(scratch_dir, input_rna_fastq)
    input_path = os.path.join(scratch_dir, input_unitigs_name)
    output_path = os.path.join(scratch_dir, output_contigs_name)
    k = args.k
    min_freq = args.min_freq
    
    os.chdir(scratch_dir)
    run_bcalm(input_bcalm, k + 1, min_freq)
    
    weights, overlap_indices, unitigs = parse_unitigs(input_path)
    
    unitigs_rc = [reverse_complement(i) for i in unitigs]
    
    begin, end = convert_kmers_to_array(unitigs, k)
    begin_rc, end_rc = convert_kmers_to_array(unitigs_rc, k)
    
    overlap_inds_list_0 = find_matching_indices(begin, begin, begin_rc, weights,
                                                args.num_jobs)
    overlap_inds_list_1 = find_matching_indices(end, end, end_rc, weights,
                                                args.num_jobs)
    
    rep_overlap_ind = overlap_inds_list_0 + overlap_inds_list_1
    
    overlap_indices_filtered = remove_low_w_tuples(overlap_indices, rep_overlap_ind)
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
    
    optimum_paths, optimum_paths_W = find_optimum_path(G, wcc_l, unitigs, weights)
    
    Contigs = [
        merge_overlaps(path, overlap_indices_filtered, unitigs, unitigs_rc, k)
        for path in tqdm(optimum_paths, desc='Making contigs from optimum paths', colour="#800080")
    ]
    Contigs += [unitigs[i] for i in wcc_s if len(unitigs[i]) > 300]
    
    optimum_paths_W += [weights[i] for i in wcc_s if len(unitigs[i]) > 300]
    
    save_to_fasta(headers=optimum_paths_W, n_sequences=Contigs, file_path=output_path, include_header=True)
    
    header_hc, seqs_hc = read_fasta(path_hc, imgt_fasta=True)
    header_lc, seqs_lc = read_fasta(path_lc, imgt_fasta=True)
    header_hv, seqs_hv = read_fasta(path_hv, imgt_fasta=True)
    header_lv, seqs_lv = read_fasta(path_lv, imgt_fasta=True)
    
    seqs_hv = [seq.upper() for seq in seqs_hv]
    seqs_lv = [seq.upper() for seq in seqs_lv]
    
    header_c, seqs_contigs = read_fasta(output_path)
    
    output_contigs_name = input_unitigs_name.replace('unitigs', 'BCR' if not args.TCR else 'TCR')
    output_bcr_path = os.path.join(scratch_dir, output_contigs_name)
    
    variable_chain_flag = args.variable_chain
    
    if variable_chain_flag:
        process_alignment(seqs_hc, header_hc, seqs_contigs, header_c, output_bcr_path, tag=variable_chain_flag, variable=False, num_jobs=args.num_jobs)
        process_alignment(seqs_hv, header_hv, seqs_contigs, header_c, output_bcr_path, tag=variable_chain_flag, variable=True, num_jobs=args.num_jobs)
        process_alignment(seqs_lc, header_lc, seqs_contigs, header_c, output_bcr_path, tag=variable_chain_flag, variable=False, num_jobs=args.num_jobs)
        process_alignment(seqs_lv, header_lv, seqs_contigs, header_c, output_bcr_path, tag=variable_chain_flag, variable=True, num_jobs=args.num_jobs)
    
    else:
        process_alignment(seqs_hc, header_hc, seqs_contigs, header_c, output_bcr_path, num_jobs=args.num_jobs)
        process_alignment(seqs_lc, header_lc, seqs_contigs, header_c, output_bcr_path, num_jobs=args.num_jobs)
        
        
    if args.merge_bcrs:
        headers, sequences =  process_bcr_sequences(output_bcr_path)
        
        output_highest_expressed_BCR = input_unitigs_name.replace('unitigs', 'BCR_HE' if not args.TCR else 'TCR_HE')

        save_to_fasta(headers=headers, n_sequences=sequences, file_path=output_highest_expressed_BCR, include_header=True, append=False)
    
if __name__ == '__main__':
    main()
    
