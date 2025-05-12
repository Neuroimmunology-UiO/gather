#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:55:33 2025

@author: Neuroimmunology UiO

"""

import re
from joblib import Parallel, delayed
from tqdm import tqdm
from itertools import product
from Bio import Align, SeqIO
import os
import subprocess
import gzip
import shutil
from collections import defaultdict
from difflib import SequenceMatcher
import numpy as np
import pandas as pd
from Bio.Seq import Seq


class IO:
    
    @staticmethod
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

    @staticmethod
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

    @staticmethod
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
                    fasta_file.write(f">{seq_id}\n")

                fasta_file.write(sequence + "\n")

    @staticmethod
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

    @staticmethod
    def concatenate_fastq_files(file_1, file_2, output_file):
        """Concatenate two FASTQ.gz files into one gzipped file using Biopython."""
        with gzip.open(output_file, 'wt') as out:
            # Read and write records from the first compressed file
            with gzip.open(file_1, 'rt') as f1:
                for record in SeqIO.parse(f1, "fastq"):
                    SeqIO.write(record, out, "fastq")
            # Read and write records from the second compressed file if provided
            if file_2:
                with gzip.open(file_2, 'rt') as f2:
                    for record in SeqIO.parse(f2, "fastq"):
                        SeqIO.write(record, out, "fastq")


class CommandExecutor:
    
    def run_bcalm(self, input_file, kmer_size, abundance_min):
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
        if result.returncode != 0:
            print("Error running bcalm:")
            raise RuntimeError("BCALM command failed with return code {}".format(result.returncode))

    def run_spades(self, spades_path, output_dir, input_file_1=None, input_file_2=None, is_rna=False):
        """
        Run the SPAdes command with specified parameters.

        :param spades_path: Path to the SPAdes executable.
        :param output_dir: Directory where SPAdes output will be stored.
        :param input_file_1: Path to the first paired-end input file or the interlaced input file.
        :param input_file_2: Path to the second paired-end input file.
        :param is_rna: Boolean flag to specify whether to run SPAdes in RNA mode (default is False).
        :return: A boolean indicating success or failure, and an error message if any.
        """
        # Determine if the input is interlaced or paired-end based on provided files
        if input_file_1 and input_file_2:
            # Use paired-end files
            command = [spades_path, "-o", output_dir, "-1", input_file_1, "-2", input_file_2]
        elif input_file_1 and not input_file_2:
            # Use interlaced file
            command = [spades_path, "-o", output_dir, "-s", input_file_1]
        else:
            return False, "Invalid input: Provide either both input_file_1 and input_file_2, or a single interlaced input_file_1."

        if is_rna:
            command.append("--rna")

        # Run the command, suppress standard output, and capture standard error
        result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        # Check if the command was successful
        if result.returncode != 0:
            error_message = result.stderr.decode('utf-8')
            return False, f"SPAdes command failed with return code {result.returncode}. Error: {error_message}"

        return True, None

    def run_blast(command_type, blast_dir, output_dir="output", **kwargs):
        """
        Run the specified BLAST command with the given parameters.

        :param command_type: The type of command to run ('makeblastdb' or 'blastn').
        :param blast_dir: The directory where the BLAST binaries are located.
        :param output_dir: The directory where the output files should be saved.
        :param kwargs: The parameters to pass to the command.
        :raises RuntimeError: If the command fails.
        """
        # Ensure output directory exists
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if command_type == 'makeblastdb':
            output_db = os.path.join(output_dir, kwargs.get("output_db"))
            command = [
                f"{blast_dir}/makeblastdb",
                "-in", kwargs.get("input_file"),
                "-dbtype", kwargs.get("dbtype"),
                "-out", output_db
            ]

        elif command_type == 'blastn':
            output_file = os.path.join(output_dir, kwargs.get("output"))
            command = [
                f"{blast_dir}/blastn",
                "-query", kwargs.get("query"),
                "-db", kwargs.get("db"),
                "-out", output_file,
                "-outfmt", str(kwargs.get("outfmt")),
                "-num_threads", str(kwargs.get("num_threads")),
                "-perc_identity", str(kwargs.get("perc_identity"))
            ]

        else:
            raise ValueError("Unknown command type: {}".format(command_type))

        # Run the command and suppress all output
        result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Check if the command was successful
        if result.returncode != 0:
            print("Error running {}:".format(command_type))
            raise RuntimeError("{} command failed with return code {}".format(command_type, result.returncode))

    def run_single_cell_assembler(
    self,
    script_path,
    seq_merged,
    seq_1,
    spades_path,
    blast_dir, 
    use_spades=False,
    min_freq=1,
    num_jobs=8):
        """
        Run the single-cell assembler command with specified parameters.

        :param script_path: Path to the single-cell assembler script (e.g., sc_asm.py).
        :param seq_merged: Path to the merged sequence file.
        :param seq_1: Path to the first sequence file.
        :param spades_path: Path to SPAdes executable or "spades.py" if in PATH. Checks SPADES_BIN env var by default.
        :param blast_dir: Directory for BLAST executables or empty string if in PATH. Checks BLAST_BIN env var by default.
        :param use_spades: Boolean flag to specify whether to use SPAdes (default is False).
        :param min_freq: Minimum frequency parameter for the assembler.
        :param num_jobs: Number of jobs (threads) to use.
        :return: A boolean indicating success or failure, and an error message if any.
        """

        # Construct the command with the required parameters
        command = [
            "python", script_path,
            "--seq_merged", seq_merged,
            "--seq_1", seq_1,
            "--min_freq", str(min_freq),
            "--num_jobs", str(num_jobs),
            "--blast_dir", blast_dir,
        ]

        # Include SPAdes arguments only if `use_spades` is True
        if use_spades:
            command.extend(["--use_spades", "--spades_path", spades_path])

        # Run the command, suppress stdout, capture stderr
        result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

        # Check for errors
        if result.returncode != 0:
            error_message = result.stderr.decode('utf-8')
            return False, f"Single-cell assembler failed (code {result.returncode}): {error_message}"

        return True, None
    
   

class AlignmentError(Exception):
    """Exception raised for errors in the alignment process."""
    pass

class Aligner:

    @staticmethod
    def count_contiguous_pipes(match_line):
        """
        Count the maximum number of contiguous '|' characters in a match line.

        :param match_line: A string representing matches ('|') and mismatches ('*').
        :return: The maximum count of contiguous matches.
        """
        max_count = 0
        current_count = 0
        
        for char in match_line:
            if char == '|':
                current_count += 1
                if current_count > max_count:
                    max_count = current_count
            else:
                current_count = 0
        
        return max_count

    @staticmethod
    def format_alignment_pretty(alignment):
        """
        Format the alignment for pretty printing and calculate contiguous matches.

        :param alignment: A single alignment object.
        :return: The count of contiguous matches.
        """
        seq1 = alignment.sequences[0]
        seq2 = alignment.sequences[1]
        
        align1 = ""
        align2 = ""
        match_line = ""

        pos1, pos2 = 0, 0
        aligned_pairs = list(zip(alignment.aligned[0], alignment.aligned[1]))

        for (start1, end1), (start2, end2) in aligned_pairs:
            # Handle gaps in seq1
            if start1 > pos1:
                align1 += seq1[pos1:start1]
                align2 += '-' * (start1 - pos1)
                match_line += ' ' * (start1 - pos1)
            
            # Handle gaps in seq2
            if start2 > pos2:
                align1 += '-' * (start2 - pos2)
                align2 += seq2[pos2:start2]
                match_line += ' ' * (start2 - pos2)
            
            align1 += seq1[start1:end1]
            align2 += seq2[start2:end2]
    
            for i in range(end1 - start1):
                match_line += "|" if seq1[start1 + i] == seq2[start2 + i] else "*"

            pos1, pos2 = end1, end2

        contiguous_pipes_count = Aligner.count_contiguous_pipes(match_line)
        return contiguous_pipes_count

    @staticmethod
    def compute_aligned_fraction(identities, gaps, mismatches):
        """
        Compute the fraction of aligned characters.

        :param identities: Number of matching bases.
        :param gaps: Number of gaps in the alignment.
        :param mismatches: Number of mismatching bases.
        :return: Aligned fraction.
        """
        return identities / (identities + gaps + mismatches)

    @staticmethod
    def local_alignment_best(seq1, seq2, match_score=3, mismatch_penalty=-3, gap_open_penalty=-10, gap_extend_penalty=-1):
        """
        Perform local alignment between two sequences and return the best alignment.

        :param seq1: First sequence.
        :param seq2: Second sequence.
        :param match_score: Score for a match.
        :param mismatch_penalty: Penalty for a mismatch.
        :param gap_open_penalty: Penalty for opening a gap.
        :param gap_extend_penalty: Penalty for extending a gap.
        :return: Best alignment score and details.
        """
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_penalty
        aligner.open_gap_score = gap_open_penalty
        aligner.extend_gap_score = gap_extend_penalty

        alignments = aligner.align(seq1, seq2)

        if not alignments:
            raise AlignmentError("No alignments found between sequences.")

        best_alignment = alignments[0]
        aligned_counts = list(best_alignment.counts())
        gaps, identities, mismatches = aligned_counts
        
        return best_alignment.score, gaps, identities, mismatches, best_alignment

    @staticmethod
    def reverse_complement(seq):
        """Generate the reverse complement of a DNA sequence."""
        complement = str.maketrans('ATCG', 'TAGC')
        return seq.translate(complement)[::-1]

    @staticmethod
    def process_combination(i, j, seq1, seq2):
        """
        Process a combination of sequence indices and perform alignments.
        
        :param i: Index of the first sequence.
        :param j: Index of the second sequence.
        :param seq1: First sequence.
        :param seq2: Second sequence.
        :return: Alignment results including score and direction.
        """
        try:
            score_forward, gaps_forward, identities_forward, mismatches_forward, best_alignment = Aligner.local_alignment_best(seq1, seq2)
            score_reverse, gaps_reverse, identities_reverse, mismatches_reverse, best_alignment = Aligner.local_alignment_best(seq1, Aligner.reverse_complement(seq2))

            if score_forward > score_reverse:
                A_forward = Aligner.compute_aligned_fraction(identities_forward, gaps_forward, mismatches_forward)
                return score_forward, A_forward, (i, j, 'forward'), best_alignment
            else:
                A_reverse = Aligner.compute_aligned_fraction(identities_reverse, gaps_reverse, mismatches_reverse)
                return score_reverse, A_reverse, (i, j, 'reverse'), best_alignment
        except AlignmentError as e:
            print(f"AlignmentError for pair (seq1 index {i}, seq2 index {j}): {e}")
            return float('-inf'), 0.0, (i, j, None)

    @staticmethod
    def process_combination_c_calls(i, j, seq1, seq2):
        """
        Process a combination of sequence indices and perform alignments for c-domains.
        
        :param i: Index of the first sequence.
        :param j: Index of the second sequence.
        :param seq1: First sequence.
        :param seq2: Second sequence.
        :return: Alignment results with specific coordinates.
        """
        try:
            score_forward, gaps_forward, identities_forward, mismatches_forward, best_alignment_f = Aligner.local_alignment_best(seq1, seq2)
            score_reverse, gaps_reverse, identities_reverse, mismatches_reverse, best_alignment_r = Aligner.local_alignment_best(seq1, Aligner.reverse_complement(seq2))

            if score_forward > score_reverse:
                A_forward = identities_forward / len(seq2)
                seq1_cordinates_f = tuple(best_alignment_f.coordinates[0])
                
                return score_forward, A_forward, (i, j, 'forward'), seq1_cordinates_f
            else:
                A_reverse = identities_reverse / len(seq2)
                seq1_cordinates_r = tuple(best_alignment_r.coordinates[0])
                return score_reverse, A_reverse, (i, j, 'reverse'), seq1_cordinates_r
        except AlignmentError as e:
            print(f"AlignmentError for pair (seq1 index {i}, seq2 index {j}): {e}")
            return float('-inf'), 0.0, (i, j, None)

    @staticmethod
    def process_combination_spades(i, j, seq1, seq2):
        """
        Process a combination of sequence indices for SPAdes alignments.
        
        :param i: Index of the first sequence.
        :param j: Index of the second sequence.
        :param seq1: First sequence.
        :param seq2: Second sequence.
        :return: Alignment results and format result.
        """
        try:
            _, _, _, _, top_alignment = Aligner.local_alignment_best(seq1, seq2)

            A_forward = Aligner.format_alignment_pretty(top_alignment)
            return A_forward, (i, j, 'forward')
        except AlignmentError as e:
            print(f"AlignmentError for pair (seq1 index {i}, seq2 index {j}): {e}")
            return 0.0, (i, j, None)

    @staticmethod
    def find_best_alignments(seqs1, seqs2, spades=False, c_domains=False, num_jobs=1):
        """
        Find the best alignments between two sets of sequences in parallel.

        :param seqs1: List of sequences from the first set.
        :param seqs2: List of sequences from the second set.
        :param spades: Boolean for SPAdes-related alignment processing.
        :param c_domains: Boolean for C-domain alignment processing.
        :param num_jobs: Number of parallel jobs.
        :return: Sorted list of results based on alignment scores.
        """
        filtered_seqs2 = [s for s in seqs2 if len(s) > 0]

        all_combinations = list(product(range(len(seqs1)), range(len(filtered_seqs2))))
        total_combinations = len(seqs1) * len(filtered_seqs2)
        
        if spades:
            results = Parallel(n_jobs=num_jobs)(
                delayed(Aligner.process_combination_spades)(i, j, seqs1[i], filtered_seqs2[j]) 
                for i, j in tqdm(all_combinations, total=total_combinations,
                                 desc='Processing pairwise alignments',
                                 bar_format="{l_bar}{bar:20}{r_bar}", colour="#4CAF50")
            )
        elif c_domains:
            results = Parallel(n_jobs=num_jobs)(
                delayed(Aligner.process_combination_c_calls)(i, j, seqs1[i], filtered_seqs2[j]) 
                for i, j in tqdm(all_combinations, total=total_combinations,
                                 desc='Processing pairwise alignments',
                                 bar_format="{l_bar}{bar:20}{r_bar}", colour="#4CAF50")
            )
        else:
            results = Parallel(n_jobs=num_jobs)(
                delayed(Aligner.process_combination)(i, j, seqs1[i], filtered_seqs2[j]) 
                for i, j in tqdm(all_combinations, total=total_combinations,
                                 desc='Processing pairwise alignments',
                                 bar_format="{l_bar}{bar:20}{r_bar}", colour="#4CAF50")
            )

        sorted_results = sorted(results, key=lambda x: x[0], reverse=True)
        return sorted_results

# Define the mapping from DNA characters to float values
DNA_TO_FLOAT = {'A': 1.0, 'T': 2.0, 'C': 3.0, 'G': 4.0}

class DBG:

    @staticmethod
    def convert_kmers_to_array(list_of_strings, k):
        """
        Convert k-mers in each string of a list to 2D numpy arrays of float values.
        
        :param list_of_strings: List of DNA sequences as strings.
        :param k: Length of k-mers to extract from the beginning and end of each sequence.
        :return: Tuple of two 2D numpy arrays: beginnings and ends.
        """
        num_strings = len(list_of_strings)
        beginnings = np.zeros((num_strings, k), dtype=float)
        ends = np.zeros((num_strings, k), dtype=float)
        
        for i, dna_str in enumerate(list_of_strings):
            beginning_str = dna_str[:k]
            end_str = dna_str[-k:]
            beginnings[i] = [DNA_TO_FLOAT[ch] for ch in beginning_str]
            ends[i] = [DNA_TO_FLOAT[ch] for ch in end_str]
        
        return beginnings, ends

    @staticmethod
    def find_matching_indices(array1, array2, array2_rc, weights, num_jobs):
        """
        Finds the indices of rows that are equal in two given 2D NumPy arrays,
        including consideration of the reverse complement.

        Parameters:
        - array1 (np.ndarray): First 2D NumPy array.
        - array2 (np.ndarray): Second 2D NumPy array.
        - array2_rc (np.ndarray): Reverse complement of second 2D NumPy array.
        - weights (np.ndarray): 1D NumPy array of weights corresponding to rows of arrays.
        - num_jobs (int): Number of parallel jobs to run.
        
        Returns:
        - list of int: Indices of the rows that have the lower weight compared to their equivalent pairs.
        """
        def compare_arrays(_array1, _array2, _num_jobs):
            rows = _array1.shape[0]
            result_indices = []

            def compare_chunk(start, end):
                chunk_results = []
                for i in range(start, end):
                    comparison_result = np.all(_array1[i, np.newaxis] == _array2, axis=1)
                    matched_indices = np.argwhere(comparison_result).flatten()
                    for m_idx in matched_indices:
                        chunk_results.append((i, m_idx))
                return chunk_results

            chunk_size = (rows + _num_jobs - 1) // _num_jobs
            chunks = [(i, min(i + chunk_size, rows)) for i in range(0, rows, chunk_size)]

            parallel_results = Parallel(n_jobs=_num_jobs)(
                delayed(compare_chunk)(start, end) for start, end in tqdm(chunks, desc='Processing unitigs', colour="#4CAF50")
            )
            
            result_indices.extend(indice for sublist in parallel_results for indice in sublist)
            return result_indices

        equal_indices = compare_arrays(array1, array2, num_jobs)
        equal_indices_rc = compare_arrays(array1, array2_rc, num_jobs)
        
        equal_indices_low_w_norc = [i if weights[j] > weights[i] else j for i, j in equal_indices if i != j]
        equal_indices_low_w_rc = [i if weights[j] > weights[i] else j for i, j in equal_indices_rc if i != j]
        
        return equal_indices_low_w_norc + equal_indices_low_w_rc

    @staticmethod
    def find_overlap_indices(array1, array2, array1_rc, array2_rc, num_jobs):
        """
        Finds the indices of rows that overlap between two given 2D NumPy arrays or their reverse complements.

        Parameters:
        - array1 (np.ndarray): First 2D NumPy array.
        - array2 (np.ndarray): Second 2D NumPy array.
        - array1_rc (np.ndarray): Reverse of first 2D NumPy array.
        - array2_rc (np.ndarray): Reverse of second 2D NumPy array.
        - num_jobs (int): Number of parallel jobs to run.
        
        Returns:
        - list of tuples: Each tuple contains the indices of matching rows and their directions.
        """
        def compare_arrays(_array1, _array2, _num_jobs):
            rows = _array1.shape[0]
            result_indices = []

            def compare_chunk(start, end):
                chunk_results = []
                for i in range(start, end):
                    comparison_result = np.all(_array1[i, np.newaxis] == _array2, axis=1)
                    matched_indices = np.argwhere(comparison_result).flatten()
                    for m_idx in matched_indices:
                        chunk_results.append((i, m_idx))
                return chunk_results

            chunk_size = (rows + _num_jobs - 1) // _num_jobs
            chunks = [(i, min(i + chunk_size, rows)) for i in range(0, rows, chunk_size)]

            parallel_results = Parallel(n_jobs=_num_jobs)(
                delayed(compare_chunk)(start, end) for start, end in tqdm(chunks, desc='Processing unitigs', colour="#4CAF50", ncols=100)
            )
            
            result_indices.extend(indice for sublist in parallel_results for indice in sublist)
            return result_indices

        equal_indices = compare_arrays(array1, array2, num_jobs)
        equal_indices_rc = compare_arrays(array1, array2_rc, num_jobs)
        equal_indices_rc_1 = compare_arrays(array1_rc, array2, num_jobs)
        equal_indices_rc_2 = compare_arrays(array1_rc, array2_rc, num_jobs)

        overlap_1 = [(i, j, 'f', 'f') for i, j in equal_indices if i != j]
        overlap_2 = [(i, j, 'f', 'r') for i, j in equal_indices_rc if i != j]
        overlap_3 = [(i, j, 'r', 'f') for i, j in equal_indices_rc_1 if i != j]
        overlap_4 = [(i, j, 'r', 'r') for i, j in equal_indices_rc_2 if i != j]
        
        return overlap_1 + overlap_2 + overlap_3 + overlap_4

    @staticmethod
    def remove_low_w_tuples(tuples_list, indices_array):
        """
        Remove tuples from a list where the first or second element is present in a provided array.

        :param tuples_list: List of tuples, each containing (int, int, str, str).
        :param indices_array: A numpy array of indices.
        :return: List of filtered tuples.
        """
        dtype = [('first', int), ('second', int), ('third', 'U1'), ('fourth', 'U1')]
        structured_array = np.array(tuples_list, dtype=dtype)

        mask_first = np.isin(structured_array['first'], indices_array)
        mask_second = np.isin(structured_array['second'], indices_array)

        mask = ~mask_first & ~mask_second

        filtered_structured_array = structured_array[mask]

        return [
            (row['first'], row['second'], row['third'], row['fourth']) 
            for row in filtered_structured_array
        ]

    @staticmethod
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

    @staticmethod
    def find_optimum_path(G, wcc, unitigs, weights):
        """
        Find the optimum path in a subgraph of G defined by the weakly connected components (WCC).

        :param G: The input graph (networkx graph).
        :param wcc: A list containing WCC sets.
        :param unitigs: an array of unitigs.
        :param weights: an array of weights corresponding to unitigs.
        :return: List of nodes representing the optimum path and their weights. 
        """
        longest_paths = []
        long_w = []

        for cluster in tqdm(wcc, desc="Finding the optimum path in WCCs", colour="#808080", ncols=100):
            longest_path_size = len(cluster)
            max_path_avg_weight = 0
            longest_cluster_path = None

            for start_node in cluster:
                paths = DBG.dfs_all_paths(G, start_node)
                for path in paths:
                    path_avg_weight = np.sum(np.fromiter((weights[node] for node in path), dtype=float)) / longest_path_size
                    if path_avg_weight > max_path_avg_weight:
                        longest_cluster_path = path
                        max_path_avg_weight = path_avg_weight

            if longest_cluster_path is not None:
                longest_paths.append(longest_cluster_path)
                long_w.append(max_path_avg_weight)

        return longest_paths, long_w

    @staticmethod
    def merging_unitigs(overlap_info, unitigs_max, unitigs, unitigs_rc, k, first_index=False):
        """
        Merge two unitigs based on overlap information and their orientations.

        :param overlap_info: Tuple containing (ind_1, ind_2, rc_1, rc_2).
        :param unitigs_max: The current maximum unitig sequence being assembled.
        :param unitigs: List of unitig sequences.
        :param unitigs_rc: List of reverse complement unitig sequences.
        :param k: The length of the k-mer overlap.
        :param first_index: Boolean indicating if the first index should be merged.
        :return: The updated maximum unitig sequence.
        """
        ind_1, ind_2, rc_1, rc_2 = overlap_info
        
        if first_index:
            if rc_1 == 'f' and rc_2 == 'f':
                unitigs_max = unitigs[ind_1] + unitigs[ind_2][k:]
            elif rc_1 == 'f' and rc_2 == 'r':
                unitigs_max = unitigs[ind_1] + unitigs_rc[ind_2][k:]
            elif rc_1 == 'r' and rc_2 == 'f':
                unitigs_max = unitigs_rc[ind_1] + unitigs[ind_2][k:]
            else:
                unitigs_max = unitigs_rc[ind_1] + unitigs_rc[ind_2][k:]
        else:
            if rc_2 == 'f':
                unitigs_max += unitigs[ind_2][k:]
            else:
                unitigs_max += unitigs_rc[ind_2][k:]

        return unitigs_max

    @staticmethod
    def merge_overlaps(indices, overlap_inds_list, unitigs, unitigs_rc, k):
        """
        Merge multiple unitigs based on overlap information and given indices.

        :param indices: List of indices to merge.
        :param overlap_inds_list: List of tuples containing overlap information.
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

            mask = (filtered_oil[:, 0] == np.int64(current_index)) & (filtered_oil[:, 1] == np.int64(next_index))
            overlap_info = filtered_oil[mask][0]

            if (i-1) == 0 and unitigs_max is None:
                unitigs_max = DBG.merging_unitigs(overlap_info, unitigs_max, unitigs, unitigs_rc, k, first_index=True)
            else:
                unitigs_max = DBG.merging_unitigs(overlap_info, unitigs_max, unitigs, unitigs_rc, k, first_index=False)

        return unitigs_max

class ProcessSeq:

    @staticmethod
    def process_alignment(seqs, headers, seqs_contigs, header_c, output_bcr_path, variable=False, num_jobs=1):
        """Processes sequence alignments and saves the processed data."""
        best_alignments = Aligner.find_best_alignments(seqs, seqs_contigs, num_jobs)
        
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
            
            if variable:
                ig_name += ', W: {}, variable_region'.format(header_seq)
                output_seq = str(Seq(seq).reverse_complement()) if most_have_none and rc_status == 'reverse' else seq
            else:
                ig_name += ', W: {}, constant_region'.format(header_seq)
                output_seq = str(Seq(seq).reverse_complement()) if most_have_none and rc_status == 'forward' else seq
            
            output_data.append((ig_name, output_seq))
        
        output_data.sort(key=lambda x: float(re.search(r'\d+\.\d+', x[0]).group()), reverse=True)
        
        for ig_name, output_seq in output_data:
            IO.save_to_fasta([ig_name], [output_seq], output_bcr_path, include_header=True, append=True)

    @staticmethod
    def process_alignment_blast(seqs_contigs, output_bcr_path, output_bcr_contigus_path, blast_output_file, blast_dir,
                                path_imgt, scratch_dir, output_path, output_spades_path=None, spades=False,
                                seqs_contigs_spades=None, num_jobs=1):
        """
        Processes sequence alignments and saves the processed data.

        :param seqs_contigs: List of contig sequences.
        :param output_bcr_path: Path to save the output BCR data.
        :param blast_output_file: The input BLAST output file path.
        :param blast_dir: Directory for BLAST operations.
        :param path_imgt: Path to the IMGT reference.
        :param scratch_dir: The scratch directory.
        :param output_path: The path for output sequences.
        :param spades: A boolean to indicate whether to use SPAdes results.
        :param seqs_contigs_spades: List of SPAdes-contig sequences, only relevant if spades is True.
        :param num_jobs: Number of parallel jobs to run.
        """

        CommandExecutor.run_blast(
             'makeblastdb',
             blast_dir=blast_dir,
             input_file=path_imgt,
             dbtype="nucl",
             output_db="imgt_db",
             output_dir=scratch_dir
        )


        # Define sequence extraction functions
        def get_spades_seq(qseqid):
            contig_index = int(qseqid.split('_')[1]) - 1
            return seqs_contigs_spades[contig_index]

        def get_gather_seq(qseqid):
            contig_index = int(qseqid.split('_')[0])
            return seqs_contigs[contig_index]

        output_data = []
        if spades and seqs_contigs_spades is not None:
            CommandExecutor.run_blast(
                    'blastn',
                    blast_dir=blast_dir,
                    query=output_spades_path,
                    db="imgt_db",
                    output=blast_output_file,
                    outfmt="6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore ",
                    num_threads=num_jobs,
                    perc_identity=90.0,
                    
                )
            
            
            column_names = [
               "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
            
            blast_results = pd.read_csv(blast_output_file, sep='\t', names=column_names)

            # Apply filtering
            filtered_blast_results = blast_results[(blast_results['evalue'] < 0.01) & (blast_results['bitscore'] > 50)].copy()

            # Function to determine the sequence direction
            def determine_direction(row):
                return 'forward' if row['sstart'] < row['send'] else 'reverse'

            filtered_blast_results.loc[:, 'direction'] = filtered_blast_results.apply(determine_direction, axis=1)

            filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['sseqid'].apply(lambda x: x.split('|')[1][:4])
            # exceptions = ['IGHV', 'IGHD', 'IGHJ']  

            # filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['chain_type'].apply(
            #     lambda x: 'IGHC' if x.startswith('IGH') or  x.startswith('IGH') and x not in exceptions else x
            #     )

            # filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['chain_type'].apply(
                # lambda x: 'IGHC' if x.startswith('IGH') and x != 'IGHV' else x)

            grouped_results = (
                filtered_blast_results
                .sort_values(by=['bitscore', 'evalue'], ascending=[False, True])
                .groupby(['qseqid', 'chain_type'], as_index=False)
                .first()
            )

            grouped_results['chain_type_full'] = grouped_results['sseqid'].apply(lambda x: x.split('|')[1])

            aggregated = grouped_results.groupby('qseqid', as_index=False).agg({
                'sseqid': 'first',
                'pident': 'first',
                'length': 'first',
                'mismatch': 'first',
                'gapopen': 'first',
                'qstart': 'first',
                'qend': 'first',
                'sstart': 'first',
                'send': 'first',
                'evalue': 'first',
                'bitscore': 'first',
                'direction': 'first',
                'chain_type': 'first',
                'chain_type_full': lambda x: ','.join(x)
                })
            
            aggregated['spades_seq'] = aggregated['qseqid'].apply(get_spades_seq)
            
            single_chain_rows = aggregated[~aggregated['chain_type_full'].str.contains(',')]

            multi_chain_rows = aggregated[aggregated['chain_type_full'].str.contains(',')]

            shared_single_chain_types = {
                single for single in single_chain_rows['chain_type_full']
                if any(single in row['chain_type_full'].split(',') for _, row in multi_chain_rows.iterrows())
            }

            # Filter out rows with single chain types that are shared in multi-chain strings
            final_aggregated = aggregated[~aggregated['chain_type_full'].isin(shared_single_chain_types)].reset_index(drop=True)
            
            final_aggregated['spades_seq_length'] = final_aggregated['spades_seq'].apply(len)
            idx = final_aggregated.groupby('chain_type_full')['spades_seq_length'].idxmax()
            final_results_spades = final_aggregated.loc[idx].reset_index(drop=True)
            sorted_final_results = final_results_spades.sort_values(by='chain_type_full').reset_index(drop=True)
            
            for ind, row in sorted_final_results.iterrows():
                contig_header = row['qseqid']
                imgt_header = row['sseqid'].split('|')
                imgt_header = [item.strip() for item in imgt_header]

                ig_name = row['chain_type_full']
                rc_status = row['direction']
                contig_index = int(contig_header.split('_')[1]) - 1
                seq = seqs_contigs_spades[contig_index]
                header_seq = contig_header.split('_')[5]

                ig_name += ', W: {}'.format(header_seq)
                output_seq = str(Seq(seq).reverse_complement()) if rc_status == 'reverse' else seq

                output_data.append((ig_name, output_seq))
        
            for ig_name, output_seq in output_data:
                IO.save_to_fasta([ig_name], [output_seq], output_bcr_contigus_path, include_header=True, append=True)
        
        
        CommandExecutor.run_blast(
                'blastn',
                blast_dir=blast_dir,
                query=output_path,
                db="imgt_db",
                output=blast_output_file,
                outfmt="6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore ",
                num_threads=num_jobs,
                perc_identity=90.0,
                
            )
        
        
        column_names = [
           "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
        
        blast_results = pd.read_csv(blast_output_file, sep='\t', names=column_names)

        # Apply filtering
        filtered_blast_results = blast_results[(blast_results['evalue'] < 0.01) & (blast_results['bitscore'] > 50)].copy()

        # Function to determine the sequence direction
        def determine_direction(row):
            return 'forward' if row['sstart'] < row['send'] else 'reverse'

        filtered_blast_results.loc[:, 'direction'] = filtered_blast_results.apply(determine_direction, axis=1)

        filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['sseqid'].apply(lambda x: x.split('|')[1][:4])

        # exceptions = ['IGHV', 'IGHD', 'IGHJ']  

        # filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['chain_type'].apply(
            # lambda x: 'IGHC' if x.startswith('IGH') or  x.startswith('IGH') and x not in exceptions else x
            # )

        grouped_results = (
            filtered_blast_results
            .sort_values(by=['bitscore', 'evalue'], ascending=[False, True])
            .groupby(['qseqid', 'chain_type'], as_index=False)
            .first()
        )

        grouped_results['chain_type_full'] = grouped_results['sseqid'].apply(lambda x: x.split('|')[1])

        aggregated = grouped_results.groupby('qseqid', as_index=False).agg({
            'sseqid': 'first',
            'pident': 'first',
            'length': 'first',
            'mismatch': 'first',
            'gapopen': 'first',
            'qstart': 'first',
            'qend': 'first',
            'sstart': 'first',
            'send': 'first',
            'evalue': 'first',
            'bitscore': 'first',
            'direction': 'first',
            'chain_type': 'first',
            'chain_type_full': lambda x: ','.join(x)
            })
        output_data = []
        aggregated['gather_seq'] = aggregated['qseqid'].apply(get_gather_seq)
        aggregated['gather_seq_length'] = aggregated['gather_seq'].apply(len)
        idx = aggregated.groupby('chain_type_full')['gather_seq_length'].idxmax()
        final_results_gather = aggregated.loc[idx].reset_index(drop=True)
        sorted_final_results = final_results_gather.sort_values(by='chain_type_full').reset_index(drop=True)

        for ind, row in sorted_final_results.iterrows():
            contig_header = row['qseqid']
            imgt_header = row['sseqid'].split('|')
            imgt_header = [item.strip() for item in imgt_header]

            ig_name = row['chain_type_full']
            rc_status = row['direction']

            contig_index = int(contig_header.split('_')[0])
            seq = seqs_contigs[contig_index]
            header_seq = contig_header.split('_')[1]

            ig_name += ', W: {}'.format(header_seq)
            output_seq = str(Seq(seq).reverse_complement()) if rc_status == 'reverse' else seq

            output_data.append((ig_name, output_seq))
        
        # Save the sequences in the sorted order
        for ig_name, output_seq in output_data:
            IO.save_to_fasta([ig_name], [output_seq], output_bcr_path, include_header=True, append=True)
        
    @staticmethod
    def process_alignment_blast_10x(seqs_contigs, output_bcr_path, output_bcr_contigus_path, blast_output_file, blast_dir,
                                path_imgt, scratch_dir, output_path, output_spades_path=None, spades=False,
                                seqs_contigs_spades=None, num_jobs=1):
        """
        Processes sequence alignments and saves the processed data.

        :param seqs_contigs: List of contig sequences.
        :param output_bcr_path: Path to save the output BCR data.
        :param blast_output_file: The input BLAST output file path.
        :param blast_dir: Directory for BLAST operations.
        :param path_imgt: Path to the IMGT reference.
        :param scratch_dir: The scratch directory.
        :param output_path: The path for output sequences.
        :param spades: A boolean to indicate whether to use SPAdes results.
        :param seqs_contigs_spades: List of SPAdes-contig sequences, only relevant if spades is True.
        :param num_jobs: Number of parallel jobs to run.
        """

        CommandExecutor.run_blast(
             'makeblastdb',
             blast_dir=blast_dir,
             input_file=path_imgt,
             dbtype="nucl",
             output_db="imgt_db",
             output_dir=scratch_dir
        )


        # Define sequence extraction functions
        def get_spades_seq(qseqid):
            contig_index = int(qseqid.split('_')[1]) - 1
            return seqs_contigs_spades[contig_index]

        def get_gather_seq(qseqid):
            contig_index = int(qseqid.split('_')[0])
            return seqs_contigs[contig_index]

        output_data = []
        if spades and seqs_contigs_spades is not None:
            CommandExecutor.run_blast(
                    'blastn',
                    blast_dir=blast_dir,
                    query=output_spades_path,
                    db="imgt_db",
                    output=blast_output_file,
                    outfmt="6 qseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore ",
                    num_threads=num_jobs,
                    perc_identity=90.0,
                    
                )
            
            
            column_names = [
               "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore" ]
            
            blast_results = pd.read_csv(blast_output_file, sep='\t', names=column_names)

            # Apply filtering
            filtered_blast_results = blast_results[(blast_results['evalue'] < 0.01) & (blast_results['bitscore'] > 50)].copy()

            # Function to determine the sequence direction
            def determine_direction(row):
                return 'forward' if row['sstart'] < row['send'] else 'reverse'

            filtered_blast_results.loc[:, 'direction'] = filtered_blast_results.apply(determine_direction, axis=1)

            filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['sseqid'].apply(lambda x: x.split('|')[1][:4])

            exceptions = ['IGHV', 'IGHD', 'IGHJ']  

            filtered_blast_results.loc[:, 'chain_type'] = filtered_blast_results['chain_type'].apply(
                lambda x: 'IGHC' if x.startswith('IGH') and x not in exceptions else x
                )

            grouped_results = (
                filtered_blast_results
                .sort_values(by=['bitscore', 'evalue'], ascending=[False, True])
                .groupby(['qseqid', 'chain_type'], as_index=False)
                .first()
            )

            grouped_results['chain_type_full'] = grouped_results['sseqid'].apply(lambda x: x.split('|')[1])

            aggregated = grouped_results.groupby('qseqid', as_index=False).agg({
                'sseqid': 'first',
                'pident': 'first',
                'length': 'first',
                'mismatch': 'first',
                'gapopen': 'first',
                'qstart': 'first',
                'qend': 'first',
                'sstart': 'first',
                'send': 'first',
                'evalue': 'first',
                'bitscore': 'first',
                'direction': 'first',
                'chain_type': 'first',
                'chain_type_full': lambda x: ','.join(x)
                })
            
            aggregated['spades_seq'] = aggregated['qseqid'].apply(get_spades_seq)
            
            single_chain_rows = aggregated[~aggregated['chain_type_full'].str.contains(',')]

            multi_chain_rows = aggregated[aggregated['chain_type_full'].str.contains(',')]

            shared_single_chain_types = {
                single for single in single_chain_rows['chain_type_full']
                if any(single in row['chain_type_full'].split(',') for _, row in multi_chain_rows.iterrows())
            }

            # Filter out rows with single chain types that are shared in multi-chain strings
            final_aggregated = aggregated[~aggregated['chain_type_full'].isin(shared_single_chain_types)].reset_index(drop=True)
            
            final_aggregated['spades_seq_length'] = final_aggregated['spades_seq'].apply(len)
            idx = final_aggregated.groupby('chain_type_full')['spades_seq_length'].idxmax()
            final_results_spades = final_aggregated.loc[idx].reset_index(drop=True)
            sorted_final_results = final_results_spades.sort_values(by='chain_type_full').reset_index(drop=True)
            
            for ind, row in sorted_final_results.iterrows():
                contig_header = row['qseqid']
                imgt_header = row['sseqid'].split('|')
                imgt_header = [item.strip() for item in imgt_header]

                ig_name = row['chain_type_full']
                rc_status = row['direction']
                contig_index = int(contig_header.split('_')[1]) - 1
                seq = seqs_contigs_spades[contig_index]
                header_seq = contig_header.split('_')[5]

                ig_name += ', W: {}'.format(header_seq)
                output_seq = str(Seq(seq).reverse_complement()) if rc_status == 'reverse' else seq

                output_data.append((ig_name, output_seq))
        
            for ig_name, output_seq in output_data:
                IO.save_to_fasta([ig_name], [output_seq], output_bcr_contigus_path, include_header=True, append=True)
        

    @staticmethod
    def classify_reads(headers, sequences):
        """Classify sequences into BCR heavy and light chains."""
        bcr_heavy_chain = ["IgG", "IgA", "IgM", "IgE", "IgD"]
        bcr_light_chain = ["IgK", "IgL"]
        search_list = ["IgG", "IgA", "IgM", "IgE", "IgD", "IgK", "IgL"]
        
        classified_reads = {
            'BCR_heavy_chain': {},
            'BCR_light_chain': {}
        }
        
        for header, sequence in zip(headers, sequences):
            matches = [SequenceMatcher(None, header[0][:4], ig.upper()).ratio() for ig in search_list]
            ind_max = np.argmax(matches)
            read_type = search_list[ind_max]
            ig_type = header[0].split(",")[0]
            region = header[-1]
            weight = float(header[2].split(",")[0])
            
            if read_type in bcr_heavy_chain:
                if region not in classified_reads['BCR_heavy_chain']:
                    classified_reads['BCR_heavy_chain'][region] = []
                classified_reads['BCR_heavy_chain'][region].append((ig_type, weight, sequence))
                
            elif read_type in bcr_light_chain:
                if region not in classified_reads['BCR_light_chain']:
                    classified_reads['BCR_light_chain'][region] = []
                classified_reads['BCR_light_chain'][region].append((ig_type, weight, sequence))
        
        return classified_reads

    @staticmethod
    def local_alignment_BCR(seq1, seq2, match_score=3, mismatch_penalty=-10, gap_open_penalty=-10, gap_extend_penalty=-10):
        """Perform a local alignment on two sequences."""
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = match_score
        aligner.mismatch_score = mismatch_penalty
        aligner.open_gap_score = gap_open_penalty
        aligner.extend_gap_score = gap_extend_penalty

        alignments = aligner.align(seq1, seq2)
        best_alignment = alignments[0]
        
        return best_alignment

    @staticmethod
    def merge_sequences(seq1, seq2, weight_1, weight_2, K):
        """Merge two sequences based on local alignment of overlapping regions."""
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

        avg_w = (weight_1 + weight_2) / 2.

        # Perform best local alignments
        best_alignment_1 = ProcessSeq.local_alignment_BCR(end_1, begin_2)
        aligned_pairs_1 = list(zip(best_alignment_1.aligned[0], best_alignment_1.aligned[1]))
        
        best_tuple_1 = max(aligned_pairs_1, key=calculate_length)
        len_align_1 = calculate_length(best_tuple_1)

        best_alignment_2 = ProcessSeq.local_alignment_BCR(end_2, begin_1)
        aligned_pairs_2 = list(zip(best_alignment_2.aligned[0], best_alignment_2.aligned[1]))
        
        best_tuple_2 = max(aligned_pairs_2, key=calculate_length)
        len_align_2 = calculate_length(best_tuple_2)

        if len_align_1 >= 20 or len_align_2 >= 20:
            if len_align_1 >= len_align_2:
                merged_contig = glue_contigs(seq1, seq2, best_tuple_1,  K)
                return len_align_1, avg_w, merged_contig
            else:
                merged_contig = glue_contigs(seq2, seq1, best_tuple_2, K)
                return len_align_2, avg_w, merged_contig
        else:
            return None, None, None

    @staticmethod
    def process_bcr_sequences_merge(output_bcr_path):
        """Process BCR sequences by merging."""
        headers, sequences = IO.read_fasta(output_bcr_path)
        classified_reads = ProcessSeq.classify_reads(headers, sequences)

        result_headers = []
        result_sequences = []
        weights_for_sorting = []
        unique_sequences = {}

        chain_types = ['BCR_heavy_chain', 'BCR_light_chain']
        regions = ['constant_region', 'variable_region']

        for chain_type in chain_types:
            constant_immunoglobins = classified_reads.get(chain_type, {}).get('constant_region', [])
            variable_immunoglobins = classified_reads.get(chain_type, {}).get('variable_region', [])
            
            if not constant_immunoglobins and not variable_immunoglobins:
                continue

            variable_seqs_set = {t[2] for t in variable_immunoglobins} if variable_immunoglobins else set()

            for region in regions:
                if region == 'constant_region':
                    immunoglobins = constant_immunoglobins
                else:
                    immunoglobins = variable_immunoglobins

                if not immunoglobins:
                    continue

                identical_idx = None

                if region == 'constant_region' and variable_seqs_set:
                    for idx, t in enumerate(immunoglobins):
                        if t[2] in variable_seqs_set:
                            identical_idx = idx
                            break
                elif region == 'variable_region' and constant_immunoglobins:
                    for idx, t in enumerate(immunoglobins):
                        if t[2] in {s[2] for s in constant_immunoglobins}:
                            identical_idx = idx
                            break

                if identical_idx is not None:
                    new_merged = immunoglobins[identical_idx][2]
                    new_merged_w = immunoglobins[identical_idx][1]
                    immunoglobulin_type = immunoglobins[identical_idx][0]
                    remove_nodes = [identical_idx]
                else:
                    new_merged = immunoglobins[0][2]
                    new_merged_w = immunoglobins[0][1]
                    immunoglobulin_type = immunoglobins[0][0]
                    remove_nodes = [0]

                while True:
                    merged_seqs = []
                    for i, (seq, weight) in enumerate([(t[2], t[1]) for t in immunoglobins]):
                        if i not in remove_nodes:
                            score, avg_weight, merged_seq = ProcessSeq.merge_sequences(new_merged, seq, new_merged_w, weight, 50)
                            if score is not None:
                                merged_seqs.append((i, score, avg_weight, merged_seq))

                    if not merged_seqs:
                        break

                    merged_seqs.sort(key=lambda x: x[2], reverse=True)
                    new_merged = merged_seqs[0][3]
                    new_merged_w = merged_seqs[0][2]
                    remove_nodes.append(merged_seqs[0][0])

                if new_merged not in unique_sequences:
                    header = f"{chain_type} {region} {immunoglobulin_type} W: {new_merged_w}"
                    result_headers.append(header)
                    result_sequences.append(new_merged)
                    weights_for_sorting.append(new_merged_w)
                    unique_sequences[new_merged] = new_merged_w

        return result_headers, result_sequences

    @staticmethod
    def process_bcr_sequences(output_bcr_path):
        """Process BCR sequences and classify them into headers and regions."""
        def find_identical_indices_np(list_a, list_b):
            array_a = np.array(list_a)
            array_b = np.array(list_b)
            comparison_matrix = array_a[:, None] == array_b
            identical_indices = np.argwhere(comparison_matrix)
            return list(map(tuple, identical_indices))

        def create_header(chain_type, immunoglobulin_type, w_1, region=None):
            if region:
                return f"{chain_type} {region} {immunoglobulin_type} W: {w_1}"
            return f"{chain_type} {immunoglobulin_type} W: {w_1}"

        try:
            headers, sequences = IO.read_fasta(output_bcr_path)
        except Exception as e:
            raise ValueError(f"Error reading FASTA data: {e}")

        try:
            classified_reads = ProcessSeq.classify_reads(headers, sequences)
        except Exception as e:
            raise ValueError(f"Error classifying reads: {e}")

        result_headers = []
        result_sequences = []

        chain_types = ['BCR_heavy_chain', 'BCR_light_chain']

        for chain_type in chain_types:
            reads = classified_reads.get(chain_type, {})
            
            constant_immunoglobins = reads.get('constant_region', [])
            variable_immunoglobins = reads.get('variable_region', [])
            
            if not constant_immunoglobins and not variable_immunoglobins:
                print(f"No constant or variable regions found for {chain_type}.")
                continue

            variable_seqs = [t[2] for t in variable_immunoglobins]
            constant_seqs = [t[2] for t in constant_immunoglobins]

            identical_seqs_inds = find_identical_indices_np(constant_seqs, variable_seqs)
            
            if identical_seqs_inds:
                if len(identical_seqs_inds) == 1 and identical_seqs_inds[0][0] != 0:
                    idx_c = identical_seqs_inds[0][0]
                    immunoglobulin_type = constant_immunoglobins[idx_c][0]
                    w_1 = constant_immunoglobins[idx_c][1]
                    header = create_header(chain_type, immunoglobulin_type, w_1)
                    result_headers.append(header)
                    result_sequences.append(constant_immunoglobins[idx_c][2])

                    immunoglobulin_type = constant_immunoglobins[0][0]
                    w_1 = constant_immunoglobins[0][1]
                    header = create_header(chain_type, immunoglobulin_type, w_1, 'constant')
                    result_headers.append(header)
                    result_sequences.append(constant_immunoglobins[0][2])
                else:
                    for idx in range(min(2, len(identical_seqs_inds))): 
                        idx_c = identical_seqs_inds[idx][0]
                        immunoglobulin_type = constant_immunoglobins[idx_c][0]
                        w_1 = constant_immunoglobins[idx_c][1]
                        header = create_header(chain_type, immunoglobulin_type, w_1)
                        result_headers.append(header)
                        result_sequences.append(constant_immunoglobins[idx_c][2])
            
            else:
                for region, immunoglobulins in {'constant': constant_immunoglobins, 'variable': variable_immunoglobins}.items():
                    if immunoglobulins:
                        immunoglobulin_type = immunoglobulins[0][0]
                        w_1 = immunoglobulins[0][1]
                        header = create_header(chain_type, immunoglobulin_type, w_1, region)
                        result_headers.append(header)
                        result_sequences.append(immunoglobulins[0][2])

        return result_headers, result_sequences

    @staticmethod
    def process_bcr_sequences_spades(output_bcr_path, num_jobs):
        """Process BCR sequences from SPAdes."""
        def find_identical_indices_np(list_a, list_b):
            array_a = np.array(list_a)
            array_b = np.array(list_b)
            comparison_matrix = array_a[:, None] == array_b
            identical_indices = np.argwhere(comparison_matrix)
            return list(map(tuple, identical_indices))
        
        def get_connected_components(match_constants_align, score_threshold=50):
            graph = defaultdict(set)
            all_nodes = set()
        
            for score, (i, j, direction) in match_constants_align:
                all_nodes.update([i, j])
                if i != j and score > score_threshold:
                    graph[i].add(j)
                    graph[j].add(i)

            def explore(node, visited, component):
                visited.add(node)
                component.add(node)
                for neighbor in graph[node]:
                    if neighbor not in visited:
                        explore(neighbor, visited, component)
        
            visited = set()
            components = []

            for node in all_nodes:
                if node not in visited:
                    component = set()
                    explore(node, visited, component)
                    components.append(component)

            return components

        def create_header(chain_type, immunoglobulin_type, region=None):
            if region:
                return f"{chain_type} {region} {immunoglobulin_type}"
            return f"{chain_type} {immunoglobulin_type}"

        try:
            headers, sequences = IO.read_fasta(output_bcr_path)
        except Exception as e:
            raise ValueError(f"Error reading FASTA data: {e}")

        try:
            classified_reads = ProcessSeq.classify_reads(headers, sequences)
        except Exception as e:
            raise ValueError(f"Error classifying reads: {e}")

        result_headers = []
        result_sequences = []

        chain_types = ['BCR_heavy_chain', 'BCR_light_chain']

        for chain_type in chain_types:
            reads = classified_reads.get(chain_type, {})

            constant_immunoglobins = reads.get('constant_region', [])
            variable_immunoglobins = reads.get('variable_region', [])
            
            combined_immunoglobins = constant_immunoglobins + variable_immunoglobins

            constant_seqs = [t[2] for t in combined_immunoglobins]
            
            match_constants_align = Aligner.find_best_alignments(constant_seqs, constant_seqs, spades=True, num_jobs=num_jobs)
            components = get_connected_components(match_constants_align)
            
            result = [tuple(sorted(component)) for component in components]
            for chain_group in result:
                longest_in_result = np.argmax([len(constant_seqs[i]) for i in chain_group])
                index_longest_seq = chain_group[longest_in_result]
            
                immunoglobulin_type = combined_immunoglobins[index_longest_seq][0]
                header = create_header(chain_type, immunoglobulin_type)
                result_headers.append(header)
                result_sequences.append(combined_immunoglobins[index_longest_seq][2])
            
        return result_headers, result_sequences

    @staticmethod
    def process_bcr_file(path_to_bcr):
        """Processes a BCR file to get sorted B-cell receptor sequences."""
        headers_BCR, seqs_BCR = IO.read_fasta(path_to_bcr, imgt_fasta=False)

        filtered_H = [
            (headers_BCR[i], seqs_BCR[i])
            for i in range(len(headers_BCR))
            if 'IGHV' in headers_BCR[i][0] and len(headers_BCR[i][0].split(',')) > 2
        ]

        weights_H = [
            float(headers_BCR[i][-1])
            for i in range(len(headers_BCR))
            if 'IGHV' in headers_BCR[i][0] and len(headers_BCR[i][0].split(',')) > 2
        ]
        
        if not filtered_H:
            filtered_H = [
                (headers_BCR[i], seqs_BCR[i])
                for i in range(len(headers_BCR))
                if 'IGH' in headers_BCR[i][0]
            ]
            weights_H = [
                float(headers_BCR[i][-1])
                for i in range(len(headers_BCR))
                if 'IGH' in headers_BCR[i][0]
            ]

        filtered_L = [
            (headers_BCR[i], seqs_BCR[i])
            for i in range(len(headers_BCR))
            if (
                ('IGKV' in headers_BCR[i][0] or 'IGLV' in headers_BCR[i][0])
                and len(headers_BCR[i][0].split(',')) > 2
            )
        ]

        weights_L = [
            float(headers_BCR[i][-1])
            for i in range(len(headers_BCR))
            if (
                ('IGKV' in headers_BCR[i][0] or 'IGLV' in headers_BCR[i][0])
                and len(headers_BCR[i][0].split(',')) > 2
            )
        ]

        if not filtered_L:
            filtered_L = [
                (headers_BCR[i], seqs_BCR[i])
                for i in range(len(headers_BCR))
                if 'IGK' in headers_BCR[i][0] or 'IGL' in headers_BCR[i][0]
            ]
            weights_L = [
                float(headers_BCR[i][-1])
                for i in range(len(headers_BCR))
                if 'IGK' in headers_BCR[i][0] or 'IGL' in headers_BCR[i][0]
            ]

        combined_H = sorted(zip(weights_H, filtered_H), key=lambda x: x[0], reverse=True)
        combined_L = sorted(zip(weights_L, filtered_L), key=lambda x: x[0], reverse=True)

        top_combined = combined_H[:2] + combined_L[:2]

        sorted_sequences = [seq for _, (_, seq) in top_combined]
        sorted_headers = [' '.join(header) for _, (header, _) in top_combined]

        return sorted_headers, sorted_sequences

    @staticmethod
    def clean_scratch_dir(scratch_dir, new_items):
        """Cleans up the specified directory except for specific file types."""
        allowed_extensions = {'.fa', '.fasta', '.fq', '.fastq', '.gz'}
        for item in new_items:
            item_path = os.path.join(scratch_dir, item)
            if os.path.isfile(item_path):
                if not item.endswith(tuple(allowed_extensions)):
                    os.remove(item_path)
            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)
