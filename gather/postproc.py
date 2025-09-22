#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 11:11:46 2024

@author:  Neuroimmunology UiO
"""

import os
import subprocess
import argparse
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import numpy as np
from gather.utils import Aligner, IO, ProcessSeq
import shutil
from pathlib import Path


def run_command(command, env=None):
    """Run a shell command."""
    subprocess.run(command, shell=True, check=True, env=env)
    
def process_bcr_heavy_chains_cds(
    df1,
    fasta_path, c_domains):

    io = IO()
    algn = Aligner()

    c_domains_seqs = [f'{domain}-seq' for domain in c_domains]

    # Read sequences and headers
    header_hc, seqs_hc = io.read_fasta(fasta_path, imgt_fasta=True)
    seqs_hc = [seq.upper() for seq in seqs_hc]

    # Extract sequences and c_call values from the dataframe
    gather_heavy_seqs = df1['sequence'].tolist()
    gather_heavy_ccalls = df1['c_call'].tolist()
    gather_heavy_ccalls = [
    x.split(',')[0] if isinstance(x, str) else (np.nan if np.isnan(x) else x)
    for x in gather_heavy_ccalls
    ]

    # Initialize new columns in df1 for C_domains_seqs
    for seq_col in c_domains_seqs:
        df1[seq_col] = None  # Initialize with None or a default value

    # Loop through each row; currently the script processes only the first, but this is a start
    for row_idx in range(len(df1)):
        gather_heavy_seqs_i = gather_heavy_seqs[row_idx]
        index_ccalls = [i for i, j in enumerate(header_hc) if gather_heavy_ccalls[row_idx] in j]
        ccalls_i = [j[4] for i, j in enumerate(header_hc) if gather_heavy_ccalls[row_idx] in j]

        best_alignments = algn.find_best_alignments(
            [gather_heavy_seqs_i], 
            [j for i, j in enumerate(seqs_hc) if i in index_ccalls], 
            spades=False, 
            c_domains=True, 
            num_jobs=8
        )

        # Prepare data to update DataFrame
        flattened_data = [
            (A0, A1, A2[0], A2[1], A2[2], A3[0], A3[1]) 
            for A0, A1, A2, A3 in best_alignments
        ]

        scores      = [i[0] for i in flattened_data]
        A_fraction   = [i[1] for i in flattened_data]
        A_cdomains   = [ccalls_i[i[3]] for i in flattened_data]
        A_cdomains_seq = [gather_heavy_seqs_i[i[-2]:i[-1]] for i in flattened_data]

        score_min = 100  # threshold for SCORE

        # Update DataFrame with A_fraction and A_cdomains_seq values
        for score, domain, frac, seq in zip(scores, A_cdomains, A_fraction, A_cdomains_seq):
            if domain in c_domains:
                frac_col_name = f'{domain}-identity'
                seq_col_name = f'{domain}-seq'

                if frac >= 0.75 and score >= score_min:
                    df1.at[row_idx, frac_col_name] = frac
                    df1.at[row_idx, seq_col_name] = seq
                else:
                    df1.at[row_idx, frac_col_name] = None
                    df1.at[row_idx, seq_col_name] = None

    return df1

def main():
    parser = argparse.ArgumentParser(description="Run IgBlastn, MakeDb, and additional processing on specified fasta file")

    # Try to find binaries in PATH
    igblastn_default = shutil.which('igblastn')
    igphyml_default = shutil.which('igphyml')

    if igblastn_default is None:
        print("Warning: 'igblastn' not found in your PATH. You should provide it using --igblastn_path.")
    if igphyml_default is None:
        print("Warning: 'igphyml' not found in your PATH. You should provide it using --igphyml_path.")

    parser.add_argument('--igblastn_path',
                        default=igblastn_default,
                        help='Path to the IgBlastn binary (default: found in PATH)')

    parser.add_argument('--igphyml_path',
                        default=igphyml_default,
                        help='Path to the IgPhyML binary (default: found in PATH)')

    parser.add_argument('--bcrs_dir', required=True, help='Directory where the BCR fasta files located')
    parser.add_argument('--chain', choices=['heavy', 'light'], help='Type of BCR chain: heavy or light or mixed')
    parser.add_argument('--output_dir', required=True, help='Directory where the output files will be saved')
    parser.add_argument('--clonality', action='store_true', help='Run the combined heavy and light chain analysis if both chains are present')
    parser.add_argument('--plot_output_file', default='clonal_plot.pdf', help='Filename for the clonal plot')
    parser.add_argument('--lineage_tree', action='store_true', help='Generate lineage trees using the provided R script')
    parser.add_argument('--mouse', action='store_true', help="If set, use mouse database, otherwise use human")
    parser.add_argument('--db_base', help='Custom base path to IgBlast database files (defaults to ~/share/igblast)')
    parser.add_argument('--germline_base', help='Custom base path to germline FASTA files (defaults to ~/share/germlines/imgt/<organism>)')

    args = parser.parse_args()

    # Final validation
    if not args.igblastn_path or not os.path.isfile(args.igblastn_path):
        print("\nERROR: IgBlastn binary not found. Please provide a valid path using --igblastn_path.")
        exit(1)

    if not args.igphyml_path or not os.path.isfile(args.igphyml_path):
        print("\nERROR: IgPhyML binary not found. Please provide a valid path using --igphyml_path.")
        exit(1)
    organism = "mouse" if args.mouse else "human"
    
    package_path = Path(__file__).resolve().parent.parent
    
    # db_base = Path(args.db_base).expanduser() if args.db_base else Path.home() / "share" / "igblast" 
    # germline_base = Path(args.germline_base).expanduser() if args.germline_base else Path.home() / "share" / "germlines" / "imgt" / organism
    
    # db_base = Path(args.db_base).expanduser() if args.db_base else Path("/share/igblast")
    # germline_base = Path(args.germline_base).expanduser() if args.germline_base else Path("/share/germlines/imgt") / organism
    
    # db_base
    if args.db_base:
        db_base = Path(args.db_base).expanduser()
    else:
        db_base = Path.home() / "share" / "igblast"
        if not db_base.exists():
            db_base = Path("/share/igblast")

    # germline_base
    if args.germline_base:
        germline_base = Path(args.germline_base).expanduser()
    else:
        germline_base = Path.home() / "share" / "germlines" / "imgt" / organism
        if not germline_base.exists():
            germline_base = Path("/share/germlines/imgt") / organism

    # Set IGDATA env
    os.environ["IGDATA"] = str(db_base)

    # IgBlast paths
    germline_db_V = db_base / "database" / f"imgt_{organism}_ig_v"
    germline_db_D = db_base / "database" / f"imgt_{organism}_ig_d"
    germline_db_J = db_base / "database" / f"imgt_{organism}_ig_j"
    germline_db_C = db_base / "database" / f"imgt_{organism}_ig_c"
    auxiliary_data = db_base / "optional_file" / f"{organism}_gl.aux"

    igblastn_path = Path(args.igblastn_path).expanduser()
    igphyml_path = Path(args.igphyml_path).expanduser()
    main_dir = Path(args.bcrs_dir).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    heavy_file_fa, light_file = ProcessSeq.collect_chains(main_dir, output_dir,
                   mouse=args.mouse,
                   heavy_name="collected_heavy_chains.fasta",
                   light_name="collected_light_chains.fasta")

    if args.chain == 'heavy':
        output_file = output_dir / "heavy_chains.fmt7"
        fasta_file = heavy_file_fa
    else:
        output_file = output_dir / "light_chains.fmt7"
        fasta_file = light_file

    # IgBlast command
    igblast_cmd = [
        str(igblastn_path),
        "-germline_db_V", str(germline_db_V),
        "-germline_db_D", str(germline_db_D),
        "-germline_db_J", str(germline_db_J),
        "-c_region_db", str(germline_db_C),
        "-auxiliary_data", str(auxiliary_data),
        "-domain_system", "imgt",
        "-ig_seqtype", "Ig",
        "-organism", organism,
        "-outfmt", "7 std qseq sseq btop",
        "-query", str(fasta_file),
        "-num_threads", "4",
        "-evalue", "0.05",
        "-out", str(output_file)
    ]

    try:
        result = subprocess.run(igblast_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during IgBlast: {e.stderr.decode('utf-8')}")
        return

    # Reference sequences based on organism
    ref_v = germline_base / "vdj" / f"imgt_{organism}_IGHV.fasta"
    ref_d = germline_base / "vdj" / f"imgt_{organism}_IGHD.fasta"
    ref_j = germline_base / "vdj" / f"imgt_{organism}_IGHJ.fasta"
    ref_lv = germline_base / "vdj" / f"imgt_{organism}_IGLV.fasta"
    ref_lj = germline_base / "vdj" / f"imgt_{organism}_IGLJ.fasta"
    ref_kv = germline_base / "vdj" / f"imgt_{organism}_IGKV.fasta"
    ref_kj = germline_base / "vdj" / f"imgt_{organism}_IGKJ.fasta"
    ref_kc = germline_base / "constant" / f"imgt_{organism}_IGKC.fasta"
    ref_lc = germline_base / "constant" / f"imgt_{organism}_IGLC.fasta"
    ref_hc = germline_base / "constant" / f"imgt_{organism}_IGHC.fasta"

    print(f"Using {organism} reference files from {germline_base}")

    # Determine references based on chain type
    if args.chain == 'heavy':
        germline_refs = [ref_v, ref_d, ref_j, ref_hc]
        output_name =  output_dir / 'heavy_chains_db-pass.tsv'
    elif args.chain == 'light':
        output_name =  output_dir / 'light_chains_db-pass.tsv'
        germline_refs = [ref_lv, ref_lj, ref_kv, ref_kj, ref_kc, ref_lc]
    else:  # mixed
        germline_refs = [ref_v, ref_d, ref_j, ref_hc, ref_lv, ref_lj, ref_kv, ref_kj, ref_kc, ref_lc]
        output_name =  output_dir / 'heavy_light_chains_db-pass.tsv'
       

    # Construct the MakeDb command
    make_db_cmd = [
    "MakeDb.py", "igblast",
    "-i", output_file,
    "-s", fasta_file,
    "-o", output_name,
    "-r"] + germline_refs + ["--extended", "--partial"]

    # Running the MakeDb command
    try:
        result = subprocess.run(make_db_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(f"Error occurred during MakeDb: {e.stderr.decode('utf-8')}")
        return
    
    if args.chain == 'heavy':
        db_pass_tsv_heavy = os.path.join(output_dir, "heavy_chains_db-pass.tsv")
        df1 = pd.read_csv(db_pass_tsv_heavy, sep='\t')
        c_domains = [
            'CH1', 'CH2', 'CH3', 'CH3-CHS', 'CH4-CHS', 'CHS', 
            'H', 'H-CH2', 'H1', 'H2', 'H3', 'H4',
            'M', 'M1', 'M2'
        ]
        if args.mouse:
            db_species = 'mouse'
        else:
            db_species = 'human'
        
        database = os.path.join(package_path, 'database', db_species, 'BCR', 'IGHC.fasta')
        df1 = process_bcr_heavy_chains_cds(df1, database, c_domains)
        identity_cols = [f"{r}-identity" for r in c_domains]

        available_identity_cols = [c for c in identity_cols if c in df1.columns]

        # Condition: any available identity < 1.0
        condition = (df1[available_identity_cols].lt(1.0)).any(axis=1)
        
        def add_suffix_first(val, suffix="-N"):
            if pd.isna(val):  # keep NaN as is
                return val
            parts = str(val).split(",", 1)  
            parts[0] = parts[0] + suffix    
            return ",".join(parts)         

        df1["c_novel"] = np.where(
            condition & df1["c_call"].notna(),
            df1["c_call"].apply(lambda x: add_suffix_first(x, "-N")),
            df1["c_call"]
        )
        df1.to_csv(db_pass_tsv_heavy, sep='\t', index=False)
            
    # If the combined analysis is requested, check and run the combined heavy and light chain analysis
    if args.clonality:
        if args.mouse:
            db_species = 'mouse'
        else:
            db_species = 'human'
        
        input_tsv_h = os.path.join(output_dir, "heavy_chains_db-pass.tsv")
        input_tsv_l = os.path.join(output_dir, "light_chains_db-pass.tsv")
        df1 = pd.read_csv(input_tsv_h, sep='\t')
        df2 = pd.read_csv(input_tsv_l, sep='\t')
        
        df1 = ProcessSeq.select_best_bcr(df1)
        df2 = ProcessSeq.select_best_bcr(df2)
        
        input_tsv_h_best = input_tsv_h.replace('pass', 'best')
        input_tsv_l_best = input_tsv_l.replace('pass', 'best')
        
        df1.to_csv(input_tsv_h_best, sep='\t', index=False)
        df2.to_csv(input_tsv_l_best, sep='\t', index=False)

        # Command and arguments for DefineClones
        define_clones_cmd = [
            "DefineClones.py", 
            "-d", input_tsv_h_best, 
            "--act", "set", 
            "--model", "ham", 
            "--norm", "len", 
            "--dist", "0.16"
        ]

        # Running the DefineClones command
        try:
            result = subprocess.run(define_clones_cmd, capture_output=True, text=True)
            # Print the output
            print("STDOUT:", result.stdout)
            print("STDERR:", result.stderr)
            if result.returncode == 0:
                print("DefineClones command executed successfully.")
            else:
                print(f"DefineClones command failed with return code {result.returncode}.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred during DefineClones: {e.stderr.decode('utf-8')}")
        
        
        db_pass_tsv_heavy = os.path.join(output_dir, "heavy_chains_db-best_clone-pass.tsv")
        db_pass_tsv_heavy_best = os.path.join(output_dir, "heavy_chains_db-best.tsv")
        ProcessSeq.harmonize_pass_with_best(db_pass_tsv_heavy, db_pass_tsv_heavy_best, db_pass_tsv_heavy)
        

        if os.path.exists(db_pass_tsv_heavy) and os.path.exists(input_tsv_l_best):
            # Define the list of colors
            colors = [
                "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02",
                "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3",
                "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8",
                "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
                "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"
            ]

            # Reading TSV files into DataFrames
            df1 = pd.read_csv(db_pass_tsv_heavy, sep='\t')
            df2 = pd.read_csv(input_tsv_l_best, sep='\t')

            
            # Create a dictionary to map cell_id to clone_id using the first DataFrame
            cell_id_to_clone_id = df1.set_index('sequence_id')['clone_id'].to_dict()

            # Use the dictionary to assign clone_id to df2 based on cell_id
            df2['clone_id'] = df2['sequence_id'].map(lambda x: int(cell_id_to_clone_id.get(x, -1)))
            
            # Merge the two DataFrames by concatenation
            merged_df = pd.concat([df1, df2], ignore_index=True)

            # Save the combined DataFrame to a new TSV file if needed
            merged_df.to_csv(os.path.join(output_dir, 'merged_PW.tsv'), sep='\t', index=False)
            
            c_calls = df1['c_novel']
            unique_c_calls = c_calls.unique()
            color_dict = {c_call: colors[i % len(colors)] for i, c_call in enumerate(unique_c_calls)}
            df1['color'] = df1['c_novel'].map(color_dict)

            # Save individual dataframes
            df1.to_csv(db_pass_tsv_heavy, sep='\t', index=False)
            
            def create_clone_plot(df_clone, color_dict, output_file):
                """Create and save a clone plot from the dataframe."""
                sequence_id = df_clone['sequence_id']
                c_calls = df_clone['c_novel'].tolist()

                G = nx.Graph()

                # Add nodes to the graph
                for sid in sequence_id:
                    G.add_node(sid)

                # Add edges between cells of the same clone
                clones = df_clone.groupby('clone_id')
                for clone, group in clones:
                    cells = group['sequence_id'].tolist()
                    for i in range(len(cells)):
                        for j in range(i + 1, len(cells)):
                            G.add_edge(cells[i], cells[j])

                # Check if LaTeX is installed
                latex_installed = shutil.which("latex") is not None

                if latex_installed:
                    # Use LaTeX rendering
                    rc('text', usetex=True)
                    rc('text.latex', preamble=r'\usepackage{amsmath}')
                    rcParams.update({
                        'font.size': 14,
                        'font.family': 'serif',
                        'text.usetex': True
                    })
                else:
                    # Fall back to mathtext
                    rcParams.update({
                        'font.size': 14,
                        'font.family': 'serif',
                        'mathtext.fontset': 'cm',
                        'mathtext.rm': 'serif',
                        'axes.unicode_minus': False,
                        'text.usetex': False
                    })

                plt.figure(figsize=(14, 10), dpi=200)

                # Layout for nodes
                pos = nx.spring_layout(G, k=0.30)

                # Node coloring
                default_color = "chartreuse"
                node_colors = [color_dict.get(sid, default_color) if not (type(sid) is float and np.isnan(sid)) else default_color for sid in c_calls]

                # Draw nodes and edges
                nx.draw_networkx_nodes(G, pos, node_size=100, node_color=node_colors, alpha=0.8, edgecolors='black')
                nx.draw_networkx_edges(G, pos, edgelist=G.edges(), edge_color='black', width=0.7)

                # Hide axis
                plt.axis('off')

                # Legend
                # Legend
                legend_labels = {v: k for k, v in color_dict.items()}
                handles = [
                    plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, label=legend_labels[color],
                               markersize=12, markerfacecolor=color, markeredgewidth=1, markeredgecolor='black')
                    for color in color_dict.values()
                    ]
                plt.legend(
                handles=handles,
                loc='center left',
                bbox_to_anchor=(1, 0.5),   # push legend to the right side
                ncol=1,
                frameon=False,
                fontsize=20,
                title_fontsize=18
                )

                # Save and close
                plt.savefig(output_file, bbox_inches="tight")
                plt.close()

            df_clone = pd.read_csv(db_pass_tsv_heavy, sep='\t')
            plot_output_file = os.path.join(output_dir, args.plot_output_file)
            create_clone_plot(df_clone, color_dict, plot_output_file)
                
            if args.lineage_tree:
                r_script_path = os.path.join(package_path, 'gather', 'generate_tree_plots.R')
                merged_pw_path = os.path.join(output_dir, 'merged_PW.tsv')
                germline_path = os.path.expanduser(germline_base / "vdj" )
                try:
                    subprocess.run([
                     "Rscript", r_script_path,
                     merged_pw_path,
                     db_pass_tsv_heavy,
                     germline_path,
                     output_dir,
                     igphyml_path
                     ], check=True)   
                except subprocess.CalledProcessError as e:
                    print(f"Error occurred while running R script for lineage trees: {e.stderr}")
                
                if os.path.exists(merged_pw_path):
                    os.remove(merged_pw_path)
        else:
            print("Both heavy and light chain data files must be available to run the combined analysis.")

if __name__ == "__main__":
    main()