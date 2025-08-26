#!/usr/bin/env python3
"""
Enhanced GraphBAN Drug Screening Pipeline with Multi-UniProt Support
Usage: python GraphBAN_screening_pipeline.py --uniprot_ids P04637,Q9Y6Y9 --zinc_data zinc250k.csv
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd
from typing import List, Tuple
from multiprocessing import Pool
from functools import partial

def parse_args():
    parser = argparse.ArgumentParser(
        description="Run GraphBAN drug screening for multiple biomarkers",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input parameters
    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument(
        "--uniprot_ids",
        required=True,
        help="Comma-separated list of UniProt IDs (e.g. P04637,Q9Y6Y9) or path to txt file"
    )
    input_group.add_argument(
        "--zinc_data",
        required=True,
        help="Path to ZINC compounds dataset (CSV with SMILES column)"
    )
    
    # Processing control
    process_group = parser.add_argument_group("Processing Options")
    process_group.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel processes (1=sequential)"
    )
    process_group.add_argument(
        "--chunk_size",
        type=int,
        default=1000,
        help="Number of compounds per prediction chunk"
    )
    
    # (保留原有其他参数...)
    
    return parser.parse_args()

def load_uniprot_ids(input_arg: str) -> List[str]:
    """Load UniProt IDs from either comma string or file"""
    if os.path.exists(input_arg):
        with open(input_arg) as f:
            return [line.strip() for line in f if line.strip()]
    return [uid.strip() for uid in input_arg.split(",")]

def process_single_target(uniprot_id: str, args):
    """Process one UniProt ID with shared arguments"""
    target_dir = os.path.join(args.output_dir, uniprot_id)
    os.makedirs(target_dir, exist_ok=True)
    
    print(f"\nProcessing {uniprot_id}...")
    try:
        # 1. Download FASTA
        fasta_path = download_uniprot_fasta(uniprot_id, target_dir)
        
        # 2. Load compounds in chunks
        for chunk in pd.read_csv(args.zinc_data, chunksize=args.chunk_size):
            # 3. Predict interactions
            results = predict_interactions(fasta_path, chunk)
            
            # 4. Apply filters
            filtered = apply_filters(results)
            
            # 5. Save per-target results
            save_results(filtered, target_dir, f"{uniprot_id}_results.csv")
            
        return (uniprot_id, True, None)
    except Exception as e:
        return (uniprot_id, False, str(e))

def main():
    args = parse_args()
    validate_args(args)
    
    # Load UniProt IDs
    uniprot_ids = load_uniprot_ids(args.uniprot_ids)
    print(f"Loaded {len(uniprot_ids)} targets to process")
    
    # Create processing function with shared args
    processor = partial(process_single_target, args=args)
    
    # Parallel/sequential execution
    if args.parallel > 1:
        with Pool(args.parallel) as pool:
            results = pool.map(processor, uniprot_ids)
    else:
        results = [processor(uid) for uid in uniprot_ids]
    
    # Generate consolidated report
    generate_consolidated_report(results, args.output_dir)
    
    # Print summary
    success = sum(1 for r in results if r[1])
    print(f"\nProcessing complete. Success: {success}/{len(uniprot_ids)}")

def generate_consolidated_report(results: List[Tuple], output_dir: str):
    """Combine results from multiple targets"""
    all_data = []
    for uniprot_id, success, error in results:
        if success:
            target_file = os.path.join(output_dir, uniprot_id, f"{uniprot_id}_results.csv")
            if os.path.exists(target_file):
                df = pd.read_csv(target_file)
                df['Target'] = uniprot_id
                all_data.append(df)
    
    if all_data:
        consolidated = pd.concat(all_data)
        save_path = os.path.join(output_dir, "consolidated_results.csv")
        consolidated.to_csv(save_path, index=False)
        print(f"Consolidated results saved to {save_path}")

# (保留原有其他函数...)

if __name__ == "__main__":
    main()
