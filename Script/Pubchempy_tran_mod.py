#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PubChem CID Fetcher
Version: 1.0
Date: 2025-07-22
Author: mqy
Description: Fetch PubChem CIDs from SMILES or compound names.
"""

import argparse
from pubchempy import get_cids

def main():
    parser = argparse.ArgumentParser(
        description='PubChem CID Fetcher (v1.0) - Convert SMILES or compound names to PubChem CIDs',
        epilog='Example: python fetch_cid.py -i input.txt -o output.txt --type smiles'
    )
    parser.add_argument('-i', '--input', dest='input_file', 
                       required=True, help='Input file containing SMILES or compound names')
    parser.add_argument('-o', '--output', dest='output_file',
                       required=True, help='Output file for CID results')
    parser.add_argument('--no-header', action='store_true',
                       help='Input file has no header row')
    parser.add_argument('--type', dest='input_type', choices=['smiles', 'name'],
                       required=True, help='Input type: "smiles" or "name"')
    args = parser.parse_args()

    # Read input file
    with open(args.input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Handle header row
    if not args.no_header:
        if len(lines) > 0 and any(c.isalpha() for c in lines[0].split()[0]):
            header = lines[0]
            input_list = lines[1:]
        else:
            input_list = lines
    else:
        input_list = lines

    # Process input and write output
    with open(args.output_file, 'w') as out_f:
        # Write appropriate header based on input type
        if args.input_type == 'smiles':
            out_f.write("SMILES\tCID(s)\n")
        else:
            out_f.write("Name\tCID(s)\n")
        
        for line in input_list:
            query = line.split()[0]  # Take first column (in case of multi-column input)
            try:
                cids = get_cids(query, args.input_type)  # 'smiles' or 'name'
                cid_list = list(cids) if cids else ['N/A']
                out_f.write(f"{query}\t{';'.join(map(str, cid_list))}\n")
            except Exception as e:
                print(f"Error processing {args.input_type} '{query}': {e}")
                out_f.write(f"{query}\tError\n")

if __name__ == '__main__':
    main()
