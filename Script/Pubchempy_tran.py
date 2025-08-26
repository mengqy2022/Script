import argparse
from pubchempy import get_cids

def main():
    parser = argparse.ArgumentParser(description='Convert SMILES to PubChem CIDs')
    parser.add_argument('-i', '--input', dest='input_file', 
                       required=True, help='Input file containing SMILES (with optional header)')
    parser.add_argument('-o', '--output', dest='output_file',
                       required=True, help='Output file for SMILES and CIDs')
    parser.add_argument('--no-header', action='store_true',
                       help='Input file has no header row')
    args = parser.parse_args()

    # Read SMILES from input file
    with open(args.input_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    # Handle header row
    if not args.no_header:
        # Check if first line looks like a header (contains non-SMILES text)
        if len(lines) > 0 and any(c.isalpha() for c in lines[0].split()[0]):
            header = lines[0]
            smiles_list = lines[1:]
        else:
            smiles_list = lines
    else:
        smiles_list = lines

    # Process SMILES and write output
    with open(args.output_file, 'w') as out_f:
        out_f.write("SMILES\tCID(s)\n")  # Consistent output header
        
        for line in smiles_list:
            # Take first column in case file has multiple columns
            smiles = line.split()[0]  
            try:
                cids = get_cids(smiles, 'smiles')
                cid_list = list(cids) if cids else ['N/A']
                out_f.write(f"{smiles}\t{';'.join(map(str, cid_list))}\n")
            except Exception as e:
                print(f"Error processing SMILES {smiles}: {e}")
                out_f.write(f"{smiles}\tError\n")

if __name__ == '__main__':
    main()
