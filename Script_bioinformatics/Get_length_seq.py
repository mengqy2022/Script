import argparse

def filter_fasta(input_file, output_file, min_length):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        sequence_id = None
        sequence = ''
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):  # This is a sequence header
                if sequence_id:  # Check if we were processing a sequence
                    if len(sequence) >= min_length:  # Check the length of the last sequence
                        outfile.write(f"{sequence_id}\n{sequence}\n")
                
                sequence_id = line  # Store the current sequence ID
                sequence = ''  # Reset the sequence holder
            else:
                sequence += line  # Append to the sequence
            
        # Check the last sequence
        if sequence_id and len(sequence) >= min_length:
            outfile.write(f"{sequence_id}\n{sequence}\n")

def main():
    parser = argparse.ArgumentParser(description='Filter sequences in a FASTA file based on length.')
    parser.add_argument('input_file', type=str, help='Path to the input FASTA file.')
    parser.add_argument('output_file', type=str, help='Path to the output FASTA file.')
    parser.add_argument('--length', type=int, default=200, help='Minimum length of sequences to keep (default: 200 bp).')

    args = parser.parse_args()

    filter_fasta(args.input_file, args.output_file, args.length)

if __name__ == "__main__":
    main()
