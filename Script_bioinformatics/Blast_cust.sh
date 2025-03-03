#!/bin/bash
# Default parameters

num_threads=20
evalue="1e-5"
max_targets=5

# Parse command-line arguments
while getopts "d:q:o:" opt; do
    case $opt in
        d) db_file="$OPTARG" ;;        # Database file
        q) query_file="$OPTARG" ;;     # Query file
        o) output_prefix="$OPTARG" ;;  # Output prefix (optional)
        *) echo "Invalid option"; exit 1 ;;
    esac
done

# Check for required arguments
if [[ -z "$db_file" || -z "$query_file" ]]; then
    echo "Usage: bash blast.sh -d <database_file> -q <query_file> [-o <output_prefix>]"
    exit 1
fi

# Ensure input files exist
if [[ ! -f "$db_file" ]]; then
    echo "Error: Database file '$db_file' does not exist."
    exit 1
fi

if [[ ! -f "$query_file" ]]; then
    echo "Error: Query file '$query_file' does not exist."
    exit 1
fi

# Set output prefix and directory
output_prefix="${output_prefix:-blast_output}"
output_dir=$(dirname "$(realpath "$output_prefix")")
db_output="${output_dir}/db_blast"

# Ensure output directory exists
mkdir -p "$output_dir"

# Function to detect sequence type (nucl or prot)
detect_sequence_type() {
    local file=$1
    # Check for protein-specific residues
    if grep -vi ">" "$file" | grep -qE "[EFILPQZ]"; then
        echo "prot"
    # Check for valid nucleotide residues
    elif grep -vi ">" "$file" | grep -qE "^[ACGTURYKMSWBDHVNacgturykmswbdhvn]+$"; then
        echo "nucl"
    else
        echo "unknown"
    fi
}

# Detect sequence types for query and database
db_type=$(detect_sequence_type "$db_file")
query_type=$(detect_sequence_type "$query_file")
if [[ "$db_type" == "unknown" || "$query_type" == "unknown" ]]; then
    echo "Error: Unable to detect sequence type for database or query."
    exit 1
fi

# If query is a CDS, ensure appropriate BLAST tool is selected
if [[ "$query_type" == "nucl" && "$db_type" == "prot" ]]; then
    echo "Query appears to be a CDS. Switching to blastx to handle translation."
    blast_tool="blastx"
elif [[ "$db_type" == "nucl" && "$query_type" == "nucl" ]]; then
    blast_tool="blastn"
elif [[ "$db_type" == "prot" && "$query_type" == "prot" ]]; then
    blast_tool="blastp"
elif [[ "$db_type" == "prot" && "$query_type" == "nucl" ]]; then
    blast_tool="tblastn"
else
    echo "Error: Could not determine the appropriate BLAST tool."
    exit 1
fi

echo "Detected database type: $db_type"
echo "Detected query type: $query_type"
echo "Using BLAST tool: $blast_tool"

# Check if the BLAST database already exists; if not, create it
if [[ ! -f "${db_output}.nsq" && ! -f "${db_output}.psq" ]]; then
    echo "Creating BLAST database at ${db_output}..."
    makeblastdb -in "$db_file" -dbtype "$db_type" -parse_seqids -out "$db_output" 2> "${output_dir}/makeblastdb.log"
    if [[ $? -ne 0 ]]; then
        echo "Error: makeblastdb failed. Check ${output_dir}/makeblastdb.log for details."
        exit 1
    fi
    echo "Database created successfully."
else
    echo "BLAST database already exists at ${db_output}. Skipping creation."
fi

# Run the selected BLAST tool with specified formats
declare -A formats=( [6]="6" [7]="7" [10]="10" [0]="0" )
for fmt in "${!formats[@]}"; do
    output_file="${output_prefix}.${formats[$fmt]}.out"
    echo "Running $blast_tool with output format $fmt..."
    $blast_tool \
        -query "$query_file" \
        -db "$db_output" \
        -evalue "$evalue" \
        -max_target_seqs "$max_targets" \
        -outfmt "$fmt" \
        -out "$output_file" \
        -num_threads "$num_threads" 2> "${output_file}.log"
    if [[ $? -ne 0 ]]; then
        echo "Error: $blast_tool failed for output format $fmt. Check ${output_file}.log for details."
    else
        echo "Output for format $fmt saved to $output_file."
    fi
done

echo "BLAST search completed with $blast_tool."
