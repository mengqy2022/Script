'''
Author: Zhigang Han
Email: hanzg@zafu.edu.cn
Version: 1.0
Date: 2023.09.26
Description: This is a script for extracting the longest CDS sequence by genome.fa and genome.gff3
Usage: python extract_and_combine_cds.py -g genome.fa -a genomic.gff3 -o longest_cds.fa
'''

import argparse
from Bio import SeqIO
from collections import defaultdict

def parse_gff3(gff3_file):
    gene_cds_info = defaultdict(list) # dictionary with defaut (list)
    with open(gff3_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                columns = line.strip().split("\t")
                feature_type = columns[2]
                chrom = columns[0]
                start = int(columns[3])
                end = int(columns[4])
                strand = columns[6]
                attributes = columns[8]

                if feature_type == "CDS":
                    parent_id = attributes.split("Parent=")[1].split(";")[0] # cds' parent id
                    gene_id = attributes.split("gene=")[1].split(";")[0] # gene id
                    gene_cds_info[gene_id].append({'chrom': chrom, 'start': start, 'end': end, 'strand': strand, 'parent_id': parent_id})

    return gene_cds_info

def extract_and_combine_cds(genome_file, gene_cds_info, output_file):
    genome_data = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    longest_sequences = {}

    for gene_id, cds_segments in gene_cds_info.items():
        longest_seq_for_gene = ""
        longest_parent_id = ""

        parent_cds_segments = defaultdict(list)
        for segment in cds_segments:
            parent_cds_segments[segment['parent_id']].append(segment) # special style with one parent_id including several cds segmentation

        for parent_id, segments in parent_cds_segments.items():
            if segments[0]['strand'] == '+':
                sorted_segments = sorted(segments, key=lambda x: x['start'])
            else:   # Sort for negative strand of DNA
                sorted_segments = sorted(segments, key=lambda x: x['end'], reverse=True)

            combined_seq = ""
            for segment in sorted_segments:
                chrom = segment['chrom']
                start = segment['start']
                end = segment['end']
                strand = segment['strand']
                sequence = genome_data[chrom].seq[start-1:end]
                if strand == '-':
                    sequence = sequence.reverse_complement()
                combined_seq += sequence
               # print(gene_id,parent_id,len(combined_seq))
            if len(combined_seq) > len(longest_seq_for_gene):
                longest_seq_for_gene = combined_seq
                longest_parent_id = parent_id
               # print(gene_id,parent_id,len(combined_seq),"_longest")

        longest_sequences[gene_id] = longest_seq_for_gene

    with open(output_file, 'w') as out_fasta:
        for gene_id, sequence in longest_sequences.items():
            out_fasta.write(f">{gene_id}\n{sequence}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract and combine CDS from genome.')
    parser.add_argument('-g', '--genome', required=True, help='Path to the genome FASTA file.')
    parser.add_argument('-a', '--annotations', required=True, help='Path to the GFF3 annotations file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output FASTA file for combined CDS.')

    args = parser.parse_args()

    gene_cds_info = parse_gff3(args.annotations)
    extract_and_combine_cds(args.genome, gene_cds_info, args.output)