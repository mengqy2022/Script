#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/03/08
# @Description : 将多条序列合并成一条序列，并输出到文件

from Bio import SeqIO
import sys

def usage():
    print('\nDescription: Combine multiple sequences into a single sequence.\n')
    print('Usage: python3 merged_fasta_one.py [INPUT FILE] [OUTPUT FILE]\n')
    print('Example: python3 merged_fasta_one.py input.fasta output.fasta\n')
    print('    >>>> mqy <<<<\n')
    sys.exit()

def merge_fasta(fasta):
    global merged_seq
    merged_seq = ""
    records = SeqIO.parse(fasta, "fasta")
    for record in records:
        merged_seq += record.seq

def output_fasta(output):
    with open(output, "w") as output_handle:
        SeqIO.write(SeqIO.SeqRecord(merged_seq, id="merged_one", description=""), output_handle, "fasta")

def main():
    if len(sys.argv)!= 3:
        usage()
        sys.exit()
    fasta = sys.argv[1]
    output = sys.argv[2]
    merge_fasta(fasta)
    output_fasta(output)

try:
    main()
    print('Done!\n')
    print('Output file: '+sys.argv[2]+'\n')
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')
except KeyboardInterrupt:
    usage()