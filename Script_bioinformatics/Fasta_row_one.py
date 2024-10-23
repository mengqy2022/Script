#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/03/08
# @Description : 将序列变成一行，并输出到文件中

import sys

def usage():
    print('\nDescription: Turn the sequence into a single line and output it to a file.\n')
    print('Usage: python3 fasta_row_one.py [INPUT FILE] > [OUTPUT FILE]\n')
    print('Example: python3 fasta_row_one.py input.fasta > output.fasta\n')
    print('    >>>> mqy <<<<\n')
    sys.exit()


dict_seq = {}

def read_fasta(fasta):
    for line in open(fasta):
        if line[0] == '>':
            key = line.strip()
            dict_seq[key] = []
        else:
            dict_seq[key].append(line.strip())


def main():
    read_fasta(sys.argv[1])
    for key,value in dict_seq.items():
        print(key)
        print(''.join(value))

try:
    main()
except IndexError:
    usage()








