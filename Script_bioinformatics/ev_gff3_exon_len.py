#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

def usage():
    print('Description: 从EVM整合的.gff3文件中得到外显子的数量和长度')
    print('Usage: python3 ev_gff3_exon_len.py [EVM.gff3] > [outfile.txt]')

def main():
   with open(sys.argv[1], 'r') as f:
     for line in f:
         line = line.rstrip("\n")  # 删除行尾的换行符
         array = line.split("\t")
         lenght= int(array[4]) - int(array[3] + 1)
         print(lenght)

try:
    main()
except IndexError:
    usage()