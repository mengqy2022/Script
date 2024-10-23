#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def descriptive():
    print('descriptive: Batch extraction of sequences [@]')

def usage():
    print('Usage: python3 extract_ids_reads.py [input_file] > [outfile]')

# read fasta
def main():
    re = {}  
    with open (sys.argv[1]) as f:
        for line in f:
            seq = []
            if line.startswith('@'):
                id = line.split(' ')[0].split('_')  #切片分割序列名称
                id = '_'.join(id[:4])  #合并切片的前5部分
            else:
                seq.append(line)
            if id not in re:
                re[id] = seq
            else:
                re[id] += seq
    # target proteinid
    df_site_pr=df_site['proteinid']
    vl=df_site_pr.values.tolist()
    vl=list(set(vl))
    #extract seq_dict
    filterseq = {}
    for i in vl:
        for k,v in re.items():
            if i in k:
                seq = v
                filterseq[k] = seq
    #save fasta
    with open(sys.argv[2],'w') as f:
        for k,v in filterseq.items():
            # f.write( k +'\n')
            # tem = v.__str__()
            # tem1 = tem.replace(',', '')
            # tem2 = tem1.replace('\n', '')
            # f.write(tem2 +'\n')
            f.write("{}\n{}".format(k, ''.join(v)))

try:
    main()
except IndexError:
    descriptive()
    usage()