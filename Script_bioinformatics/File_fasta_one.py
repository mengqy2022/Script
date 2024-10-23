#!/usr/bin/env python

from Bio import SeqIO
 
# 输入的FASTA文件路径
fasta_file = 'input.fasta'
 
# 使用SeqIO读取FASTA文件
with open(fasta_file) as handle:
    first_seq = next(SeqIO.parse(handle, "fasta"))
 
# 打印或保存第一条序列
print(first_seq.seq)
 
# 如果需要保存到新的FASTA文件中
# with open("first_seq.fasta", "w") as out_handle:
#     SeqIO.write(first_seq, out_handle, "fasta")