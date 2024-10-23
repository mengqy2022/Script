#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/10/12
# @File    : Fasta_seq_renames.py


import argparse
from Bio import SeqIO
from datetime import datetime

class FastaRenamer:
    def __init__(self, fasta_file, rename_file=None, prefix='gene', output_file=None):
        self.fasta_file = fasta_file
        self.rename_file = rename_file
        self.prefix = prefix
        self.output_file = output_file
        self.fa_dict = {}

    def read_fasta(self):
        try:
            with open(self.fasta_file, 'r') as fa:
                for seq in SeqIO.parse(fa, "fasta"):
                    self.fa_dict[str(seq.id)] = str(seq.seq).strip()
        except FileNotFoundError as e:
            print(f"错误: 文件未找到 - {e}")

    def rename_num(self):
        rename_fa = {}
        total_count = len(self.fa_dict)
        for count, (seq_id, sequence) in enumerate(self.fa_dict.items(), start=1):
            seq_id_new = f"{self.prefix}_{str(count).zfill(len(str(total_count)))}"
            rename_fa[seq_id_new] = sequence
        return rename_fa

    def rename_id(self):
        rename = {}
        try:
            with open(self.rename_file, "r") as f:
                for line in f:
                    old_id, new_id = line.strip().split("\t")
                    rename[old_id] = new_id.strip()
        except FileNotFoundError as e:
            print(f"错误: 文件未找到 - {e}")
            return {}

        rename_dict = {}
        for old_id, new_id in rename.items():
            if old_id in self.fa_dict:
                rename_dict[new_id] = self.fa_dict[old_id]
            else:
                print("原始的fasta文件不存在ID:", old_id)
        return rename_dict

    def execute(self):
        self.read_fasta()
        
        if self.rename_file:
            res = self.rename_id()
        else:
            res = self.rename_num()

        output_data = "\n".join([f">{i}\n{seq}" for i, seq in res.items()])

        if self.output_file:
            with open(self.output_file, "w") as p:
                p.write(output_data + "\n")
        else:
            print(output_data)
        
        print('Done!\n')
        print('Time: ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        print('Output file: ' + (self.output_file if self.output_file else "屏幕输出") + '\n')
        print('\t>>> mqy <<<<\n')
        print('\t如果有任何问题请及时联系\n')
        print('\t邮箱：<15877464851@163.com>\n')

def main():
    parser = argparse.ArgumentParser(description='此脚本用于对fasta文件中的序列进行重命名。')
    parser.add_argument('-a', '--fasta', help='输入待命名的fasta文件', type=str, required=True)
    parser.add_argument('-l', '--list', help='输入重命名文件，第一列为原始序列ID，第二列为新序列ID', type=str, required=False)
    parser.add_argument('-p', '--prefix', help='序列前缀，默认为gene', type=str, required=False, default='gene')
    parser.add_argument('-o', '--out', help='结果输出文件，默认输出到屏幕', type=str, required=False)
    args = parser.parse_args()

    renamer = FastaRenamer(
        fasta_file=args.fasta,
        rename_file=args.list,
        prefix=args.prefix,
        output_file=args.out
    )
    renamer.execute()

if __name__ == '__main__':
    main()
