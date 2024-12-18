#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/09/29
# @File    : Phylogenetic_rename_fa.py


import argparse
import os
from Bio import SeqIO

class FastaRenamer:
    def __init__(self, folder, separator, output_folder):
        self.folder = folder
        self.separator = separator
        self.output_folder = output_folder
        self.names = []

    def extract_names_from_folder(self):
        """从指定文件夹提取文件名"""
        try:
            file_names = os.listdir(self.folder)
            for file_name in file_names:
                if os.path.isfile(os.path.join(self.folder, file_name)):
                    # 第一次分割
                    first_name_part = file_name.split(self.separator)[0]
                    # 如果包含"."，再进行一次分割
                    if '.' in first_name_part:
                        first_name_part = first_name_part.split('.')[0]
                    # 如果包含"-"，再进行一次分割
                    if '-' in first_name_part:
                        first_name_part = first_name_part.split('-')[0]
                    # 如果包含空格，再进行一次分割
                    if ' ' in first_name_part:
                        first_name_part = first_name_part.split(' ')[0]
                    self.names.append(first_name_part)
        except FileNotFoundError:
            print(f"文件夹 '{self.folder}' 不存在。")

    def process_fasta_file(self, file_path, base_name):
        """处理单个 fasta 文件，将序列重命名并保存到输出文件"""
        fa_dict = {}
        
        # 读取 fasta 文件并将内容存入字典
        with open(file_path, 'r') as fa:
            for seq in SeqIO.parse(fa, "fasta"):
                sequence = str(seq.seq).strip()
                seq_id = str(seq.id)
                fa_dict[seq_id] = sequence

        # 保存到新的 fasta 输出文件，序列名称随着数量递增
        # 新创建输出文件用于写入所有的序列
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
            
        output_file = os.path.join(self.output_folder, f"{base_name}.fa")
        with open(output_file, "w") as p:
            for idx, (seq_id, sequence) in enumerate(fa_dict.items(), start=1):
                # 使用基本名称和递增数量进行重命名
                seq_id_new = f"{base_name}_{idx}"  # example: baseName_1
                p.write(">" + seq_id_new + "\n" + sequence + "\n")  # 写入新的序列名称和序列内容
        
        print(f"文件 '{output_file}' 已保存。")

    def run(self):
        """运行整个流程"""
        self.extract_names_from_folder()

        for i, file_name in enumerate(os.listdir(self.folder)):
            if file_name.endswith('.fasta') or file_name.endswith('.fa') or file_name.endswith('.faa'):
                file_path = os.path.join(self.folder, file_name)
                if i < len(self.names):
                    base_name = self.names[i]
                    self.process_fasta_file(file_path, base_name)

def main():
    parser = argparse.ArgumentParser(description='此脚本用于对fasta文件中的序列进行重命名',
                                    epilog="\t更详细的信息请访问:https://mengqy2022.github.io/genomics/phylogenetic/\n")
    parser.add_argument('-f', '--folder', help='指定文件夹以读取fasta文件', type=str, required=True)
    parser.add_argument('-s', '--separator', help='指定分隔符，默认为"_"', type=str, required=False, default='_')
    parser.add_argument('-o', '--output_folder', help='结果输出文件夹', type=str, required=True)
    args = parser.parse_args()

    renamer = FastaRenamer(args.folder, args.separator, args.output_folder)
    renamer.run()
    print("重命名完成，结果已保存到" + args.output_folder)

if __name__ == '__main__':
    main()
