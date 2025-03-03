#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/10/08
# @File    : Phylogenetic_singleortho_separa


import os
import argparse
from Bio import SeqIO
from collections import defaultdict

class FastaProcessor:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder
        self.sequences_dict = defaultdict(list)

        # 创建输出文件夹（如果不存在）
        os.makedirs(self.output_folder, exist_ok=True)

    def process_fasta_files(self):
        # 遍历文件夹中的所有FASTA文件
        for filename in os.listdir(self.input_folder):
            if filename.endswith('.fa') or filename.endswith('.fasta') or filename.endswith('.faa'):  # 只处理FASTA文件
                file_path = os.path.join(self.input_folder, filename)
                self._parse_fasta_file(file_path)

        # 输出分组后的序列到文件
        self._write_output_files()

    def _parse_fasta_file(self, file_path):
        # 使用SeqIO读取FASTA文件
        for record in SeqIO.parse(file_path, "fasta"):
            # 假设序列名称格式为 "part1_part2"
            name_parts = record.id.split('_', 3)
            if len(name_parts) > 1:
                first_part = name_parts[0]  # 第一部分
                self.sequences_dict[first_part].append(record)  # 将记录对象添加到对应的列表中

    def _write_output_files(self):
        # 输出分组后的序列到文件
        for first_part, records in self.sequences_dict.items():
            output_file_path = os.path.join(self.output_folder, f"{first_part}.fasta")

            # 将序列写入文件
            with open(output_file_path, "w") as output_file:
                SeqIO.write(records, output_file, "fasta")

            print(f"已生成文件: {output_file_path}")

def main():
    # 使用argparse解析命令行参数
    parser = argparse.ArgumentParser(description="Processing FASTA files and grouping output by name",
                                        epilog="\t更详细的信息请访问:https://mengqy2022.github.io/genomics/phylogenetic/\n")
    parser.add_argument('-i', '--input_folder', type=str, help='Path to the input folder containing FASTA files.', required=True)
    parser.add_argument('-o', '--output_folder', type=str, help='Output File Path.', required=True)

    args = parser.parse_args()

    # 创建FastaProcessor对象并处理文件
    fasta_processor = FastaProcessor(args.input_folder, args.output_folder)
    fasta_processor.process_fasta_files()
    
    print("任务完成，已保存到：" + args.output_folder)

if __name__ == "__main__":
    main()
