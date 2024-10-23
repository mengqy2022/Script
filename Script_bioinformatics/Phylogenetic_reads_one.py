#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/09/29
# @File    : Phylogenetic_reads_one

import argparse
import os
from Bio import SeqIO

class FastaMerger:
    def __init__(self, input_folder, output_folder):
        self.input_folder = input_folder
        self.output_folder = output_folder
        
        # 确保输出文件夹存在
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def merge_sequences(self):
        for filename in os.listdir(self.input_folder):
            if filename.endswith('.fa') or filename.endswith('.fasta') or filename.endswith('.faa'):
                file_path = os.path.join(self.input_folder, filename)
                merged_sequence = self.process_fasta(file_path)

                # 生成输出文件路径，与输入文件同名
                output_file_path = os.path.join(self.output_folder, filename)
                
                # 保存合并后的序列到文件
                with open(output_file_path, 'w') as output_file:
                    output_file.write(merged_sequence)
                print(f"文件 '{output_file_path}' 已保存。")

    def process_fasta(self, file_path):
        merged_sequence = ''

        # 读取FASTA文件并合并序列
        for record in SeqIO.parse(file_path, "fasta"):
            merged_sequence += str(record.seq)  # 提取序列字符串并直接连接
        
        return f">{os.path.basename(file_path).split('.')[0]}\n{merged_sequence}\n"  # 返回新的FASTA格式字符串

if __name__ == "__main__":
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description='Merge sequences within each FASTA file into a single sequence.')
    
    # 定义命令行参数
    parser.add_argument('-i', '--input_folder', type=str, help='Path to the input folder containing FASTA files.', required=True)
    parser.add_argument('-o', '--output_folder', type=str, help='Path to the output folder for the merged sequences.', required=True)

    # 解析命令行参数
    args = parser.parse_args()

    # 创建FastaMerger实例并执行合并操作
    merger = FastaMerger(args.input_folder, args.output_folder)
    merger.merge_sequences()

    print("任务完成，合并后的序列已保存到：" + args.output_folder)
