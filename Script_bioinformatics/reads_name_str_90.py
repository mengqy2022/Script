#!/usr/bin/env python3
# **************************************************
# 贝叶斯进化树需要序列名称在99以下，所以需要对序列名称进行截断
# 可以处理单个文件和多个文件
# Date   : 2024-10-12
# Author : 孟庆瑶
# **************************************************


import argparse
import os

class FastaProcessor:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def process_fasta(self):
        with open(self.input_file, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        processed_lines = []
        for line in lines:
            if line.startswith('>'):  # 判断是否为序列名称行
                if len(line.strip()) > 90:
                    # 如果序列名称长于 90 字符，则进行截断
                    line = line[:90] + '\n'
            processed_lines.append(line)

        self.save_processed_fasta(processed_lines)

    def save_processed_fasta(self, processed_lines):
        output_path = f"{self.output_file}_{os.path.basename(self.input_file)}"
        output_path = output_path.replace('.fasta', '_processed.fasta')
        with open(output_path, 'w', encoding='utf-8') as file:
            file.writelines(processed_lines)
        print(f"\n处理后的 FASTA 文件已保存至: {output_path}")

def main():
    parser = argparse.ArgumentParser(description='处理 FASTA 文件，使序列名称不超过 90 个字符.')
    parser.add_argument('-i','--input_files', type=str, nargs='+', help='输入的 FASTA 文件路径.', required=True)
    parser.add_argument('-o','--output_file', type=str, default='output', help='输出文件的基本名称. 默认值为 "output".', required=True)
    args = parser.parse_args()

    for input_file in args.input_files:
        fasta_processor = FastaProcessor(input_file, args.output_file)
        fasta_processor.process_fasta()

    print('\n运行结束！')
    
if __name__ == "__main__":
    main()
