#!/usr/bin/env python3
# **************************************************
# 将序列格式文件转换为以制表符分割的csv文件或excel文件
# Date   : 2024-10-11
# Author : 孟庆瑶
# **************************************************

from Bio import SeqIO
import pandas as pd
import argparse
import sys
import os

class ArgumentParser:
    def __init__(self):
        self.parser = argparse.ArgumentParser(description='将序列格式文件转换为以制表符分割的csv文件或excel文件')
        self.add_arguments()

    def add_arguments(self):
        self.parser.add_argument('-i','--input_file', type=str, help='输入的fasta文件名', required=True)
        self.parser.add_argument('-o','--output_prefix', type=str, help='输出文件前缀', required=True)
        self.parser.add_argument('-f','--format', choices=['csv', 'excel'], default='csv', 
                                 help='输出文件格式（默认是csv）')

    def parse(self):
        return self.parser.parse_args()

class FastaConverter:
    def __init__(self, input_file, output_prefix, output_format):
        self.input_file = input_file
        self.output_prefix = output_prefix
        self.output_format = output_format

    def read_fasta(self):
        fa_dict = {}
        try:
            with open(self.input_file, 'r') as fa:
                for seq in SeqIO.parse(fa, "fasta"):
                    fa_dict[str(seq.id)] = str(seq.seq).strip()
        except FileNotFoundError:
            print(f"错误：未找到文件 {self.input_file}")
            sys.exit(1)
        except Exception as e:
            print(f"读取文件时发生错误: {e}")
            sys.exit(1)
        return fa_dict

    def create_dataframe(self, fast_dict):
        return pd.DataFrame(list(fast_dict.items()), columns=['IDS', 'sequence'])

    def save_file(self, dataset):
        output_file = f"{self.output_prefix}.{self.output_format}"
        try:
            if self.output_format == 'csv':
                dataset.to_csv(output_file, index=False, sep='\t')
            elif self.output_format == 'excel':
                dataset.to_excel(output_file, index=False, engine='openpyxl')
            print('\n运行成功！')
            print(f'输出文件: {output_file}\n')
        except Exception as e:
            print(f"保存文件时发生错误: {e}")
            sys.exit(1)

    def run(self):
        fast_dict = self.read_fasta()
        dataset = self.create_dataframe(fast_dict)
        self.save_file(dataset)

def main():
    arg_parser = ArgumentParser()
    args = arg_parser.parse()

    converter = FastaConverter(args.input_file, args.output_prefix, args.format)
    converter.run()

if __name__ == '__main__':
    main()
