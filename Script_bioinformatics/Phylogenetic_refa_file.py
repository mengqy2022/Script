#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/11/27
# @File    : extract_phylogenetic_data.py

import argparse
from Bio import SeqIO

import re

def extract_sequences(input_file, separator):
    with open(input_file, "r") as file:
        for line in file:
            # 检查是否是以'>'开头的行（FASTA格式中的ID行）
            if line.startswith('>'):
                # 去除开头的'>'及换行符
                record_id = line[1:].strip()

                # 检查是否包含自定义的分隔符
                if separator in record_id:
                    # 分割ID和描述信息
                    id_and_desc = record_id.split(separator)

                    new_id = id_and_desc[0]
                    # 提取描述信息
                    description = id_and_desc[1].split(' ')

                    description_new = description[:-1]
                    description_new = '_'.join(description_new)
                    # 使用正则表达式将所有符号（除了'>的部分）替换为'_'
                    new_id = re.sub(r'[^a-zA-Z0-9]', '_', new_id)
                    description_new = re.sub(r'[^a-zA-Z0-9]', '_', description_new)
                    # 确保字符串格式没有多余空格，并保持'>'
                    print(">" + new_id + "_" + description_new)
                else:
                    # 如果没有分隔符，直接处理整行
                    new_id = record_id
                    description_new = record_id.split(' ')
                    description_new = '_'.join(description_new)
                    # 使用正则表达式将所有符号（除了'>的部分）替换为'_'
                    new_id = re.sub(r'[^a-zA-Z0-9]', '_', new_id)
                    description_new = re.sub(r'[^a-zA-Z0-9]', '_', description_new)
                    # 确保字符串格式没有多余空格，并保持'>'
                    print(">" + description_new)

            else:
                # 处理序列行
                sequence = line.strip()
                print(sequence)



def main():
    parser = argparse.ArgumentParser(description='Extract specific sequence IDs and descriptions from a FASTA file.')
    parser.add_argument('input_file', type=str, help='Input FASTA file')
    parser.add_argument('--separator', type=str, default='|', help='Custom separator for ID and description')

    args = parser.parse_args()
    
    extract_sequences(args.input_file, args.separator)

if __name__ == "__main__":
    main()
