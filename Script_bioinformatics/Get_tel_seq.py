#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/10/30

import argparse
import os
import re

def read_fasta(file_path):
    sequences = {}
    current_tag = None

    try:
        with open(file_path, 'r') as fasta_file:
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    current_tag = line[1:]  # 获取标签
                    sequences[current_tag] = ""
                elif current_tag:  # 确保当前标签不为空
                    sequences[current_tag] += line  # 连接序列
    except FileNotFoundError:
        print(f"错误: 找不到文件 {file_path}")
        raise
    except Exception as e:
        print(f"读取文件时发生错误: {e}")
        raise

    return sequences

def extract_telomere_sequences(sequences):
    telomere_sequences = {}
    non_telomere_sequences = {}
    telomere_start_pattern = r'[CA]*CCCCAAAACCCC'
    telomere_end_pattern = r'GGGGTTTTGGGG[GT]*'

    for tag, sequence in sequences.items():
        sequence_upper = sequence.upper()  # 先转换为大写，避免重复调用
        if re.search(telomere_start_pattern, sequence_upper, re.IGNORECASE) or \
           re.search(telomere_end_pattern, sequence_upper, re.IGNORECASE):
            telomere_sequences[tag] = sequence_upper
        else:
            non_telomere_sequences[tag] = sequence_upper

    return telomere_sequences, non_telomere_sequences

def write_output(sequences, output_file):
    try:
        # 先创建文件，即使没有序列也要存在
        with open(output_file, 'w') as out_f:
            if sequences:  # 如果有序列，写入序列
                for tag, sequence in sequences.items():
                    out_f.write(f">{tag}\n{sequence}\n")
            else:
                out_f.write("")  # 确保文件在没有序列时仍被创建
    except Exception as e:
        print(f"写入文件时发生错误: {e}")
        raise

#  驱动函数
def process_fasta(fasta_file, output_telomere_file='telomere_output.fasta', output_non_telomere_file='non_telomere_output.fasta'):
    sequences = read_fasta(fasta_file)
    telomere_sequences, non_telomere_sequences = extract_telomere_sequences(sequences)

    # 只在输出文件存在时写入
    if telomere_sequences:
        write_output(telomere_sequences, output_telomere_file)

    if non_telomere_sequences:
        write_output(non_telomere_sequences, output_non_telomere_file)

def main():
    parser = argparse.ArgumentParser(description='从FASTA文件中提取端粒和非端粒序列。')
    
    parser.add_argument('-f', '--fasta_file', type=str, help='输入FASTA文件路径', required=True)
    parser.add_argument('-ot', '--output_telomere_file', type=str, help='输出端粒序列的FASTA文件路径', default='telomere_output.fasta')
    parser.add_argument('-ont', '--output_non_telomere_file', type=str, help='输出非端粒序列的FASTA文件路径', default='non_telomere_output.fasta')
    
    print("""
                >>> mqy <<<
    该脚本用于从FASTA文件中提取端粒和非端粒序列。 
    输入为一个FASTA格式的文件，输出为两个文件：
    1. 端粒序列的FASTA文件
    2. 非端粒序列的FASTA文件

    使用方法：
    python Get_tel_seq.py -f <输入FASTA文件路径> [-ot <输出端粒序列文件路径>] [-ont <输出非端粒序列文件路径>]
    默认情况下，输出文件分别为'telomere_output.fasta'和'non_telomere_output.fasta'
    """)

    args = parser.parse_args()
    # 打印脚本说明


    process_fasta(args.fasta_file, args.output_telomere_file, args.output_non_telomere_file)

if __name__ == "__main__":
    main()
