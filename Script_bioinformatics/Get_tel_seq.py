#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/10/30

"""
        >>> 该脚本用于从FASTA文件中提取端粒和非端粒序列 <<<

        输入为一个FASTA格式的文件，输出将生成以下文件：
        1. 包含端粒序列的FASTA文件
        2. 包含非端粒序列的FASTA文件
        3. 统计信息文件
        4. 多余碱基分布统计文件
        5. 一端有端粒序列的FASTA文件（可选）
        6. 两端有端粒序列的FASTA文件（可选）


        使用方法：
        python Get_tel_seq.py -f <输入FASTA文件路径> 
                             [-ot <输出端粒序列文件路径>] 
                             [-ont <输出非端粒序列文件路径>] 
                             [-oe <输出一端有端粒序列的FASTA文件路径>] 
                             [-ob <输出两端有端粒序列的FASTA文件路径>] 
                             [-os <输出统计信息文件路径>] 

        默认情况下，输出文件将分别命名为：
        - 端粒序列文件：'telomeres.fasta'
        - 非端粒序列文件：'non_telomeres.fasta'
        - 统计信息文件：'stats.txt'
"""

import argparse
import re
from datetime import datetime

def read_fasta(file_path):
    sequences = {}
    try:
        with open(file_path, 'r') as fasta_file:
            current_tag = None
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    current_tag = line[1:]  # 去除前导的>
                    sequences[current_tag] = ""
                elif current_tag:
                    sequences[current_tag] += line
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
    one_end_telomere_sequences = {}
    both_ends_telomere_sequences = {}

    telomere_start_pattern = re.compile(r'[CA]*CCCCAAAACCCC', re.IGNORECASE)
    telomere_end_pattern = re.compile(r'GGGGTTTTGGGG[GT]*', re.IGNORECASE)

    for tag, sequence in sequences.items():
        sequence_upper = sequence.upper()
        has_start = telomere_start_pattern.search(sequence_upper) is not None
        has_end = telomere_end_pattern.search(sequence_upper) is not None

        if has_start or has_end:
            telomere_sequences[tag] = sequence_upper
            if has_start and has_end:
                both_ends_telomere_sequences[tag] = sequence_upper
            else:
                one_end_telomere_sequences[tag] = sequence_upper  # 只在一端有端粒
        else:
            non_telomere_sequences[tag] = sequence_upper

    return telomere_sequences, non_telomere_sequences, one_end_telomere_sequences, both_ends_telomere_sequences

def write_statistics(non_telomere_sequences, telomere_sequences, one_end_telomere_sequences, both_ends_telomere_sequences, output_stat_file, original_sequence_count, original_total_length):
    try:
        with open(output_stat_file, 'w') as stat_f:
            stat_f.write("类型\t数量\t总长度\n")
            stat_f.write(f"原始序列\t{original_sequence_count}\t{original_total_length}\n") 
            stat_f.write(f"非端粒序列\t{len(non_telomere_sequences)}\t{sum(map(len, non_telomere_sequences.values()))}\n")
            stat_f.write(f"端粒序列\t{len(telomere_sequences)}\t{sum(map(len, telomere_sequences.values()))}\n")
            stat_f.write(f"一端端粒序列\t{len(one_end_telomere_sequences)}\t{sum(map(len, one_end_telomere_sequences.values()))}\n")
            stat_f.write(f"两端端粒序列\t{len(both_ends_telomere_sequences)}\t{sum(map(len, both_ends_telomere_sequences.values()))}\n")
    except Exception as e:
        print(f"写入统计信息时发生错误: {e}")
        raise

def write_output(sequences, file_path):
    try:
        with open(file_path, 'w') as output_file:
            for tag, seq in sequences.items():
                output_file.write(f">{tag}\n{seq}\n")
    except Exception as e:
        print(f"写入输出文件时发生错误: {e}")
        raise

def process_fasta(fasta_file, output_telomere_file='telomere_output.fasta', 
                  output_non_telomere_file='non_telomere_output.fasta',
                  output_one_end_file=None,
                  output_both_ends_file=None,
                  output_stat_file='stats.txt'):
    
    sequences = read_fasta(fasta_file)
    telomere_sequences, non_telomere_sequences, one_end_telomere_sequences, both_ends_telomere_sequences = extract_telomere_sequences(sequences)

    # 计算原始序列的数量和总长度
    original_sequence_count = len(sequences)
    original_total_length = sum(len(seq) for seq in sequences.values())

    write_output(telomere_sequences, output_telomere_file) if telomere_sequences else None
    write_output(non_telomere_sequences, output_non_telomere_file) if non_telomere_sequences else None
    write_output(one_end_telomere_sequences, output_one_end_file) if output_one_end_file and one_end_telomere_sequences else None
    write_output(both_ends_telomere_sequences, output_both_ends_file) if output_both_ends_file and both_ends_telomere_sequences else None

    # 写入统计信息
    write_statistics(non_telomere_sequences, telomere_sequences, one_end_telomere_sequences, both_ends_telomere_sequences, output_stat_file, original_sequence_count, original_total_length)

def main():
    print(__doc__)

    parser = argparse.ArgumentParser(description='从FASTA文件中筛选具有端粒和不具有端粒序列。')
    parser.add_argument('-f', '--fasta_file', type=str, help='输入FASTA文件路径', required=True)
    parser.add_argument('-ot', '--output_telomere_file', type=str, help='输出端粒序列的FASTA文件路径', default='telomere_output.fasta')
    parser.add_argument('-ont', '--output_non_telomere_file', type=str, help='输出非端粒序列的FASTA文件路径', default='non_telomere_output.fasta')
    parser.add_argument('-oe', '--output_one_end_file', type=str, help='输出一端有端粒序列的FASTA文件路径')
    parser.add_argument('-ob', '--output_both_ends_file', type=str, help='输出两端有端粒序列的FASTA文件路径')
    parser.add_argument('-os', '--output_stat_file', type=str, help='输出统计信息的文件路径', default='stats.txt')
    
    args = parser.parse_args()
    
    process_fasta(args.fasta_file, args.output_telomere_file,
                  args.output_non_telomere_file, args.output_one_end_file,
                  args.output_both_ends_file, args.output_stat_file)

    now = datetime.now()
    print('Done!\n')
    print('结束时间: ' + now.strftime('%Y-%m-%d %H:%M:%S') + '\n')
    print(f"端粒序列文件路径: {args.output_telomere_file}")
    print(f"非端粒序列文件路径: {args.output_non_telomere_file}")
    print(f"统计信息文件路径: {args.output_stat_file}\n")
    if args.output_one_end_file:
        print(f"一端有端粒序列文件路径: {args.output_one_end_file}")
    if args.output_both_ends_file:
        print(f"两端有端粒序列文件路径: {args.output_both_ends_file}\n")
    print('\t>>> mengqingyao <<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()
