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
        4. 多余碱基统计摘要文件（可选）
        5. 多余碱基统计文件（可选）
        6. 一端有端粒序列的FASTA文件（可选）
        7. 两端有端粒序列的FASTA文件（可选）

        使用方法：
        python Get_tel_seq.py -f <输入FASTA文件路径> 
                             [-ot <输出端粒序列文件路径>] 
                             [-ont <输出非端粒序列文件路径>] 
                             [-oe <输出一端有端粒序列的FASTA文件路径>] 
                             [-ob <输出两端有端粒序列的FASTA文件路径>] 
                             [-os <输出统计信息文件路径>] 
                             [-odr <输出多余碱基统计文件路径>] 
                             [-ods <输出多余碱基统计摘要文件路径>]

        默认情况下，输出文件将分别命名为：
        - 端粒序列文件：'telomeres.fasta'
        - 非端粒序列文件：'non_telomeres.fasta'
        - 统计信息文件：'reads_stats.txt'
"""

import argparse
import re
from datetime import datetime
import statistics

def read_fasta(file_path):
    """读取FASTA文件并返回序列字典"""
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
        raise FileNotFoundError(f"错误: 找不到文件 {file_path}")
    except Exception as e:
        raise RuntimeError(f"读取文件时发生错误: {e}")

    return sequences

def extract_telomere_sequences(sequences):
    """提取端粒序列及其多余碱基统计"""
    telomere_sequences = {}
    non_telomere_sequences = {}
    one_end_telomere_sequences = {}
    both_ends_telomere_sequences = {}
    one_end_excess_counts = {}

    telomere_start_pattern = re.compile(r'[CA]*CCCCAAAACCCC', re.IGNORECASE)
    telomere_end_pattern = re.compile(r'GGGGTTTTGGGG[GT]*', re.IGNORECASE)

    for tag, sequence in sequences.items():
        sequence_upper = sequence.upper()
        has_start = telomere_start_pattern.search(sequence_upper) is not None
        has_end = telomere_end_pattern.search(sequence_upper) is not None

        total_length = len(sequence_upper)

        if has_start or has_end:
            telomere_sequences[tag] = sequence_upper
            if has_start and has_end:
                both_ends_telomere_sequences[tag] = sequence_upper
            else:
                one_end_telomere_sequences[tag] = sequence_upper
                if has_start:
                    upstream_extra = telomere_start_pattern.search(sequence_upper).start()
                else:
                    upstream_extra = 0

                end_matches = list(re.finditer(telomere_end_pattern, sequence_upper))

                if end_matches:
                    downstream_extra = total_length - end_matches[-1].end()
                else:
                    downstream_extra = 0

                one_end_excess_counts[tag] = (upstream_extra, downstream_extra)
        else:
            non_telomere_sequences[tag] = sequence_upper

    return (telomere_sequences, non_telomere_sequences, 
            one_end_telomere_sequences, both_ends_telomere_sequences, 
            one_end_excess_counts)

def write_statistics(non_telomere_sequences, telomere_sequences, one_end_telomere_sequences, 
                     both_ends_telomere_sequences, one_end_excess_counts, 
                     output_stat_file, original_sequence_count, original_total_length):
    """写入统计信息到文件"""
    try:
        with open(output_stat_file, 'w') as stat_f:
            stat_f.write("类型\t数量\t总长度\n")
            stat_f.write(f"原始序列\t{original_sequence_count}\t{original_total_length}\n") 
            stat_f.write(f"非端粒序列\t{len(non_telomere_sequences)}\t{sum(map(len, non_telomere_sequences.values()))}\n")
            stat_f.write(f"端粒序列\t{len(telomere_sequences)}\t{sum(map(len, telomere_sequences.values()))}\n")
            stat_f.write(f"一端端粒序列\t{len(one_end_telomere_sequences)}\t{sum(map(len, one_end_telomere_sequences.values()))}\n")
            stat_f.write(f"两端端粒序列\t{len(both_ends_telomere_sequences)}\t{sum(map(len, both_ends_telomere_sequences.values()))}\n")

    except Exception as e:
        raise RuntimeError(f"写入统计信息时发生错误: {e}")

def write_excess_counts(one_end_excess_counts, output_excess_file):
    """写入一端端粒序列的多余碱基统计到文件"""
    if output_excess_file is not None:  # 只有在指定了输出文件的情况下才写入
        try:
            with open(output_excess_file, 'w') as excess_f:
                excess_f.write("IDS\t上端多余碱基\t下端多余碱基\n")
                for tag, (upstream, downstream) in one_end_excess_counts.items():
                    if upstream != 0 or downstream != 0:  # 只写入上端或下端有多余碱基的情况
                        excess_f.write(f"{tag}\t{upstream}\t{downstream}\n")
        except Exception as e:
            raise RuntimeError(f"写入多余碱基统计文件时发生错误: {e}")

def write_excess_summary(one_end_excess_counts, output_summary_file):
    """计算并写入多余碱基统计摘要"""
    if output_summary_file is not None:  # 只有在指定了输出文件的情况下才写入
        try:
            upstream_values = [upstream for upstream, _ in one_end_excess_counts.values() if upstream > 0]
            downstream_values = [downstream for _, downstream in one_end_excess_counts.values() if downstream > 0]

            upstream_count = len(upstream_values)
            downstream_count = len(downstream_values)

            upstream_max = max(upstream_values) if upstream_values else 0
            downstream_max = max(downstream_values) if downstream_values else 0

            upstream_min = min(upstream_values) if upstream_values else 0
            downstream_min = min(downstream_values) if downstream_values else 0

            upstream_avg = int(statistics.mean(upstream_values)) if upstream_values else 0
            downstream_avg = int(statistics.mean(downstream_values)) if downstream_values else 0

            with open(output_summary_file, 'w') as summary_f:
                summary_f.write("类型\t数量\t最大值\t最小值\t平均值\n")
                summary_f.write(f"上端多余碱基\t{upstream_count}\t{upstream_max}\t{upstream_min}\t{upstream_avg}\n")
                summary_f.write(f"下端多余碱基\t{downstream_count}\t{downstream_max}\t{downstream_min}\t{downstream_avg}\n")
        except Exception as e:
            raise RuntimeError(f"写入多余碱基统计摘要文件时发生错误: {e}")

def write_output(sequences, file_path):
    """写入输出序列到FASTA文件"""
    try:
        with open(file_path, 'w') as output_file:
            for tag, seq in sequences.items():
                output_file.write(f">{tag}\n{seq}\n")
    except Exception as e:
        raise RuntimeError(f"写入输出文件时发生错误: {e}")

def process_fasta(fasta_file, output_telomere_file='telomere_output.fasta', 
                  output_non_telomere_file='non_telomere_output.fasta',
                  output_one_end_file=None,
                  output_both_ends_file=None,
                  output_stat_file='reads_stats.txt',
                  output_excess_file=None,  
                  output_excess_summary_file=None): 
    
    sequences = read_fasta(fasta_file)
    telomere_sequences, non_telomere_sequences, one_end_telomere_sequences, both_ends_telomere_sequences, one_end_excess_counts = extract_telomere_sequences(sequences)

    original_sequence_count = len(sequences)
    original_total_length = sum(len(seq) for seq in sequences.values())

    write_output(telomere_sequences, output_telomere_file) if telomere_sequences else None
    write_output(non_telomere_sequences, output_non_telomere_file) if non_telomere_sequences else None
    write_output(one_end_telomere_sequences, output_one_end_file) if output_one_end_file and one_end_telomere_sequences else None
    write_output(both_ends_telomere_sequences, output_both_ends_file) if output_both_ends_file and both_ends_telomere_sequences else None

    write_statistics(non_telomere_sequences, telomere_sequences, one_end_telomere_sequences, 
                     both_ends_telomere_sequences, one_end_excess_counts, 
                     output_stat_file, original_sequence_count, original_total_length)

    # 在此检查是否有文件路径被指定，才写入多余碱基统计信息和摘要
    write_excess_counts(one_end_excess_counts, output_excess_file)
    write_excess_summary(one_end_excess_counts, output_excess_summary_file)

def main():
    """主函数，处理命令行参数和程序入口"""
    print(__doc__)

    parser = argparse.ArgumentParser(description='从FASTA文件中提取端粒和非端粒序列。')
    parser.add_argument('-f', '--fasta_file', type=str, help='输入FASTA文件路径', required=True)
    parser.add_argument('-ot', '--output_telomere_file', type=str, help='输出端粒序列的FASTA文件路径', default='telomere_output.fasta')
    parser.add_argument('-ont', '--output_non_telomere_file', type=str, help='输出非端粒序列的FASTA文件路径', default='non_telomere_output.fasta')
    parser.add_argument('-oe', '--output_one_end_file', type=str, help='输出一端有端粒序列的FASTA文件路径')
    parser.add_argument('-ob', '--output_both_ends_file', type=str, help='输出两端有端粒序列的FASTA文件路径')
    parser.add_argument('-os', '--output_stat_file', type=str, help='输出统计信息文件路径', default='reads_stats.txt')
    parser.add_argument('-odr', '--output_excess_file', type=str, help='输出多余碱基统计信息文件路径')
    parser.add_argument('-ods', '--output_excess_summary_file', type=str, help='输出多余碱基统计摘要文件路径')

    args = parser.parse_args()
    
    process_fasta(args.fasta_file, args.output_telomere_file,
                  args.output_non_telomere_file, args.output_one_end_file,
                  args.output_both_ends_file, args.output_stat_file, 
                  args.output_excess_file, args.output_excess_summary_file)  

    now = datetime.now()
    print('Done!\n')
    print('结束时间: ' + now.strftime('%Y-%m-%d %H:%M:%S') + '\n')
    print(f"端粒序列文件路径: {args.output_telomere_file}")
    print(f"非端粒序列文件路径: {args.output_non_telomere_file}")
    print(f"统计信息文件路径: {args.output_stat_file}" + '\n')
    if args.output_excess_file:
        print(f"多余碱基统计文件路径: {args.output_excess_file}")
    if args.output_excess_summary_file:
        print(f"多余碱基统计摘要文件路径: {args.output_excess_summary_file}")
    if args.output_one_end_file:
        print(f"一端有端粒序列文件路径: {args.output_one_end_file}")
    if args.output_both_ends_file:
        print(f"两端有端粒序列文件路径: {args.output_both_ends_file}")
    print('\n\t>>> mengqingyao <<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()
