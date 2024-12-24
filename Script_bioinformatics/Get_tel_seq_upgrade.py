#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqingyao
# @Time    : 2024/12/18
# @version : 1.0
# @File    : Get_tel_seq_upgrade_1.py

import argparse
import re
import datetime
import statistics

def print_colored(text, color):
    color_codes = {
        'purple': '\033[95m',
        'green': '\033[92m',
        'red': '\033[91m',
        'reset': '\033[0m'
    }
    print(f"{color_codes.get(color, '')}{text}{color_codes['reset']}")

current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

print_colored("\n        >>> 该脚本用于从FASTA文件中提取端粒和非端粒序列 <<<\n", 'purple')
print("        输入为一个FASTA格式的文件，输出将生成以下文件：")
print("        1. 包含端粒序列的FASTA文件")
print("        2. 包含非端粒序列的FASTA文件")
print("        3. 统计信息文件")
print("        4. 多余碱基统计摘要文件（可选）\n")
print_colored('            包含一个推荐去除的序列文件：gc_g_q3.txt', 'purple')
print_colored('                (1) 两端都具有端粒的序列不进行额外处理，即使它是其序列中的一小部分；', 'green')
print_colored('                (2) 对单端端粒的序列，我们对其进行多余碱基的计算，并求出多余碱基的均值和Q3值；', 'green')
print_colored('                (3) 若多余碱基的均值大于Q3，则挑选出多余碱基大于q3+1.5*(IQR)的序列IDs；', 'green')
print_colored('                (4) 再统计所有单端序列和多余碱基大于均值序列的GC含量，并得到单端序列GC含量Q3值;', 'green')
print_colored('                (5) 获得单端序列GC含量Q3值后，挑选出GC含量中大于Q3值的序列; \n', 'green')
print('        5. 多余碱基统计文件（可选）')
print('        6. 一端有端粒序列的FASTA文件（可选）')
print('        7. 两端有端粒序列的FASTA文件（可选）')
print_colored('\n        默认情况下，输出文件将分别命名为：\n', 'purple')
print_colored('        - 端粒序列文件：\'telomeres.fasta\'', 'green')
print_colored('        - 非端粒序列文件：\'non_telomeres.fasta\'', 'green')
print_colored('        - 统计信息文件：\'reads_stats.txt\'\n', 'green')
print_colored(f"          当前日期: {current_date}\n", 'green')

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

def calculate_gc_content(sequence):
    """计算给定DNA序列的GC含量和Q3值"""
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    gc_content = (g_count + c_count) / len(sequence) if len(sequence) > 0 else 0
    return gc_content

def calculate_gc_q3(sequences):
    """计算给定序列的GC含量并返回GC Q3值"""
    gc_contents = [calculate_gc_content(seq) for seq in sequences.values()]
    
    if gc_contents:
        q3_value = statistics.quantiles(gc_contents, n=4)[2]  # 获取Q3值
    else:
        q3_value = 0
    
    return q3_value  # 确保只返回Q3值

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
        has_start = telomere_start_pattern.search(sequence_upper)
        has_end = telomere_end_pattern.search(sequence_upper)

        if has_start or has_end:
            telomere_sequences[tag] = sequence_upper
            if has_start and has_end:
                both_ends_telomere_sequences[tag] = sequence_upper
            else:
                one_end_telomere_sequences[tag] = sequence_upper
                upstream_extra = has_start.start() if has_start else 0
                downstream_extra = len(sequence_upper) - (has_end.end() if has_end else len(sequence_upper))
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
            stat_f.write("类型\t数量\t总长度(bd)\t总长度(Mb)\n")
            stat_f.write(f"原始序列\t{original_sequence_count}\t{original_total_length}\t{original_total_length/1000000:.2f}\n") 
            stat_f.write(f"非端粒序列\t{len(non_telomere_sequences)}\t{sum(map(len, non_telomere_sequences.values()))}\t{sum(map(len, non_telomere_sequences.values()))/1000000:.2f}\n")
            stat_f.write(f"端粒序列\t{len(telomere_sequences)}\t{sum(map(len, telomere_sequences.values()))}\t{sum(map(len, telomere_sequences.values()))/1000000:.2f}\n")
            stat_f.write(f"一端端粒序列\t{len(one_end_telomere_sequences)}\t{sum(map(len, one_end_telomere_sequences.values()))}\t{sum(map(len, one_end_telomere_sequences.values()))/1000000:.2f}\n")
            stat_f.write(f"两端端粒序列\t{len(both_ends_telomere_sequences)}\t{sum(map(len, both_ends_telomere_sequences.values()))}\t{sum(map(len, both_ends_telomere_sequences.values()))/1000000:.2f}\n")
    except Exception as e:
        raise RuntimeError(f"写入统计信息时发生错误: {e}")

def write_excess_counts(one_end_excess_counts, output_excess_file):
    """写入一端端粒序列的多余碱基统计到文件"""
    if output_excess_file:
        try:
            with open(output_excess_file, 'w') as excess_f:
                excess_f.write("IDS\t上端多余碱基\t下端多余碱基\n")
                for tag, (upstream, downstream) in one_end_excess_counts.items():
                    if upstream or downstream:
                        excess_f.write(f"{tag}\t{upstream}\t{downstream}\n")
        except Exception as e:
            raise RuntimeError(f"写入多余碱基统计文件时发生错误: {e}")

def write_excess_summary(one_end_excess_counts, output_summary_file, output_recommended_removal_file='recommended_removal.txt'):
    """计算并写入多余碱基统计摘要"""
    if output_summary_file:
        try:
            upstream_values = [upstream for upstream, _ in one_end_excess_counts.values() if upstream > 0]
            downstream_values = [downstream for _, downstream in one_end_excess_counts.values() if downstream > 0]

            def stats(values):
                """计算统计数据"""
                count = len(values)
                if count == 0:
                    return {'count': 0, 'max': 0, 'min': 0, 'avg': 0, 'q1': 0, 'q3': 0, 'upper_limit': 0}
                max_val = max(values)
                min_val = min(values)
                avg = int(statistics.mean(values))
                q1 = statistics.quantiles(values, n=4)[0]
                q3 = statistics.quantiles(values, n=4)[2]
                upper_limit = q3 + 1.5 * (q3 - q1)
                return {'count': count, 'max': max_val, 'min': min_val, 'avg': avg, 'q1': q1, 'q3': q3, 'upper_limit': upper_limit}

            upstream_stats = stats(upstream_values)
            downstream_stats = stats(downstream_values)

            with open(output_summary_file, 'w') as summary_f:
                summary_f.write("类型\t数量\t最大值\t最小值\t平均值\tQ1\tQ3\n")
                summary_f.write(f"上端多余碱基\t{upstream_stats['count']}\t{upstream_stats['max']}\t{upstream_stats['min']}\t{upstream_stats['avg']}\t{upstream_stats['q1']}\t{upstream_stats['q3']}\n")
                summary_f.write(f"下端多余碱基\t{downstream_stats['count']}\t{downstream_stats['max']}\t{downstream_stats['min']}\t{downstream_stats['avg']}\t{downstream_stats['q1']}\t{downstream_stats['q3']}\n")

            if output_recommended_removal_file:
                with open(output_recommended_removal_file, 'w') as id_f:
                    # 检查平均值是否大于Q3
                    if upstream_stats['avg'] > upstream_stats['q3'] and downstream_stats['avg'] > downstream_stats['q3']:
                        # 上游多余碱基大于平均值的IDs
                        upstream_ids = [tag for tag, (upstream, _) in one_end_excess_counts.items() if upstream > upstream_stats['upper_limit']]
                        id_f.write("\n".join(upstream_ids) + "\n")
                        # 下游多余碱基大于平均值的IDs
                        downstream_ids = [tag for tag, (_, downstream) in one_end_excess_counts.items() if downstream > downstream_stats['upper_limit']]
                        id_f.write("\n".join(downstream_ids) + "\n")

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

def write_recommended_gc_sequences(one_end_gc_q3, sequences, recommended_removal_file, output_gc_file):
    """写入推荐去除的GC序列到文件"""
    try:
        with open(recommended_removal_file, 'r') as rec_file:
            recommended_ids = {line.strip() for line in rec_file if line.strip()}  # 读取推荐去除的ID

        with open(output_gc_file, 'w') as gc_output:
            for tag in recommended_ids:
                if tag in sequences:
                    gc_content = calculate_gc_content(sequences[tag])
                    if gc_content > one_end_gc_q3:
                        gc_output.write(f">{tag}\n{sequences[tag]}\n")  # 输出大于Q3值的序列
    except Exception as e:
        raise RuntimeError(f"写入推荐去除的GC序列时发生错误: {e}")

def process_fasta(fasta_file, output_telomere_file='telomere_output.fasta', 
                  output_non_telomere_file='non_telomere_output.fasta',
                  output_one_end_file=None,
                  output_both_ends_file=None,
                  output_stat_file='reads_stats.txt',
                  output_excess_file=None,  
                  output_excess_summary_file=None,
                  output_recommended_removal_file='base_avg_g_q3.txt',
                  output_gc_sequences_file='gc_g_q3.txt'):  # 设置输出文件路径

    sequences = read_fasta(fasta_file)
    telomere_sequences, non_telomere_sequences, one_end_telomere_sequences, both_ends_telomere_sequences, one_end_excess_counts = extract_telomere_sequences(sequences)

    original_sequence_count = len(sequences)
    original_total_length = sum(len(seq) for seq in sequences.values())

    write_output(telomere_sequences, output_telomere_file) if telomere_sequences else None
    write_output(non_telomere_sequences, output_non_telomere_file) if non_telomere_sequences else None
    if output_one_end_file:
        write_output(one_end_telomere_sequences, output_one_end_file)
    if output_both_ends_file:
        write_output(both_ends_telomere_sequences, output_both_ends_file)

    write_statistics(non_telomere_sequences, telomere_sequences, one_end_telomere_sequences, 
                     both_ends_telomere_sequences, one_end_excess_counts, 
                     output_stat_file, original_sequence_count, original_total_length)

    write_excess_counts(one_end_excess_counts, output_excess_file)
    write_excess_summary(one_end_excess_counts, output_excess_summary_file, output_recommended_removal_file)

    # 计算单端序列的GC Q3值
    if one_end_telomere_sequences:
        one_end_gc_q3 = calculate_gc_q3(one_end_telomere_sequences)
        print_colored(f" >>> 单端端粒序列的GC Q3值为: {one_end_gc_q3:.2f} <<<\n", 'red')
        write_recommended_gc_sequences(one_end_gc_q3, sequences, output_recommended_removal_file, output_gc_sequences_file)

def main():
    """主函数，处理命令行参数和程序入口"""
    #print(script_description)

    parser = argparse.ArgumentParser(description='从FASTA文件中筛选具有端粒和不具有端粒序列。',
                                     epilog=print_colored('\tVersion : 1.0\n\t更详细的信息请访问: https://mengqy2022.github.io/genomics/telomere/\n','green'))
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

    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    print('结束时间: ' + current_date + '\n')
    print(f"端粒序列文件路径: {args.output_telomere_file}")
    print(f"非端粒序列文件路径: {args.output_non_telomere_file}")
    print(f"统计信息文件路径: {args.output_stat_file}\n")
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