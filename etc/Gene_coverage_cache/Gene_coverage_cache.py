#!/usr/bin/env python
import os
import argparse
import concurrent.futures
from Bio import SeqIO

# 定义支持的基因组和CDS文件后缀
GENOME_EXTENSIONS = ['.fa', '.fasta', '.fna']
CDS_EXTENSIONS = ['.cds', '.cds.fasta', '.ffn']

# 缓存读取的序列长度
SEQUENCE_LENGTHS_CACHE = {}

def read_fasta_lengths(fasta_file):
    """
    读取FASTA格式文件，计算并返回每个序列的长度。

    参数:
    fasta_file (str): FASTA格式文件的路径。

    返回:
    dict: 键为序列ID，值为序列长度的字典。
    """
    if fasta_file not in SEQUENCE_LENGTHS_CACHE:
        lengths = {}
        try:
            with open(fasta_file, 'r') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    lengths[record.id] = len(record.seq)
        except Exception as e:
            print(f"读取文件 {fasta_file} 时出错: {e}")
            raise
        SEQUENCE_LENGTHS_CACHE[fasta_file] = lengths
    return SEQUENCE_LENGTHS_CACHE[fasta_file]

def calculate_cds_proportions(genome_lengths, cds_lengths):
    """
    计算每个CDS（编码序列）长度占对应基因组长度的比例。

    参数:
    genome_lengths: 字典，键为基因组ID，值为基因组的长度。
    cds_lengths: 字典，键为CDS ID，值为CDS的长度。

    返回:
    proportions: 字典，键为CDS ID，值为CDS长度占基因组长度的比例。
    """
    proportions = {}
    for cds_id, cds_len in cds_lengths.items():
        if cds_id in genome_lengths:
            genome_len = genome_lengths[cds_id]
            proportions[cds_id] = cds_len / genome_len
        else:
            print(f"警告: 没有找到与CDS ID {cds_id}匹配的基因组序列。")
    return proportions

def write_results_to_file(output_file, results):
    """
    将分析结果写入到文件中。

    参数:
    output_file: 字符串类型，表示输出文件的路径。
    results: 字典类型，包含分析结果。字典的键是文件名，值是另一个字典，
             该字典的键是CDs的ID，值是相应的比例。

    返回值:
    无返回值。该函数将结果以制表符分隔的形式写入到指定的输出文件中。
    """
    try:
        with open(output_file, 'w') as out:
            out.write("File Name\tProportion\n")
            for file_name, proportions in results.items():
                for cds_id, proportion in proportions.items():
                    out.write(f"{file_name}\t{proportion:.4f}\n")
    except Exception as e:
        print(f"写入结果到文件 {output_file} 时出错: {e}")
        raise

def process_file_pair(genome_file, cds_file):
    """
    处理基因组文件和CDS文件对，计算CDS在基因组中的比例。

    参数:
    genome_file: str - 基因组的FASTA文件路径。
    cds_file: str - CDS的FASTA文件路径。

    返回:
    tuple - 包含基因组文件名和CDS比例的元组。
    """
    genome_lengths = read_fasta_lengths(genome_file)
    cds_lengths = read_fasta_lengths(cds_file)
    return (os.path.basename(genome_file), calculate_cds_proportions(genome_lengths, cds_lengths))

def process_file_pairs(directory, output_file):
    """
    处理给定目录中的基因组文件和其对应的编码序列文件对。
    并行处理所有文件对，然后将结果写入指定的输出文件。
    """
    results = {}
    file_pairs = []
    for filename in os.listdir(directory):
        base_name, ext = os.path.splitext(filename)
        if ext.lower() in GENOME_EXTENSIONS:
            genome_file = os.path.join(directory, filename)
            for cds_ext in CDS_EXTENSIONS:
                cds_file = os.path.join(directory, f"{base_name}{cds_ext}")
                if os.path.exists(cds_file):
                    file_pairs.append((genome_file, cds_file))
    
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = {file_name: proportions for file_name, proportions in executor.map(process_file_pair, *zip(*file_pairs))}
    
    write_results_to_file(output_file, results)

def main():
    parser = argparse.ArgumentParser(description="计算CDS长度相对于基因组序列的比例。")
    parser.add_argument("--directory", "-d", help="包含FASTA文件的目录", required=True)
    parser.add_argument("--output_file", "-o", help="输出结果的文件", required=True)
    args = parser.parse_args()

    if not os.path.exists(args.directory):
        print(f"目录 {args.directory} 不存在。")
        return

    process_file_pairs(args.directory, args.output_file)

    print('\nRunning successfully!\n')
    print('Output file: ', args.output_file)
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()