#!/usr/bin/env python
import os
import argparse
from Bio import SeqIO

# 定义支持的基因组和CDS文件后缀
GENOME_EXTENSIONS = ['.fa', '.fasta', '.fna']
CDS_EXTENSIONS = ['.cds', '.cds.fasta', '.ffn']

def read_fasta_lengths(fasta_file):
    """
    读取FASTA文件并返回一个字典，该字典将序列ID映射到长度。
    :param fasta_file: 输入的FASTA文件路径
    :return: 字典，键为序列ID，值为序列长度
    """
    lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        lengths[record.id] = len(record.seq)
    return lengths

def calculate_cds_proportions(genome_lengths, cds_lengths):
    """
    计算CDS长度相对于基因组长度的比例。
    :param genome_lengths: 基因组序列长度的字典
    :param cds_lengths: CDS序列长度的字典
    :return: 字典，键为CDS ID，值为CDS长度与对应基因组长度的比例
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
    将结果写入制表符分隔的文件。
    :param output_file: 输出文件路径
    :param results: 结果字典，键为文件名，值为比例字典
    """
    with open(output_file, 'w') as out:
        out.write("File Name\tProportion\n")  # 写入标题行
        for file_name, proportions in results.items():
            for cds_id, proportion in proportions.items():
                out.write(f"{file_name}\t{proportion:.4f}\n")  # 写入文件名和比例

def process_file_pairs(directory, output_file):
    """
    处理目录中的所有文件对，并将结果写入文件。
    :param directory: 包含FASTA文件的目录路径
    :param output_file: 输出结果的文件路径
    """
    results = {}  # 初始化结果字典
    for filename in os.listdir(directory):
        base_name, ext = os.path.splitext(filename)
        if ext.lower() in GENOME_EXTENSIONS:  # 如果是基因组文件
            genome_file = os.path.join(directory, filename)
            for cds_ext in CDS_EXTENSIONS:  # 遍历所有CDS后缀
                cds_file = os.path.join(directory, f"{base_name}{cds_ext}")
                if os.path.exists(cds_file):  # 如果找到了匹配的CDS文件
                    print(f"正在处理配对: {filename} 和 {cds_file}")
                    genome_lengths = read_fasta_lengths(genome_file)
                    cds_lengths = read_fasta_lengths(cds_file)
                    proportions = calculate_cds_proportions(genome_lengths, cds_lengths)
                    results[filename] = proportions  # 存储结果
                    break  # 找到匹配后停止搜索

    write_results_to_file(output_file, results)  # 将结果写入文件

def main():
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(description="计算CDS长度相对于基因组序列的比例。")
    parser.add_argument("--directory", "-d", help="包含FASTA文件的目录", required=True)
    parser.add_argument("--output_file", "-o", help="输出结果的文件", required=True)
    args = parser.parse_args()

    process_file_pairs(args.directory, args.output_file)  # 调用处理函数

    print('\nRunning successfully!')
    print('Output file: ',args.output_file)
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()  # 运行主函数