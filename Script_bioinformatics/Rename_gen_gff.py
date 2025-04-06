#!/usr/bin/env python3
"""
脚本功能：
1. 重命名fasta文件中的序列ID为Contig1, Contig2...格式
2. 同步更新gff3文件中对应的旧ID
"""

import argparse
from Bio import SeqIO
import gzip
import re

def rename_fasta(input_fasta, output_fasta):
    """
    重命名fasta文件中的序列ID
    返回新旧ID的映射字典
    """
    id_mapping = {}
    with open(output_fasta, 'w') as out_handle:
        for i, record in enumerate(SeqIO.parse(input_fasta, 'fasta'), 1):
            old_id = record.id
            new_id = f"Contig{i}"
            id_mapping[old_id] = new_id
            record.id = new_id
            record.description = ''  # 清空描述信息
            SeqIO.write(record, out_handle, 'fasta')
    return id_mapping

def update_gff3(input_gff3, output_gff3, id_mapping):
    """
    根据ID映射关系更新gff3文件
    """
    with open(input_gff3, 'r') as in_handle, open(output_gff3, 'w') as out_handle:
        for line in in_handle:
            if line.startswith('#'):
                # 保留注释行不变
                out_handle.write(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                out_handle.write(line)
                continue
            
            old_seqid = fields[0]
            if old_seqid in id_mapping:
                fields[0] = id_mapping[old_seqid]
            
            # 处理第9列属性中可能包含的旧ID
            attributes = fields[8]
            for old_id, new_id in id_mapping.items():
                attributes = attributes.replace(old_id, new_id)
            fields[8] = attributes
            
            out_handle.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description='重命名fasta序列ID并更新gff3文件中的对应ID')
    parser.add_argument('-f', '--fasta', required=True, help='输入的fasta文件')
    parser.add_argument('-g', '--gff3', required=True, help='输入的gff3文件')
    parser.add_argument('-o', '--output_fasta', required=True, help='输出的重命名后的fasta文件')
    parser.add_argument('-p', '--output_gff3', required=True, help='输出的更新后的gff3文件')
    
    args = parser.parse_args()
    
    print("开始处理fasta文件...")
    id_mapping = rename_fasta(args.fasta, args.output_fasta)
    print(f"成功重命名 {len(id_mapping)} 条序列")
    
    print("开始更新gff3文件...")
    update_gff3(args.gff3, args.output_gff3, id_mapping)
    print("gff3文件更新完成")
    
    print("所有操作完成！")

if __name__ == '__main__':
    main()