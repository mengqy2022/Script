#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

# 创建一个参数解析器
parser = argparse.ArgumentParser(description='处理 KEGG 注释和基因列表文件',
                                epilog= '\t更详细的信息请访问: https://mengqy2022.github.io/comparative%20genomics/Comparative-genomics/\n')

# 添加参数
parser.add_argument('-keg','--kaas_keg', type=str, help='输入的 KAAS KEG 文件路径', required=True)
parser.add_argument('-ko','--gene_list', type=str, help='输入的基因列表文件路径', required=True)
parser.add_argument('-g','--gene_anno_result', type=str, help='输出的基因注释结果文件路径')
parser.add_argument('-p','--pathway_anno_result', type=str, help='输出的通路注释结果文件路径')

# 解析命令行参数
args = parser.parse_args()

# 读取基因与蛋白的对应关系列表
gene = {}
with open(args.gene_list, 'r') as gene_list:
    for line in gene_list:
        line = line.strip().split('\t')
        gene[line[0]] = line[1]

# 整理 KAAS 注释列表，得到表格样式，添加基因、蛋白和 KEGG 的对应关系等
kegg = {}
with open(args.gene_anno_result, 'w') as gene_anno:
    gene_anno.write('gene_id\tprotein_id\tko4_id\tko4_gene\tEC\tko3_id\tko3_pathway\tko2_id\tko2_name\tko1_id\tko1_name\n')
    
    with open(args.kaas_keg, 'r') as kaas:
        for line in kaas:
            line = line.strip()
            if line[0] == 'A' and len(line) > 1:
                ko1_id = line[1:6]
                ko1_name = line[7:len(line)]
            elif line[0] == 'B' and len(line) > 1:
                ko2_id = line[3:8]
                ko2_name = line[9:len(line)]
            elif line[0] == 'C' and len(line) > 1:
                ko3_id = line[5:10]
                ko3_pathway = line[11:len(line)]
                if ' [' in ko3_pathway:
                    ko3_pathway = ko3_pathway.split(' [')[0]
            elif line[0] == 'D' and len(line) > 1:
                ko_detail = line[7:len(line)].split('; ', 1)
                protein_id = ko_detail[0]
                ko4_id = ko_detail[1][0:6]
                ko_detail = ko_detail[1][8:len(ko_detail[1])]
                if ' [' in ko_detail:
                    ko_detail = ko_detail.split(' [')
                    ko4_gene = ko_detail[0]
                    EC = '[' + ko_detail[1]
                else:
                    ko4_gene = ko_detail
                    EC = ''
                if protein_id in gene:
                    gene_anno.write(f'{gene[protein_id]}\t{protein_id}\t{ko4_id}\t{ko4_gene}\t{EC}\t{ko3_id}\t{ko3_pathway}\t{ko2_id}\t{ko2_name}\t{ko1_id}\t{ko1_name}\n')

# 统计注释到各通路的基因数量，以及编辑通路图链接
pathway = {}
with open(args.gene_anno_result, 'r') as gene_anno:
    gene_anno.readline()  # 跳过表头
    for line in gene_anno:
        line = line.strip().split('\t')
        if line[5] not in pathway:
            ko_class = '\t'.join([line[9], line[10], line[7], line[8], line[5], line[6]])
            pathway[line[5]] = [ko_class, [line[0]], [line[2]]]
        else:
            if line[2] not in pathway[line[5]][2]:
                pathway[line[5]][2].append(line[2])
            if line[0] not in pathway[line[5]][1]:
                pathway[line[5]][1].append(line[0])

# 生成通路注释结果文件
with open(args.pathway_anno_result, 'w') as pathway_anno:
    pathway_anno.write('ko1_id\tko1_name\tko2_id\tko2_name\tko3_id\tko3_pathway\tgene_number\tpathway_link\n')
    for key, values in pathway.items():
        gene_number = len(values[1])
        pathway_link = f'https://www.genome.jp/kegg-bin/show_pathway?ko{key}/reference%3dwhite/'
        for ko in values[2]:
            pathway_link += f'{ko}%09orange/'

        pathway_anno.write(f'{values[0]}\t{gene_number}\t{pathway_link}\n')

print(f'处理完成，输出结果保存在: {args.gene_anno_result} 和 {args.pathway_anno_result}')
