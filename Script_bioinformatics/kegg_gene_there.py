#!/usr/bin/Python3.6
# -*- coding: utf-8 -*-

#  sys模块主要针对于python解释器相关的变量和方法，这点于os模块不同。
import sys

#  sys.argv：获取命令行参数，第一个参数是自己程序名，常用来启动程序时给予某些参数。
if len(sys.argv) != 5:
	print('\nuseage:\npython3 kegg_trans.py <kaas.keg> <gene.list> <gene_anno.result> <pathway_anno.result>\nexample:\npython3 kegg_trans.py q00001.keg gene.list.txt gene_anno.txt pathway_anno.txt\n')
	exit()

# open() 函数的作用是打开一个文件，并返回一个 file对象(即文件对象)。
# open 是一个动作，可以理解为我们打开文档的点击动作。
# file 对象是一个实物，可以理解为我们打开的具体文档，例如记事本、表格、Word 或其他具体的文档。
# 读取 gene 和 protein 的对应关系列表
# strip() 方法用于移除字符串头尾指定的字符（默认为空格）或字符序列。
# 注意：该方法只能删除开头或是结尾的字符，不能删除中间部分的字符。
# Python中Split函数可以对字符串按照指定规则分割，并将分割后的字段作为list返回的函数。
gene = {}
gene_list = open(sys.argv[2], 'r')
for line in gene_list:
	line = line.strip().split('\t')
	gene[line[0]] = line[1]

gene_list.close()

#整理 KAAS 注释列表，得到表格样式，添加基因、蛋白和 KEGG 的对应关系等
kegg = {}
gene_anno = open(sys.argv[3], 'w')
gene_anno.write('gene_id\tprotein_id\tko4_id\tko4_gene\tEC\tko3_id\tko3_pathway\tko2_id\tko2_name\tko1_id\tko1_name\n')
kaas = open(sys.argv[1], 'r')

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

kaas.close()
gene_anno.close()

#统计注释到各通路的基因数量，以及编辑通路图链接
gene_anno = open(sys.argv[3], 'r')
gene_anno.readline()
pathway = {}
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

gene_anno.close()

pathway_anno = open(sys.argv[4], 'w')
pathway_anno.write('ko1_id\tko1_name\tko2_id\tko2_name\tko3_id\tko3_pathway\tgene_number\tpathway_link\n')
for key,values in pathway.items():
	gene_number = len(values[1])
	pathway_link = f'https://www.genome.jp/kegg-bin/show_pathway?ko{key}/reference%3dwhite/'
	for ko in values[2]:
		pathway_link += f'{ko}%09orange/'
	
	pathway_anno.write(f'{values[0]}\t{gene_number}\t{pathway_link}\n')

pathway_anno.close()
