#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com上
# @Time      : 2023/4/20 8:57
# @Project   : kegg_moules.py
# @File      : python_moules_test

# gene = {}
# gene_list = open("K_Module.txt", 'r')
# for line in gene_list:
# 	line = line.strip().split('\t')
# 	gene[line[1]] = line[0]
# gene_list.close()

import sys

def descriptive():
    print('\n描述: 将KEGG网页的模块分类自制文件输入，变成矩阵形式。\n')

def usage():
    print('Usage: python3 Module_name_class.py [input_file] [outfile]')

def main():
	kaas = open(sys.argv[1], 'rt')  #  只读数据
	moules_anno = open(sys.argv[2], 'wt')  #  只写数据
	moules_anno.write('B_class\tC_class\tm_id\tm_name\n')

	for line in kaas:
		line = line.strip()
		if line[0] == 'A' and len(line) > 1:
			B_class = line[4:len(line)]
			#print(B_class)
		elif line[0] == 'B' and len(line) > 1:
			C_class  = line[5:len(line)]
			#print(ko2_name)
		elif line[0] == 'C' and len(line) > 1:
			m_id = line[6:12]
			#print(m_id)
			m_name = line[13:len(line)]
			#print(m_name)
			moules_anno.write(f'{B_class}\t{C_class}\t{m_id}\t{m_name}\n')

	kaas.close()
	moules_anno.close()

try:
	main()
except IndexError:
	descriptive()
	usage()