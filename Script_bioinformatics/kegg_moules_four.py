#!/usr/bin/Python3.6
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com
# @Time      : 2023/4/20 8:57
# @Project   : kegg_moules.py
# @File      : python_moules_test

#  sys模块主要针对于python解释器相关的变量和方法，这点于os模块不同。
import sys

#  sys.argv：获取命令行参数，第一个参数是自己程序名，常用来启动程序时给予某些参数。
if len(sys.argv) != 5:
	print('\nuseage:\npython3 kegg_moules.py <kaas.keg> <gene.list> <moules_anno.result>\nexample:\npython3 kegg_moules.py q00002.keg gene.list.txt moules_anno.txt\n')
	exit()

gene = {}
gene_list = open(sys.argv[2], 'r')
for line in gene_list:
	line = line.strip().split('\t')
	gene[line[1]] = line[0]

gene_list.close()

kegg = {}
moules_anno = open(sys.argv[3], 'w')
moules_anno.write('gene_total\tko1_name\tko2_name\tm_id\tm_name\tmap_id\n')
kaas = open(sys.argv[1], 'r')

for line in kaas:
	line = line.strip()
	if line[0] == 'B' and len(line) > 1:
		ko1_name = line[6:len(line)]
	elif line[0] == 'C' and len(line) > 1:
		ko2_name = line[5:len(line)]
	elif line[0] == 'D' and len(line) > 1:
		m_id = line[7:13]
		m_name = line[15:len(line)]
		if ' [' in m_name:
			map_id = m_name.split(' [')[1]
			m_name = m_name.split(' [')[0]
			map_id = map_id.split('PATH:')[1]
		if m_id in gene:
			moules_anno.write(f'{gene[m_id]}\t{ko1_name}\t{ko2_name}\t{m_id}\t{m_name}\t{map_id}\n')

kaas.close()
moules_anno.close()