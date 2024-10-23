#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com
# @Time      : 2023/4/20 8:57
# @Project   : kegg_moules.py
# @File      : python_moules_test

import sys

def descriptive():
	print('descriptive: Categorizing KEGG modules into matrices')

def usage():
	print(f'Usage: python3 {sys.argv[0]} [input_file] [outfile]')

def main():

	moules_anno = open(sys.argv[2], 'wt')
	moules_anno.write('B_class\tC_class\tm_id\tm_name\n')
	kaas = open(sys.argv[1], 'rt')	

	for line in kaas:
		line = line.strip()
		if line[0] == 'B' and len(line) > 1:
			B_class = line[4:len(line)]
			#print(B_class)
		elif line[0] == 'C' and len(line) > 1:
			C_class  = line[5:len(line)]
			#print(ko2_name)
		elif line[0] == 'D' and len(line) > 1:
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