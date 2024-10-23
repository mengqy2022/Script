#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

script_name = sys.argv[0]

def usage():
    print('Usage: python {} [fasta_file] [namelist_file] [outfile_name]'.format(script_name))


def main():
    outf = open(sys.argv[3], 'w')
    dict = {}
    with open(sys.argv[1], 'r') as fastaf:
        for line in fastaf:
            if line.startswith('>'):
                name = line
                dict[name] = ''
            else:
                dict[name] += line.replace('\n', '')  # 读取整个fasta文件构成字典

    with open(sys.argv[2], 'r') as listf:
        for row in listf:
            row = row.strip()
            for key in dict.keys():  # 选取包含所需名称的基因名和序列
                if row in key:
                    outf.write(key)
                    outf.write(dict[key] + '\n')
    outf.close()


try:
    main()
except IndexError:
    usage()
