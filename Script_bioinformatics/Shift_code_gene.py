#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys

def usage():
    print('Description: 将序列和开放阅读框整理，输出列表')
    print('Usage: python3 shift_code_gene.py [input file] > [out file]')

def main():
    blastx=open(sys.argv[1], 'r')  #读文件
    list=[]
    for i in blastx:
        i = i.lstrip() # 删除左边的空格
        if i.startswith('>'):
            demo=i.split()[0]
            print(demo)
        elif i.startswith('Frame'):
            i = i.split('\n')[0]
            print(i)

try:
    main()
except IndexError:
    usage()