#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from datetime import datetime

def usage():
    print("------------------------------------------------------------------------------------------")
    print('Description: 将blast结果与GO信息关联。\n')
    print('Usage: python3 go_ann_setp_1.py [diamond_blast.txt] [diamond_blast.add_go.txt] [idmapping.tb]\n')
    print('Usage: python3 go_ann_setp_1.py [INPUT.TXT] [OUTP.TXT] [idmapping.tb]\n')
    print(" >>>> mqy <<<<\n")
    print('\t更详细的信息请访问: https://mengqy2022.github.io/genomics/GO/\n')

def main():
    input = sys.argv[1]
    output = sys.argv[2]
    go_db = sys.argv[3]
    
    #唯一蛋白id
    protein = {}
    
    xls = open(input, 'r')
    xls.readline()
    for line in xls:
            #  NR id 位于第 6 列（即“Hit_name”）
            id = line.split('\t')[5]
            if id not in protein:
                    protein[id] = 1
    
    xls.close()
    print(len(protein))
    
    #选择go
    tmp = open('go_tmp.xls', 'w')
    go = open(go_db, 'r')
    for line in go:
            line = line.strip().split('\t')
            #  NR id 位于第 4 列（即“RefSeq”）
            #  GO将用;分割
            nr = line[3].split('; ')
            for id in nr:
                    #  GO:0046782
                    #  len(line) > 8 列数大于8
                    if id in protein and len(line) > 8:
                            #  GO id 位于第 8 列
                            print(f'{id}\t{line[7]}', file = tmp)
                            #  RefSeq \t GO id
    
    tmp.close()
    
    #合并
    go = {}
    tmp = open('go_tmp.xls', 'r')
    for line in tmp:
            line = line.strip().split('\t')
            if len(line) == 2:
                    go[line[0]] = line[1]
            elif len(line) == 1:
                    go[line[0]] = ''
    
    tmp.close()
    
    xls2 = open(output, 'w')
    xls = open(input, 'r')
    print(f'{xls.readline().strip()}\tGO', file = xls2)
    
    for line in xls:
            line = line.strip()
            id = line.split('\t')[5]
            if id in go:
                    print(f'{line.strip()}\t{go[id]}', file = xls2)
    
    xls.close()
    xls2.close()
    os.system('rm go_tmp.xls')

    now = datetime.now()
    print('Done!\n')
    print('Time: '+ now.strftime('%Y-%m-%d %H:%M:%S')+'\n')
    print('Output file: '+sys.argv[2]+'\n')
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

try:
    main()
except IndexError:
    usage()