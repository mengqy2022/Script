#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from datetime import datetime

def usage():
    print("------------------------------------------------------------------------------------------")
    print('Description: 将关联后的结果简化，并添加GO_trem。\n')
    print('Usage: python3 go_ann_setp_2.py [diamond_blast.add_go.txt] [go_anno.txt] [go.obo]\n')
    print('Usage: python3 go_ann_setp_2.py [INPUT.TXT] [OUTP.TXT] [go.obo]\n')
    print(" >>>> mqy <<<<\n")
    print('\t更详细的信息请访问: https://mengqy2022.github.io/genomics/GO/\n')


def main():
    input = sys.argv[1]
    output = sys.argv[2]
    go_db = sys.argv[3]
    
    #
    go = {}
    go_db = open(go_db, 'r')
    for line in go_db:
            line = line.strip()
            if line[0:4] == 'id: ':
                    id = line.split('id: ')[1]
            if line[0:6] == 'name: ':
                    name = line.split('name: ')[1]
            if line[0:11] == 'namespace: ':
                    namespace = line.split('namespace: ')[1]
                    go[id] = [name, namespace]
    
    go_db.close()
    
    #
    output = open(output, 'w')
    print('Query_name\tGO_id\tGO_name\tGO_ontology', file = output)
    
    input = open(input, 'r')
    input.readline()
    for line in input:
            line = line.strip().split('\t')
            if len(line) == 14:
                    #  line[13]是GO号
                    for id in line[13].split('; '):
                            if id in go:
                                    print(f'{line[0]}\t{id}\t{go[id][0]}\t{go[id][1]}', file = output)
    
    input.close()
    output.close()

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