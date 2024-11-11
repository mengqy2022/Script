#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from datetime import datetime

def usage():
    print("\nDescription: 将blast --outfmt 0 变换格式")
    print("--------------------------------------------------------------------")
    print(f'Usage: python3 {sys.argv[0]} [infile.fa] [outfile.fa] \n')
    print(" >>>> mqy <<<<\n")
    print('\t更详细的信息请访问: https://mengqy2022.github.io/genomics/GO/\n')

def main():
    if len(sys.argv) != 2 or sys.argv[1] == "-h" or sys.argv[1] == "--help" or sys.argv[1] == "--help" or sys.argv[1] == "-help":
        usage()
        sys.exit(1)

    xls = open(sys.argv[2], 'w')
    print('Query_name\tQuery_description\tQuery_length\tQuery_start\tQuery_end\tHit_name\tHit_description\tHit_length\tHit_start\tHit_end\tGap\tIdentity\tEvalue', file = xls)

    #  读取blast 0 格式比对结果
    blast = open(sys.argv[1], 'r')

    n = 1
    frame = 0
    id = ''
    Query_name = ''
    Query_description = ''
    Query_length = ''
    Query_start = []
    Query_end = []
    Hit_name = ''
    Hit_description = ''
    Hit_length = ''
    Hit_start = []
    Hit_end = []
    Gap = ''
    Identity = ''
    Evalue = ''

    for line in blast:
        line = line.strip()
        #  Query= APKFPFOI_1_pseudo_00001 APKFPFOI_1 [7907:8158](-)
        if 'Query=' in line:
            #  id不能重复 使用 !=
            #  APKFPFOI_1_pseudo_00001 APKFPFOI_1 [7907:8158](-)
            #  APKFPFOI_1_pseudo_00001
            if id and line.split('Query= ')[1].split()[0] != id:
                if frame > 0:
                    Query_start = min(Query_start)
                    Query_end = max(Query_end)
                    print(f'{Query_name}\t{Query_description}\t{Query_length}\t{Query_start}\t{Query_end}\t{Hit_name}\t{Hit_description}\t{Hit_length}\t{min(Hit_start)}\t{max(Hit_end)}\t{Gap}\t{Identity}\t{Evalue}', file = xls)
                elif frame < 0:
                    Query_start = max(Query_start)
                    Query_end = min(Query_end)
                    print(f'{Query_name}\t{Query_description}\t{Query_length}\t{Query_start}\t{Query_end}\t{Hit_name}\t{Hit_description}\t{Hit_length}\t{min(Hit_start)}\t{max(Hit_end)}\t{Gap}\t{Identity}\t{Evalue}', file = xls)
            Query_name = line.split('Query= ')[1].split()[0]
            description = line.split('Query= ')[1].split(maxsplit = 1)
            Query_description = description[1] if len(description) > 1 else ''
            id = Query_name
            Hit_start = []
            Hit_end = []
            Query_start = []
            Query_end = []
            n = 1
        if 'Length=' in line:
            if n == 1:
                Query_length = int(line.split('Length=')[1])
            if n == 2:
                Hit_length = int(line.split('Length=')[1])
        if '>' in line:
            n += 1
        if '>' in line and n == 2:
            Hit_name = line.split('>')[1].split()[0]
            description = line.split('>')[1].split(maxsplit = 1)
            Hit_description = description[1] if len(description) > 1 else ''
        if 'Identities = ' in line and n == 2:
            Identity = line.split('Identities = ')[1].split(', ')[0]
        if 'Expect = ' in line and n == 2:
            Evalue = line.split('Expect = ')[1]
        if 'Query  ' in line and n == 2:
            Query_start.append(int(line.split('Query  ')[1].split()[0]))
            Query_end.append(int(line.split()[-1]))
        if 'Sbjct  ' in line and n == 2:
            Hit_start.append(int(line.split('Sbjct  ')[1].split()[0]))
            Hit_end.append(int(line.split()[-1]))
        if 'Gaps = ' in line and n == 2:
            Gap = line.split('Gaps = ')[1]
        if 'No hits found' in line:
            id = ''
        if 'Frame = ' in line and n == 2:
            frame = int(line.split('Frame = ')[1])
    
    if id:
        if frame > 0:
            Query_start = min(Query_start)
            Query_end = max(Query_end)
            print(f'{Query_name}\t{Query_description}\t{Query_length}\t{Query_start}\t{Query_end}\t{Hit_name}\t{Hit_description}\t{Hit_length}\t{min(Hit_start)}\t{max(Hit_end)}\t{Gap}\t{Identity}\t{Evalue}', file = xls)
        elif frame < 0:
            Query_start = max(Query_start)
            Query_end = min(Query_end)
            print(f'{Query_name}\t{Query_description}\t{Query_length}\t{Query_start}\t{Query_end}\t{Hit_name}\t{Hit_description}\t{Hit_length}\t{min(Hit_start)}\t{max(Hit_end)}\t{Gap}\t{Identity}\t{Evalue}', file = xls)
    
    blast.close()
    xls.close()

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