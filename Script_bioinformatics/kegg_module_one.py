#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com
# @Time      : 2023/4/20 8:57

import sys

class KeggModuleProcessor:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.B_class = ""
        self.C_class = ""

    def descriptive(self):
        print('\n描述: 将KEGG网页的模块分类自制文件输入，变成矩阵形式。\n')
        print('\t更详细的信息请访问: https://mengqy2022.github.io/comparative%20genomics/Comparative-genomics/\n')

    def usage(self):
        print('Usage: python3 Kegg_module_name_class_one.py [input_file] [outfile]')

    def process_file(self):
        with open(self.input_file, 'rt') as kaas, open(self.output_file, 'wt') as moules_anno:
            moules_anno.write('B_class\tC_class\tm_id\tm_name\n')

            for line in kaas:
                line = line.strip()
                if line[0] == 'A' and len(line) > 1:
                    self.B_class = line[3:len(line)]
                elif line[0] == 'B' and len(line) > 1:
                    self.C_class = line[4:len(line)]
                elif line[0] == 'C' and len(line) > 1:
                    m_id = line[6:12]
                    m_name = line[13:len(line)]
                    moules_anno.write(f'{self.B_class}\t{self.C_class}\t{m_id}\t{m_name}\n')

def main():
    try:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        processor = KeggModuleProcessor(input_file, output_file)
        processor.process_file()
    except IndexError:
        processor = KeggModuleProcessor("", "")
        processor.descriptive()
        processor.usage()

if __name__ == '__main__':
    main()
