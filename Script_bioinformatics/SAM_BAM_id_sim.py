#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import pysam
import re

def usage():
    print('Usage: python3 fasta_id_sim.py [input_file] > [outfile]')
    print('descriptive: Unified processing simplifies sequence ids')

#  Pysam 是一个 python 模块，可以轻松读取和操作存储在 SAM/BAM 文件中的映射短读取序列数据。

# 任何一门编程语言中，文件的输入输出、数据库的连接断开等，都是很常见的资源管理操作。
# 但资源都是有限的，在写程序时，必须保证这些资源在使用过后得到释放，不然就容易造成资源泄露，轻者使得系统处理缓慢，严重时会使系统崩溃。

# 例如，前面在介绍文件操作时，一直强调打开的文件最后一定要关闭，否则会程序的运行造成意想不到的隐患。
# 但是，即便使用 close() 做好了关闭文件的操作，如果在打开文件或文件操作过程中抛出了异常，还是无法及时关闭文件。

# 为了更好地避免此类问题，不同的编程语言都引入了不同的机制。
# 在 Python 中，对应的解决方式是使用 with as 语句操作上下文管理器（context manager），它能够帮助我们自动分配并且释放资源。 

def main():
    if len(sys.argv) != 2 or sys.argv[1] == "-h" or sys.argv[1] == "--help" or sys.argv[1] == "--help" or sys.argv[1] == "-help":
        usage()
        sys.exit(1)
    with pysam.FastxFile(sys.argv[1]) as fh:
        for r in fh:
            #  .不能放最后面
            new_name = re.split(r'[:;| ]',r.name)[0]
            print(">"+new_name)
            print(r.sequence)

try:
    main()
except IndexError:
    usage()