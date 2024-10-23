#!/usr/bin/python
# -*- coding: utf-8 -*-

import filetype #  filetype库精准判断文件类型
import gzip
import argparse

#  创建一个解析对象
parser = argparse.ArgumentParser(description='filter reads from fastq')
#  然后向该对象中添加要关注的命令行参数和选项，每一个add_argument方法对应一个要关注的参数或选项；
parser.add_argument('--fastq', '-q', dest='fastq', help='input a fastq file')
parser.add_argument('--idlist', '-i', dest='idlist', help='input idlist file')
parser.add_argument('--outfile', '-o', dest='outfile', help='input outfile name')
#  解析，解析成功之后即可使用。
args = parser.parse_args()

#  创建空的字符串
fastqdict = {}
kind1 = filetype.guess(args.fastq)  #  识别文件类型
kind2 = filetype.guess(args.idlist)

if kind1 is None:
    with open(args.fastq,'r') as fastq:
        for line in fastq:
            #  用于检查字符串是否以指定的前缀开始
            if line.startswith('@'):
                #  strip 去除字符串两端的空白字符
                #  split  可以用于将一个字符串按照指定的分隔符进行分割，返回一个包含分割后字符串的列表。
                fastqid = line.strip().split()[0][1:]
                fastqdict[fastqid] = ''
            else:
                fastqdict[fastqid] += line
elif kind1.extension == 'gz':
    #  r表示只读，b表示二进制与此对应的是w表示可写，t表示文本方式打开。
    with gzip.open(args.fastq, 'rb') as fastq:
        for line in fastq:
            if line.decode().startswith('@'):
                fastqid = line.decode().strip().split()[0][1:]
                fastqdict[fastqid] = ''
            else:
                #  += 将序列和分数更到fastqid对应的值，形成字典
                fastqdict[fastqid] += line.decode()

if kind1 is None:
    outfile = open(args.outfile, 'w')
elif kind1.extension == 'gz':
    outfile = gzip.open(args.outfile, 'wb')

if kind2 is None:
    with open(args.idlist, 'r') as idfile:
        for line in idfile:
            readsid = line.strip()
            #  是被字典的键
            for key in fastqdict.keys():
                if readsid == key:
                    res = '@' + key + '\n' + fastqdict[key]
                    #  字典的写入可以将键和值一起写入
                    if kind1 is None:
                        outfile.write(res)
                    elif kind1.extension == 'gz':
                        #  将目标字符串str编写为目标二进制数据bytes类型
                        outfile.write(res.encode())
elif kind2.extension == 'gz':
    with gzip.open(args.idlist, 'rb') as idfile:
        for line in idfile:
            readsid = line.decode().strip()
            for key in fastqdict.keys():
                if readsid == key:
                    res = '@' + key + '\n' + fastqdict[key]
                    if kind1 is None:
                        outfile.write(res)
                    elif kind1.extension == 'gz':
                        outfile.write(res.encode())

outfile.close()

