# /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import pyfastx
import argparse
import os
from multiprocessing import Pool,freeze_support

def removePolyA(seq):
    reverse = seq[::-1]  # 将输入序列反转，方便从末尾开始处理
    Astate, Astretch, Vstretch, trimPos = False, 0, 0, 0  # 初始化状态和计数变量
    for pos in range(0, len(reverse), 1):  # 遍历反转后的序列的每个位置
        base = reverse[pos]  # 获取当前位置的碱基
        if not Astate:  # 初始设定不处于PolyA位置，如果检测到6个连续A则进入Astate
            if base in 'Aa':  # 如果当前碱基是 'A' 或 'a'
                Astretch += 1  # 记录连续PolyA长度
                if Astretch == 6:  # 如果连续PolyA长度达到 6
                    Astate = True  # 进入PolyA位置
                    lastA = pos  # 记录最后一个 'A' 的位置
            else:
                Astretch = 0  # 重置连续PolyA长度为 0
            if pos > 20:  # 如果20个碱基中没有检测到连续6个A，则终止，该序列没有PolyA
                break  # 结束循环，不再处理后面的碱基
        if Astate:  # 如果处于PolyA位置，找到polyA的终止
            if not (base in 'Aa'):  # 如果当前碱基不是 'A' 或 'a'
                Vstretch += 1  # 记录连续非PolyA长度
                Astretch = 0  # 重置连续PolyA长度为 0
            else:
                Astretch += 1  # 记录连续PolyA长度
                if Astretch >= 3:  # 如果连续PolyA长度达到 3
                    Vstretch = 0  # 重置连续非PolyA长度为 0
                    lastA = pos  # 更新最后一个 'A' 的位置
            if Vstretch >= 3:  # 如果连续非PolyA长度达到 3
                trimPos = lastA  # 记录需要截断的位置
                break  # 结束循环，不再处理后面的碱基
        if pos == len(reverse)-1:  # 如果PolyA一直判定到序列结束，说明该序列全部为PolyA
                trimPos = pos  # 截断位置就是序列长度
    if (len(reverse)-trimPos) < 50 : #如果trim后剩余的碱基数小于50，则放弃该reads
        return None,None # 返回空值
    else:
        reverseTrim = reverse[trimPos:]  # 根据截断位置截取序列的一部分
        seqTrim = reverseTrim[::-1]  # 将截取的序列反转回来
        return seqTrim,trimPos  # 返回去除PolyA后的序列和trim的长度


def process_file(args):
    inFile, output_dir = args
    baseName = os.path.splitext(os.path.basename(inFile))[0]
    outFileName = os.path.join(output_dir, baseName + '_trimmedPolyA.fa')
    with open(outFileName, 'w', newline='\n') as out:
        for name, seq in pyfastx.Fastx(inFile):
            seqTrim, trimPos = removePolyA(seq)
            if seqTrim is not None:
                out.write('>' + name + '\n' + seqTrim + '\n')

# 自定义函数，用于将逗号分割的字符串转换为列表，因为TBtools传入的文件名是逗号分隔的
def parse_file_list(file_list_str):
    return file_list_str.split(',')

if __name__ == '__main__':
    freeze_support()
    parser = argparse.ArgumentParser()
    parser.add_argument('--inFiles', '-i', type=parse_file_list,  help='files to be adapter and poly trimmed')
    parser.add_argument('--processes', '-p', type=int, default=1, help='number of processes to use')
    parser.add_argument('--outputDir', '-o', type=str, default='output', help='output directory for the processed files')
    args = parser.parse_known_args() 

    inFiles = args[0].inFiles
    num_processes = args[0].processes
    output_dir = args[0].outputDir

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 使用指定数量的进程并行处理输入的多个文件
    with Pool(num_processes) as pool:
        pool.map(process_file, [(inFile, output_dir) for inFile in inFiles])

