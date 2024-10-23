#!/bin/bash/
# **************************************************
# 将序列格式文件转换为以制表符分割的csv文件
# Date   : 2024-03-13
# Author : 孟庆瑶
# **************************************************

from Bio import SeqIO
import pandas as pd
import sys

def usage():
    print('\n将序列格式文件转换为以制表符分割的csv文件\n')
    print(f'Usage: python3 {sys.argv[0]} [FILE NAME] [OUTPUT PERFIX]\n')
    print(f'Example: python3 {sys.argv[0]} reads.fasta output\n')
    print('注意：输入文件必须是fasta格式，输出文件将以.csv为后缀名')
    #print('注意：目前只支持序列名称header中包含“_”的情况，必须为“metabolism_EZids_class”格式')
    print('注意：输出文件将会保存在当前目录下，请不要在文件名中出现空格！')
    print('注意：输出文件将会覆盖原有文件！\n')
    print('    >>>> mqy <<<<\n')
    sys.exit()

def read_fasta(fasta):
    #  创建空的字典
    fa_dict = {}
    with open(fasta,'r') as fa:
        for seq in SeqIO.parse(fa, "fasta"):
            sequence = str(seq.seq).strip()
            seq_id = str(seq.id)
            fa_dict[seq_id] = sequence
    return fa_dict

def create_dataframe(fast_dict):
    dataset = pd.DataFrame.from_dict(fast_dict, orient='index', columns=['sequence'])
    dataset.index.name = 'IDS'
    dataset.reset_index(inplace=True)
    #dataset[["metabolism","EZ_ids","class"]] = dataset["IDS"].str.split("_",expand=True)
    #dataset = dataset[["metabolism","EZ_ids","class","sequence"]]
    #dataset.drop(columns=["id"],inplace=True)
    dataset.to_csv(sys.argv[2]+'.csv', index=False,sep='\t')
    #return dataset

def main():
    if __name__ == '__main__':
        #  检查参数
        # 脚本 1 参数 2 参数 3
        if len(sys.argv)!= 3:
            usage()
        create_dataframe(read_fasta(sys.argv[1]))
        print('\nRunning successfully!')
        print('Output file: '+sys.argv[2]+'.csv\n')
        print('\t>>> mqy <<<<\n')
        print('\t如果有任何问题请及时联系\n')
        print('\t邮箱：<15877464851@163.com>\n')
        sys.exit()

try:
    main()
except IndexError:
    usage()
