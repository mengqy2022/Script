#!/bin/bash/
# **************************************************
# 处理blast结果文件，保留最佳的hit，并转换为excel文件
# Date   : 2024-03-13
# Author : 孟庆瑶
# **************************************************

import pandas as pd
import sys

def usage():
    print('\n >>>> 处理blast结果文件，保留最佳的hit，并转换为excel文件 <<<<\n')
    print('Usage: python3 metabolic_blast_ftm6.py [ASSOCIATED FILES] [BLAST RESULT FILE NAME]\n')
    print('Example: python3 metabolic_blast_ftm6.py associated_files.txt blast_result.txt\n')
    print('    associated_files.txt: 包含基因组注释id与基因组序列id和ncbi id的映射关系的文件\n')
    print('    blast_result.txt: blast结果文件[outfmt 6]\n')
    print('    >>>> mqy <<<<\n')
    sys.exit()

def main(associated_file, blast_file,out_file):
    # 读取基因组注释id与基因组序列id和ncbi id的映射关系
    df = pd.read_csv(associated, sep='\t',header=None)
    df[["prokka_ids","3"]] = df[2].str.split("_",expand=True)
    df.rename(columns={0:'ncbi_ids',1:'tree_ids'},inplace=True)
    df = df[['ncbi_ids','tree_ids','prokka_ids']]

    # 读取blast结果文件
    df1 = pd.read_csv(blast_file, sep='\t',header=None)
    df1[["prokka_ids","ids_count"]] = df1[1].str.split("_",expand=True)
    df1.rename(columns={0: "query_id", 10: "e_evalue"}, inplace=True)
    df1 = df1[["query_id","prokka_ids","ids_count","e_evalue"]]
    
    #  合并两个表
    merged_df = pd.merge(df, df1, on='prokka_ids')
    merged_df = merged_df.sort_values(by=['prokka_ids', 'query_id','e_evalue'])
    merged_df.drop_duplicates(subset=['prokka_ids','query_id'], inplace=True)
    merged_df['cds_ids'] = pd.concat([merged_df['prokka_ids'], merged_df['ids_count']], axis=1).agg('_'.join, axis=1)
    merged_df.drop(['ids_count','prokka_ids'], axis=1, inplace=True)
    merged_df.to_excel(out_file, index=False)
    print('merged_df.xlsx has been generated.')

if __name__ == '__main__': 
    if len(sys.argv) != 3:
        usage()
    associated = sys.argv[1]
    blast_file = sys.argv[2]
    out_file = blast_file.split('.')[0]+'.xlsx'
    main(associated, blast_file, out_file)
    print('Done!\n')
    print('Output file: '+sys.argv[2]+'.xlsx\n')
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')
    sys.exit()


