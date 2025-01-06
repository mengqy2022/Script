#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/12/25
# @Description : 处理比对结果，去除对应的查询序列ID

import pandas as pd
import argparse
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description='处理比对结果，去除对应的查询序列ID')
    parser.add_argument('-i', '--input_file', type=str, help='输入的比对结果文件路径', required=True)
    parser.add_argument('-l', '--unique_ids', type=str, help='输出的结果文件路径', default='unique_query_ids.txt')
    parser.add_argument('-s', '--search_file', type=str, help='用于搜索的文件路径', required=True)
    parser.add_argument('-r', '--result_file', type=str, help='保存搜索结果的文件路径', required=True)
    return parser.parse_args()

def filter_sequences(input_file, unique_ids, search_file, result_file):
    try:
        # 读取比对结果
        df = pd.read_csv(input_file, sep="\t", header=None)
        df.columns = ['query_id', 'subject_id', 'identity', 'length', 'mismatches', 
                      'gap_opens', 'start', 'end', 'subject_start', 'subject_end', 
                      'evalue', 'score']
        
        # 创建一个集合来保存需要去除的查询序列IDs
        ids_to_remove = {query_id for query_id, group in df.groupby('query_id') if len(group) > 1 and group['identity'].iloc[1] >= 70}

        # 筛选出不在ids_to_remove中的查询序列
        filtered_df = df[~df['query_id'].isin(ids_to_remove)]
        
        # 输出过滤后的结果
        filtered_df.to_csv(unique_ids, sep="\t", index=False, header=False)
        unique_query_ids = filtered_df['query_id'].unique()

        with open(unique_ids, 'w') as f:
            for query_id in unique_query_ids:
                f.write(f"{query_id}\n") 

        # 从搜索文件中提取包含唯一query_id的数据
        search_df = pd.read_csv(search_file, sep="\t", header=None)
        search_df.columns = ['contig', 'type', 'feature', 'start', 'end', 
                             'score', 'strand', 'phase', 'attributes']

        # 根据'query_id'在attributes列中搜索匹配项
        matching_rows = search_df[search_df['attributes'].str.contains('|'.join(unique_query_ids), na=False)]
        
        # 保存匹配结果到新的文件
        matching_rows.to_csv(result_file, sep="\t", index=False, header=False)
        print(f"匹配结果已保存到: {result_file}")
    
    except FileNotFoundError as e:
        print(f"错误：文件未找到 - {e}")
    except pd.errors.EmptyDataError:
        print("错误：输入文件为空")
    except Exception as e:
        print(f"发生了一个错误: {e}")

def main():
    args = parse_arguments()
    filter_sequences(args.input_file, args.unique_ids, args.search_file, args.result_file)

if __name__ == '__main__':
    main()
