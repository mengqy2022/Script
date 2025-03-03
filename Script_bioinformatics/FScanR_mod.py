#！/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/02/01
# @Author  : mqy
# @File    : FScanR_mod.py

import argparse
import pandas as pd

def print_colored(text, color):
    color_codes = {
        'purple': '\033[95m',
        'green': '\033[92m',
        'red': '\033[91m',
        'reset': '\033[0m'
    }
    print(f"{color_codes.get(color, '')}{text}{color_codes['reset']}")

def main():
    parser = argparse.ArgumentParser(description='Modification of programmed ribosomal translocation (PRF) blast results.')
    
    parser.add_argument('-f', '--input', type=str, help='Path to the input TSV file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file.', required=False, default='RPF.txt')

    args = parser.parse_args()
    file_path = args.input
    output_path = args.output

    try:
        df = pd.read_csv(file_path, delimiter='\t', encoding='utf-8', dtype={'qseqid': str, 'sseqid': str})

        # 对 'qseqid' 进行分组，并提取每个分组中第一个 'sseqid' 值
        grouped_df = df.groupby('qseqid').agg({'sseqid': 'first'}).reset_index()

        # 重命名 'sseqid' 列以便后续使用
        grouped_df.rename(columns={'sseqid': 'first_sseqid'}, inplace=True)

        # 将提取的 'first_sseqid' 与原始 DataFrame 合并，以便筛选
        merged_df = df.merge(grouped_df, on='qseqid', how='inner')

        # 筛选出 'sseqid' 是每个 'qseqid' 组第一个 'sseqid' 的记录
        first_sseqid_df = merged_df[merged_df['sseqid'] == merged_df['first_sseqid']]

        # 对 'first_sseqid' 进行分组，并计算每组的记录数量
        first_sseqid_counts = first_sseqid_df.groupby('first_sseqid').size()

        # 筛选数量大于1的 'first_sseqid' 组合
        filtered_first_sseqids = first_sseqid_counts[first_sseqid_counts > 1].index

        # 筛选出原始 DataFrame 中 'sseqid' 是 filtered_first_sseqids 中的记录
        final_filtered_df = df[df['sseqid'].isin(filtered_first_sseqids)]

        # 创建结果缓冲区
        result_buffer = []

        # 添加列名
        result_buffer.append("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqframe\tsframe\n")

        for sseqid, group in final_filtered_df.groupby('sseqid'):
            # 将数量大于1的 'sseqid' 的所有记录添加到结果缓冲区
            result_buffer.append(group.to_csv(sep='\t', index=False, header=False))

        # 将结果缓冲区的内容一次性写入文件
        with open(output_path, 'w') as f:
            f.writelines(result_buffer)

    except Exception as e:
        print(f"读取文件时出错: {e}")


    except FileNotFoundError:
        print_colored(f"文件 {file_path} 未找到，请检查文件路径。", 'red')
    except pd.errors.EmptyDataError:
        print_colored(f"文件 {file_path} 内容为空。", 'red')
    except pd.errors.ParserError as e:
        print_colored(f"解析文件 {file_path} 时出错，请检查文件格式。详细错误信息: {e}", 'red')
    except Exception as e:
        print_colored(f"发生错误: {e}", 'red')

if __name__ == "__main__":
    main()
