#!/usr/bin/env python3
# **************************************************
# 将数据转换为数据框
# Date   : 2024-09-13
# Author : 孟庆瑶
# **************************************************

# 在服务器里处理一下文件
# 0格式的blast文件
#cat last.blastx | grep "Query=" -A 10  | grep -v -e "Length" -e "Sequences" -e "Score" -e "--" -e ">" -e "Lambda" -e "\*\*\*"  -e "0.318" > last_mod.blastx

import argparse
import pandas as pd
import re

def process_and_write(input_file, output_file):
    """
    读取文件，提取数据并转换为数据框，然后将数据框写入输出文件。
    """
    # 读取文件内容
    with open(input_file, 'r') as file:
        content = file.read()

    # 使用正则表达式提取每个Query及其对应的数据
    pattern = re.compile(r'Query= (TRINITY_\w+)\s+(.*?)(?=Query= |$)', re.DOTALL)
    matches = pattern.findall(content)

    # 初始化一个空的列表来存储数据
    data = []

    # 遍历每个匹配项
    for query, details in matches:
        # 提取每个Query的详细信息
        lines = details.strip().split('\n')
        for line in lines:
            if line.strip():
                parts = line.split(maxsplit=2)
                if len(parts) == 3:
                    accession, score, e_value = parts
                    # 提取完整的Accession信息
                    full_accession = line.split()[0] + ' ' + line.split()[1] + ' ' + line.split()[2]
                    data.append([query, full_accession, score, e_value])

    # 将数据转换为数据框
    df = pd.DataFrame(data, columns=['Query', 'Accession', 'Score', 'E-value'])

    # 将数据框写入输出文件
    df.to_csv(output_file, sep='\t',index=False)

def main():
    # 创建ArgumentParser对象
    parser = argparse.ArgumentParser(description='提取文件中的数据并转换为数据框，然后将数据框写入输出文件。')
    
    # 添加输入文件参数
    parser.add_argument('-i','--input_file', type=str, help='输入文件路径', required=True)
    
    # 添加输出文件参数
    parser.add_argument('-o','--output_file', type=str, help='输出文件路径', required=True)
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用处理函数
    process_and_write(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
