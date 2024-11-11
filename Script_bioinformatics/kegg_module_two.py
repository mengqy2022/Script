#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Email     : 15877464851@163.com
# @Time      : 2024/10/24

import os
import argparse
import pandas as pd

# 创建一个解析器
parser = argparse.ArgumentParser(description='处理指定文件夹下的所有 .ko 文件',
                                epilog= '\t更详细的信息请访问: https://mengqy2022.github.io/comparative%20genomics/Comparative-genomics/\n')

# 添加一个参数，指定文件夹路径
parser.add_argument('-f', '--folder_path', type=str, help='要读取的文件夹路径', required=True)

# 添加一个可选参数，设置输出文件的名称
parser.add_argument('-o', '--output_name', type=str, default='output.ko', help='自定义输出文件的名称 (默认: output.ko)')

# 解析命令行参数
args = parser.parse_args()

# 遍历文件夹及其子文件夹
for root, dirs, files in os.walk(args.folder_path):
    for file in files:
        # 检查文件后缀是否为 .ko
        if file.endswith('.ko'):
            # 构建文件的完整路径
            file_path = os.path.join(root, file)
            print(f"处理文件: {file_path}")

            # 读取 .ko 文件中的数据
            # 假设 .ko 文件是以制表符分隔的
            try:
                # 读取数据
                df = pd.read_csv(file_path, sep='\t', header=None, names=['pro_ids', 'ko_ids'])
                
                # 去除有一列为空的数据
                df = df.dropna()
                
                # 设置输出文件的名称，保持在同一文件夹中
                output_file_name = os.path.splitext(file)[0] + '_' + args.output_name  # 例如: 原文件名_output.ko
                output_file_path = os.path.join(root, output_file_name)

                # 将处理后的数据保存到新的文件中，使用制表符分隔
                df.to_csv(output_file_path, index=False, header=False, sep='\t')

                print(f"已保存处理后的文件到: {output_file_path}")

            except Exception as e:
                print(f"处理 {file_path} 时发生错误: {e}")
