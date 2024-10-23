#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2024/09/29
# @File    : Phylogenetic_reads_one

import os
import pandas as pd
import argparse

class NameExtractor:
    def __init__(self, folder, separator, output_file):
        self.folder = folder
        self.separator = separator
        self.output_file = output_file
        self.names = []
        self.original_names = []

    def extract_names_from_folder(self):
        """从指定文件夹提取文件名并保存到数据框中"""
        try:
            file_names = os.listdir(self.folder)
            for file_name in file_names:
                if os.path.isfile(os.path.join(self.folder, file_name)):
                    
                    # 去掉文件扩展名
                    remove_file_name = os.path.splitext(file_name)[0]
                    # 保存原文件名
                    self.original_names.append(remove_file_name)

                    # 处理文件名
                    first_name_part = file_name.split(self.separator)[0]
                    if '.' in first_name_part:
                        first_name_part = first_name_part.split('.')[0]
                    if '-' in first_name_part:
                        first_name_part = first_name_part.split('-')[0]
                    if ' ' in first_name_part:
                        first_name_part = first_name_part.split(' ')[0]
                    self.names.append(first_name_part)

            # 创建一个DataFrame并保存为用户指定的CSV文件
            df = pd.DataFrame({
                'rawnames': self.original_names,
                'renames': self.names
            })
            df.to_csv(self.output_file, index=False, encoding='utf-8', sep='\t')
            print(f"文件名已保存到 '{self.output_file}' 中。")

        except FileNotFoundError:
            print(f"文件夹 '{self.folder}' 不存在。")

def main():
    parser = argparse.ArgumentParser(description="从指定文件夹提取文件名并处理")
    parser.add_argument('-f', '--folder', type=str, help='要提取文件名的文件夹路径', required=True)
    parser.add_argument('-s', '--separator', help='指定分隔符，默认为"_"', type=str, required=False, default='_')
    parser.add_argument('-o', '--output_file', type=str, help='输出文件名，包括文件扩展名（如 .csv）', required=True)
    args = parser.parse_args()

    extractor = NameExtractor(args.folder, args.separator, args.output_file)
    extractor.extract_names_from_folder()

if __name__ == "__main__":
    main()
