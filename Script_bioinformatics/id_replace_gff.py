#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

if len(sys.argv) != 4:
    print(f"Usage: python {sys.argv[0]} <id_list_file> <gff_file> <out_file>")
    print(f"\n\tid_list——file: [ord_ID   parent_ID   new_ID]")
    sys.exit(1)

id_list_file = sys.argv[1]
gff_file = sys.argv[2]
out_file = sys.argv[3]

# 定义一个空字典，用于存储ID列表文件中的信息
id_dict = {}
# 读取ID列表文件，将其中的每行ID、Parent和新ID信息以'id_str':[parent_str, new_id_str]的形式存储在字典中
with open(id_list_file, "r") as id_file:
    for line in id_file:
        id_str, parent_str, new_id_str = line.strip().split()
        id_dict[id_str] = [parent_str, new_id_str]

# 打开gff输入文件和gff输出文件
with open(gff_file, "r") as f, open(out_file, "w") as out_file:
    # 逐行读取gff输入文件
    for line in f:
        # 将读取的行按照制表符分隔成列表
        cols = line.strip().split('\t')

        # 对于mRNA和CDS行信息，进行重新组合
        if cols[2] in ["mRNA", "CDS"]:
            # 初始化ID信息为空字符串
            id_str = ''
            # 循环遍历该行的最后一个字段，解析出ID信息
            for ele in cols[-1].split(";"):
                if ele.startswith('ID='):
                    id_str = ele.split("=")[1]
                    break
            # 如果ID信息存在于id_dict中，更新该行信息的ID和Parent信息，并将行信息转换为字符串写入输出文件中
            if id_str in id_dict:
                parent_str, new_id_str = id_dict[id_str]
                cols[-1] = f"ID={new_id_str};Parent={new_id_str}"
                out_file.write("\t".join(cols) + "\n")

        # 对于UTR行信息，进行重新组合
        elif cols[2] == "UTR":
            # 初始化ID信息为空字符串
            id_str = ''
            # 循环遍历该行的最后一个字段，解析出ID信息
            for ele in cols[-1].split(";"):
                if ele.startswith('Parent='):
                    id_str = ele.split("=")[1]
                    break
            # 如果ID信息存在于id_dict中，更新该行信息的Parent信息，并将行信息转换为字符串写入输出文件中
            if id_str in id_dict:
                parent_str, new_id_str = id_dict[id_str]
                cols[-1] = f"Parent={new_id_str}"
                out_file.write("\t".join(cols) + "\n")

        # 对于其他行信息，直接写入输出文件中
        else:
            out_file.write(line)