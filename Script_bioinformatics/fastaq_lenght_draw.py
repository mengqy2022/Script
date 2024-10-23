#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
对fastq文件中的序列进行处理
1.获取序列的id和序列信息
2.统计每个id对应的序列的长度
3.对序列长度进行统计
"""
import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt


def read_fastq_seq(file_name):
    with open(file_name, 'r') as f:
        raw_index = -1
        for line in f:
            raw_index +=1
            if raw_index % 4 == 0:
                seq_id.append(line.rstrip())
            elif raw_index % 4 == 1:
                seq.append(line.rstrip())
            else:
                continue
    seq_dic = dict(zip(seq_id, seq))  # 将序列编号和序列信息两个列表合并成字典
    for v in seq_dic.values():  # 统计序列的长度
        seq_len.append(len(v))
    return seq_id, seq, seq_len, seq_dic


def count_reads_rate(seq, seq_len):
    # 统计总的reads数
    total_reads = len(seq)
    # print(total_reads)
    # 统计每个长度的reads数量，并计算其比值
    # 新建一个字典用于存储长度对应的reads信息

    for read_seq in seq_len:
        read_len_dict.setdefault(read_seq,0)
        read_len_dict[read_seq] += 1
    # print(read_len_dict)
    return read_len_dict


def write_count(seq_id, seq, seq_len, read_len_dict):
    # pandas 写入序列编号、序列信息、序列长度
    data1 = pd.DataFrame({"Seq_ID": seq_id})
    data2 = pd.DataFrame({"Seq_Info": seq})
    data3 = pd.DataFrame({"Seq_Len": seq_len})
    # 统计不同长度序列的条数
    k_len = []
    v_count = []
    for k_len, v_count in read_len_dict.items():
        k_len.append()
        v_count.append()
        return k_len, v_count
    data4 = pd.DataFrame({ "Read_Len":k_len})
    data5 = pd.DataFrame({ "Read_Len":v_count})

    # 将序列编号，序列信息，以及序列长度写入*.xlsx文件中
    writer = pd.ExcelWriter(abs_path + '\\' + "test3.xlsx")
    data1.to_excel(writer, sheet_name="data", startcol=0, index=False)
    data2.to_excel(writer, sheet_name="data", startcol=1, index=False)
    data3.to_excel(writer, sheet_name="data", startcol=2, index=False)
    # 将序列统计后的数据写入 reads_len_count中
    data4.to_excel(writer,sheet_name="reads_len_count", startcol=0, index=False)
    data5.to_excel(writer, sheet_name="reads_len_count", startcol=1, index=False)
    writer.save()  # 数据保存为excel文件


def count_bar(seq_len):
    """
    根据上一步获得的序列长度信息，对其进行sort/uniq，matplotlib 处理并绘制直方图
    首先对数据进行排序统计对相同长度进行计数
    数据清洗后进行画bar图
    提取上一步获得第三列长度数据进行清洗并统计每个长度的个数
    """
    len_count = Counter(seq_len)
    # matplotlib绘图
    x = []
    y = []
    for k, v in len_count.items():
        x.append(k)
        y.append(v)
    # print(x, y)
    plt.bar(x, y)
    plt.xlabel("Sequence Length")
    plt.ylabel("Sequence Numbers")
    # plt.savefig()
    plt.show()


if __name__ == "__main__":
    abs_path = os.getcwd()  # 获取当前目录路径
    # print(abs_path)
    file_name = abs_path + '\\' + 'HMP_MOCK_SRR2726667_8.25M.1.trim.100.fastq.gz'  # 获取当前目录下的文件信息
    seq_id = []  # 新建列表存储fasta文件中序列编号信息
    seq = []  # 新建列表存储对应fasta文件中序列信息
    seq_len = []  # 新建列表存储对应序列的长度信息
    read_len_dict = {}
    read_fastq_seq(file_name)
    write_count(seq_id, seq, seq_len,read_len_dict)
    count_reads_rate(seq, seq_len)
    count_bar(seq_len)

