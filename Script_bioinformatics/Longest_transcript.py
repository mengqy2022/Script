#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse

def get_info(arg1, arg2):
    # 首先按行读取gff3文件
    gene_info = {}
    mRNA_length = {}
    mRNA_compare_length = {}
    mRNA_compare_length_2 = {}
    CDS_info = {}
    CDS_info_compare = {}

    # 读取gff3文件
    with open(arg1, "r") as fr:
        for line in fr:
            if line.startswith("#"):
                continue
            Info = line.strip().split("\t")
            if Info[2] == "gene":
                # print(Info[8])
                gene_info[Info[8].split(";")[1].replace("Name=", "")] = [Info[0], Info[3], Info[4],
                                                                         int(Info[4]) - int(Info[3])]
            if Info[2] == "mRNA":
                # state length
                mRNA_length[Info[8].split(";")[1].replace("Name=", "")] = int(Info[4]) - int(Info[3])
            if Info[2] == "CDS":
                CDS_id = Info[8].split(";")[1].replace("Parent=", "")
                CDS_id_length = int(Info[4]) - int(Info[3])
                if CDS_id not in CDS_info:
                    CDS_info[CDS_id] = CDS_id_length
                else:
                    CDS_info[CDS_id] += CDS_id_length

    # compare longest mRNA
    for k, v in mRNA_length.items():
        if k.split(".")[0] not in mRNA_compare_length:
            mRNA_compare_length[k.split(".")[0]] = [(k, v)]  #  VALUES : ID + LENGTH
        else:
            mRNA_compare_length[k.split(".")[0]].append((k, v))

    # print(len(mRNA_compare_length))
    # print(mRNA_compare_length)


    # compare longest CDS
    for k, v in mRNA_length.items():
        if k.split(".")[0] not in mRNA_compare_length_2:
            mRNA_compare_length_2[k.split(".")[0]] = v
            mRNA_compare_length_2[k] = v
        else:
            if v > mRNA_compare_length_2[k.split(".")[0]]:
                mRNA_compare_length_2[k] = v
                del mRNA_compare_length_2[k.split(".")[0]]
    mRNA_compare_length_3 = {k: v for k, v in mRNA_compare_length_2.items() if "." in k}

    # compare longest CDS
    for k, v in CDS_info.items():
        if k.split(".")[0] not in CDS_info_compare:
            CDS_info_compare[k.split(".")[0]] = v
        else:
            if v > CDS_info_compare[k.split(".")[0]]:
                CDS_info_compare[k.split(".")[0]] = v
            else:
                continue

    # 读取序列文件 protein.fa or CDS.fa or cdns.fa
    sequence = {}
    with open(arg2, "r") as seqr:
        for line in seqr:
            if line.startswith(">"):
                ID = line.strip().split(" ")[0].replace(">", "")
                sequence[ID] = ""
            else:
                seq = line.strip()
                sequence[ID] += seq

    # print(sequence)
    # 输出结果
    # 输出基因的物理位置信息
    Info_w = open("Gene_Info.txt", "w")
    for k, v in gene_info.items():
        Info_w.write(v[0] + "\t" + v[1] + "\t" + v[2] + "\t" + k + "\t" + str(v[3]) + "\n")
    Info_w.close()

    # 输出最长转录本序列
    Longest_trans_w = open("Longest_transcript.fa", "w")
    for key in mRNA_compare_length_2.keys():
        if key in sequence.keys():
            Longest_trans_w.write(">" + key + "\n" + sequence[key] + "\n")
    Longest_trans_w.close()

if __name__ == "__main__":
    # create ArgumentParser
    parser = argparse.ArgumentParser(description="This Program will get gene information and longest transcription sequence by GFF3 and sequence file")
    # add argument
    parser.add_argument("--gff", "-g",nargs='+', required=True, help="Input GFF3 file")
    parser.add_argument("--fasta", "-f",nargs='+', required=True, help="Input fasta sequence file")
    # parse argument
    args = parser.parse_args()

    if args.gff and args.fasta:
        get_info(args.gff, args.fasta)

    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')