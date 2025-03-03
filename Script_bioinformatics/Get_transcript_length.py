#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import os
import re
from collections import defaultdict
import argparse

def get_attribute_label(line: str, label: str):
    match = re.search(rf'{label} "(.*?)"', line)
    return match.group(1) if match else None

def remove_id_version(id: str, remove: bool = False):
    return re.sub(r'\.\d+', '', id) if remove else id

def extract_transcript_info(gtf_file: str, output_file: str, remove_version: bool = True):
    try:
        transcript_lengths = defaultdict(int)
        gene_max_lengths = defaultdict(int)

        with open(gtf_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                split_line = line.split("\t")
                if split_line[2] == "exon":
                    start, end = int(split_line[3]), int(split_line[4])
                    exon_length = end - start + 1
                    transcript_id = get_attribute_label(split_line[8], "transcript_id")
                    gene_id = get_attribute_label(split_line[8], "gene_id")
                    gene_name = get_attribute_label(split_line[8], "gene_name") or gene_id

                    transcript_id = remove_id_version(transcript_id, remove_version)
                    gene_id = remove_id_version(gene_id, remove_version)

                    transcript_lengths[(gene_id, gene_name, transcript_id)] += exon_length
                    if transcript_lengths[(gene_id, gene_name, transcript_id)] > gene_max_lengths[gene_id]:
                        gene_max_lengths[gene_id] = transcript_lengths[(gene_id, gene_name, transcript_id)]

        with open(output_file, "w", encoding="utf-8") as out:
            out.write("gene_id\tgene_name\ttranscript_id\ttranscript_length\tis_longest\n")
            for (gene_id, gene_name, transcript_id), length in transcript_lengths.items():
                is_longest = "Yes" if length == gene_max_lengths[gene_id] else "No"
                out.write(f"{gene_id}\t{gene_name}\t{transcript_id}\t{length}\t{is_longest}\n")

    except Exception as e:
        print_error(f"ERROR: An error occurred while processing the file: {str(e)}")

def print_error(*args, exit_code: int = 1):
    err_infor = " ".join(args)
    err_infor = f"\n\033[31m{err_infor}\033[0m\n"
    print(err_infor, file=sys.stderr)
    sys.exit(exit_code)

def check_files(file: str, parser=None):
    if not os.path.isfile(file):
        if parser is not None:
            parser.print_help()
        print_error(f"ERROR: 未找到文件: {file}")

if __name__ == "__main__":
    desc = """从GTF文件中提取基因ID、基因名称、转录本ID和转录本长度。\n\t还标记转录本是否是该基因中最长的一个。"""
    epilog = """\t版本: 0.0.1
\t电子邮件: zhengshimao007@163.com
\t创建日期: 2025.02.24
\t更新日期: -
\t作者: Shimao Zheng
"""
    parser = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--gtf_file", required=True, help="输入GTF文件")
    parser.add_argument("-o", "--output_file", required=True, help="输出文件")
    parser.add_argument("--remove_version", action="store_true",
                        help="从基因和转录本ID中移除版本号 [默认: False]", default=False)
    args = parser.parse_args()

    check_files(file=args.gtf_file, parser=parser)
    extract_transcript_info(args.gtf_file, args.output_file, args.remove_version)
