#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqingyao
# @Time    : 2024/12/18
# @version : 1.0
# @File    : sbatch_add.py

import argparse
import os

def convert_days_to_time_format(days):

    total_hours = days * 24
    hours = total_hours
    minutes = 0
    seconds = 0
    return f"{hours:02}:{minutes:02}:{seconds:02}"

def insert_header(input_file, output_file, threads, partition, days):

    time_limit = convert_days_to_time_format(days)

    header_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={os.path.splitext(os.path.basename(input_file))[0]}",
        f"#SBATCH --output={os.path.splitext(os.path.basename(input_file))[0]}.log",
        f"#SBATCH --error={os.path.splitext(os.path.basename(input_file))[0]}.error",
        f"#SBATCH --nodes={threads}",
        f"#SBATCH --ntasks={threads}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --time={time_limit}",
        ""
    ]

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 写入头信息
        for line in header_lines:
            outfile.write(line + '\n')

        non_empty_lines = []
        lines = infile.readlines()

        for line in lines:
            if line.strip():
                non_empty_lines.append(line.strip())

        if threads >= 2:
            for i, line in enumerate(non_empty_lines):
                outfile.write(line + ' &\n')
                if (i + 1) % threads == 0:
                    outfile.write("wait\n")
        else:
            for line in non_empty_lines:
                outfile.write(line + ' &\n')

def main():
    parser = argparse.ArgumentParser(description='在文件的第一行插入头信息',epilog=print('\n\tVersion : 1.0\n\texampe: python sbatch_add.py -i input.sh -t 2 -p apollo -d 7\n\tsbatch input_sbatch.sh\n'))
    parser.add_argument('-i', '--input_file', type=str, help='输入文件路径', required=True)
    parser.add_argument('-t', '--threads', type=int, default=1, help='节点数（默认为1）')
    parser.add_argument('-p', '--partition', type=str, default='apollo', help='分区名称（默认为apollo）')
    parser.add_argument('-d', '--days', type=int, default=7, help='运行时间（天数，默认为7天）')
    args = parser.parse_args()
    output_file = f"{os.path.splitext(args.input_file)[0]}_sbatch.sh"
    insert_header(args.input_file, output_file, args.threads, args.partition, args.days)
    print('\n\t 脚本执行完毕，输出文件为：', output_file)
    print('\n\t>>> mengqingyao <<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()
