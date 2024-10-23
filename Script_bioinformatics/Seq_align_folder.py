#!/usr/bin/env python3
import argparse
import os
import subprocess
from glob import glob

def run_blast(query_file, output_file):
    """使用blastp与本地nr数据库进行比对"""
    cmd = [
        "blastp",
        "-query", query_file,
        "-db", "/data/mengqy/database/nr_db/nr",  # 确保这里是你nr数据库的实际路径
        "-out", output_file,
        "-evalue", "e-5",  # 示例e-value阈值，可以根据需要调整
        "-max_target_seqs", "5",  # 返回的最佳命中数，根据需要调整
        "-outfmt", "6"  # 输出格式，这里使用XML格式，可根据需要调整
    ]
    #  用于创建和控制子进程。
    #  它提供了一种在Python程序中调用其他外部命令、执行系统命令和与系统进程进行交互的方法。
    #  stdout: 标准输出流（默认为None，表示将输出传递给父进程）。
    #  check：如果子进程返回非零状态码，则抛出CalledProcessError异常（默认为False）。
    subprocess.run(cmd, check=True)

def run_mafft(input_files, output_file):
    """运行MAFFT多序列比对"""
    cmd = ["mafft", input_files, ">", output_file]  # 示例命令，请根据实际需求调整
    subprocess.run(cmd, check=True)

def main():
    parser = argparse.ArgumentParser(description="Compare sequence files in a folder using BLAST or MAFFT.")
    parser.add_argument("-f", "--folder", help="Path to the folder containing sequence files.",required=True)
    parser.add_argument("-p", "--program", choices=["blast", "mafft"], default="blast", help="Choose the alignment program. Defaults to blast.",required=True)
    parser.add_argument("-o", "--output_dir", default="./outputs", help="Directory to save the output files. Defaults to './outputs'.",required=True)
    args = parser.parse_args()

    if not os.path.isdir(args.folder):
        parser.error(f"The folder '{args.folder}' does not exist.")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)  # 如果输出目录不存在，则创建它

    # 获取文件夹中所有序列文件
    sequence_files = glob(os.path.join(args.folder, "*.fasta"))  # 假设文件是FASTA格式
    if not sequence_files:
        print("No sequence files found in the specified folder.")
        return

    # 构建输出文件名
    base_name = os.path.basename(args.folder)
    output_basename = f"{base_name}_{args.program}_result.fasta"
    output_file = os.path.join(args.output_dir, output_basename)

    try:
        if args.program == "blast":
            run_blast(sequence_files, output_file)
            print(f"BLAST results saved to {output_file}")
        elif args.program == "mafft":
            run_mafft(sequence_files, output_file)
            print(f"MAFFT results saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running {args.program}: {e}")

    print('\nRunning successfully!')
    print(f"\nProcessing complete. Results saved to {args.output_dir}\n")
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()