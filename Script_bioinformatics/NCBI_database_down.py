#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Email     : 15877464851@163.com
# @Time      : 2024/11/20

import requests
import json
import subprocess
import os
import glob
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm 

# 定义两个不同的数据库URL
NR_METADATA_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/nr-prot-metadata.json"
NT_METADATA_URL = "https://ftp.ncbi.nlm.nih.gov/blast/db/nt-nucl-metadata.json"

ADDITIONAL_URLS = {
    "nr_fasta": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
    "nt_fasta": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz",
    "swissprot": "https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz",
}

DATABASE_URLS = {
    "nr": NR_METADATA_URL,
    "nt": NT_METADATA_URL,
}

def download_database(metadata_url, output_dir):
    response = requests.get(metadata_url)

    if response.status_code == 200:
        data = response.json()
        files = data.get("files", [])

        # 设置输出目录
        os.makedirs(output_dir, exist_ok=True)
        os.chdir(output_dir)  # 切换到输出目录

        # 将文件链接写入文本文件
        with open("file_links.txt", "w") as f:
            for file_url in files:
                f.write(file_url + "\n")
        print(f"{metadata_url} 的文件链接已成功写入 file_links.txt。")

        # 使用wget下载文件，添加重试参数
        try:
            print(f"开始下载 {metadata_url} 的文件...")
            process = subprocess.Popen(
                [
                    "wget",
                    "-c",
                    "--retry-connrefused",  # 重试被拒绝的连接
                    "--tries", "0",          # 设置为无限次重试
                    "-t", "0",               # 指定重试连接次数
                    "-i", "file_links.txt"
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True
            )

            total_files = len(files)
            downloaded_files = 0

            # 使用tqdm显示进度条
            with tqdm(total=total_files, desc="下载进度", unit="file") as pbar:
                while True:
                    output = process.stdout.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        sys.stdout.write(output)
                        sys.stdout.flush()
                        # 检测文件下载完成
                        if "100%" in output:
                            downloaded_files += 1
                            pbar.update(1)  # 更新进度条

            process.wait()

            # 检查下载的文件数量是否与预期相等
            if downloaded_files == total_files:
                print("所有文件下载完成。")
                return True
            else:
                print(f"下载完成的文件数量与预期不符：{downloaded_files}/{total_files}。")
                return False

        except Exception as e:
            print(f"下载过程发生异常: {e}")
            return False
    else:
        print(f"请求失败，状态码：{response.status_code}")
        return False

def download_additional_file(file_url, output_dir):
    # 设置输出目录
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)  # 切换到输出目录

    print(f"开始下载 {file_url}...")
    process = subprocess.Popen(
        [
            "wget",
            "-c",
            "--retry-connrefused",  # 重试被拒绝的连接
            "--tries", "0",          # 设置为无限次重试
            file_url
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    # 使用tqdm显示下载进度条
    with tqdm(total=100, desc="额外文件下载进度", unit="%") as pbar:
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                sys.stdout.write(output)
                sys.stdout.flush()
                # 检测文件下载完成
                if "100%" in output:
                    pbar.update(100)  # 更新进度条到100%
    
    process.wait()

    if process.returncode == 0:
        print(f"{file_url} 下载完成。")
        return True
    else:
        print(f"{file_url} 下载过程中出现错误。")
        return False

def extract_files(output_dir):
    os.chdir(output_dir)  # 切换到输出目录
    downloaded_files = glob.glob("*.tar.gz")
    
    if downloaded_files:
        print("开始解压下载的文件...")
        for tar_file in downloaded_files:
            subprocess.run(["tar", "-xzvf", tar_file])
            print(f"{tar_file} 解压完成。")

        # 删除原始的.tar.gz文件
        for tar_file in downloaded_files:
            os.remove(tar_file)
            print(f"已删除文件: {tar_file}")
    else:
        print("没有找到需要解压的文件。")

def main():
    parser = argparse.ArgumentParser(description="Download and extract NCBI BLAST databases.",
                                     epilog="更新数据库信息，如果需要构建diamond数据库，请下载fasta序列文件。")
    parser.add_argument(
        "--blast_db",
        required=False,
        choices=DATABASE_URLS.keys(),
        help="选择要下载的数据库，支持：'nr' 或 'nt'。"
    )
    parser.add_argument(
        "--fasta",
        required=False,
        choices=ADDITIONAL_URLS.keys(),
        help="选择要下载的额外文件，支持：'nr_fasta', 'nt_fasta', 'swissprot'。"
    )

    parser.add_argument(
        "-od","--output_dir",
        type=str,
        required=True,
        help="输出数据库文件的目录."
    )
    parser.add_argument(
        "-nw", "--num_wget",
        type=int,
        default=3,
        help="设置同时运行的 wget 命令数量，默认 3."
    )
    
    args = parser.parse_args()

    # 下载选择的数据库
    if args.blast_db:
        with ThreadPoolExecutor(max_workers=args.num_wget) as executor:
            print(f"调度下载数据库: {args.blast_db}")
            executor.submit(download_database, DATABASE_URLS[args.blast_db], args.output_dir)

    # 下载额外的文件
    if args.fasta:
        print(f"调度下载额外的文件: {args.fasta}")
        download_additional_file(ADDITIONAL_URLS[args.fasta], args.output_dir)

    # 所有下载完成后解压和删除原始文件
    extract_files(args.output_dir)

if __name__ == "__main__":
    main()