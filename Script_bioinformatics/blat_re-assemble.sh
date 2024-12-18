#!/bin/bash

# 使用说明
usage() {
    echo "Usage: $0 [input_directory] [contig_file] [file_extension]"
    echo "  [input_directory]: Directory containing the files for processing."
    echo "  [contig_file]: Path to the Contig file used for BLAT comparison."
    echo "  [file_extension]: File extension to process (e.g., .fq.gz or .fp.fq.gz)."
    echo " >>> Example: $0 /mnt/e/Euplotes/Euplotes_amieti/00_艾美核geneo ../Contig6391.txt .fq.gz <<<"
    exit 1
}

# 检查参数数量
if [ "$#" -ne 3 ]; then
    usage
fi

# 定义输入目录、contig 文件和文件扩展名
INPUT_DIR="$1"
CONTIG_FILE="$2"
FILE_EXT="$3"

# 确保输入目录存在
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory '$INPUT_DIR' does not exist."
    exit 1
fi

# 记录开始时间
start_time=$(date +%s)

# 遍历输入目录下的所有指定扩展名文件
for fq_file in "$INPUT_DIR"/*"$FILE_EXT"; do
    # 检查文件是否存在
    if [ ! -f "$fq_file" ]; then
        echo "No files with extension '$FILE_EXT' found in '$INPUT_DIR'."
        exit 1
    fi

    # 提取文件名基础部分（去掉扩展名）
    base_name=$(basename "$fq_file" "$FILE_EXT")

    # 检查是否存在同名的 .fp.* 文件
    if [ -f "$INPUT_DIR/${base_name}.fp${FILE_EXT}" ]; then
        echo "Skipping '${base_name}.fp${FILE_EXT}' as '${base_name}${FILE_EXT}' is present."
        continue
    fi

    # 输出路径定义
    fa_file="${base_name}.fa"
    blat_out="${base_name}.blat"
    mod_blat_out="${base_name}_mod.blat"
    mod_ids="${base_name}_mod.ids"
    mod_fq="${base_name}_mod${FILE_EXT}"

    # 将fastq转为fasta格式
    seqkit fq2fa "$fq_file" > "$fa_file"

    # 用 BLAT 比对
    blat "$CONTIG_FILE" "$fa_file" -minIdentity=99 -out=blast8 "$blat_out"

    # 额外筛选：保留第四列大于等于50和第五列小于等于2的数据，同时去重
    awk 'BEGIN {OFS="\t"} {if (!tophit[$1] && $4 >= 50 && $5 <= 2) {tophit[$1]; print}}' "$blat_out" > "$mod_blat_out"

    # 提取 ID
    cut -f 1 "$mod_blat_out" > "$mod_ids"

    # 用 seqkit grep 提取相关序列
    seqkit grep -f "$mod_ids" "$fq_file" -o "$mod_fq"

done

# 记录结束时间
end_time=$(date +%s)

# 计算并输出总运行时间
runtime=$((end_time - start_time))
echo "Total execution time: $runtime seconds"
