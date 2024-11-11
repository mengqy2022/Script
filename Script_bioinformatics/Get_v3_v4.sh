#!/bin/bash
# Author: mengqingyu
# Date: 2024-11-09
# Function: Extract sequences from fasta files based on the presence of two target sequences.
# File: Get_v3_v4.sh
# Usage: ./Get_v3_v4.sh -d input_directory -o output_directory


# 用法提示
usage() {
    echo "Extract sequences from fasta files based on the presence of two target sequences."
    echo "Usage: $0 -d input_directory -o output_directory"
    exit 1
}

# 初始化参数
input_dir=""
output_dir=""

# 解析参数
while getopts ":d:o:" opt; do
    case ${opt} in
        d)
            input_dir="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# 确保输入和输出目录不为空
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    echo "Input directory and output directory are required." >&2
    usage
fi

# 检查 seqkit 是否存在
if ! command -v seqkit &> /dev/null; then
    echo "seqkit could not be found, please install it before running this script." >&2
    exit 1
fi

# 定义目标序列
target1="ACTCCTACGGGAGGCAGCAG"
target2="GGACTACHVGGGTWTCTAAT"

# 处理指定目录下的所有 .fasta 文件
for file in "$input_dir"/*.fasta; do
    # 提取文件名
    base=$(basename "$file" .fasta)

    # 查找第一个序列的位置，取所有匹配的start位置并排序
    start_output=$(seqkit locate -d -p "$target1" "$file")
    starts=($(echo "$start_output" | awk 'NR>1 {print $5}' | sort -n))
    echo "Start positions for $file: ${starts[@]}"
    echo "Full start output:"
    echo "$start_output"

    # 查找第二个序列的位置，取所有匹配的end位置并排序
    end_output=$(seqkit locate -d -p "$target2" "$file")
    ends=($(echo "$end_output" | awk 'NR>1 {print $6}' | sort -n))
    echo "End positions for $file: ${ends[@]}"
    echo "Full end output:"
    echo "$end_output"

    # 确保start和end位置有效，都是正值
    if [[ ${#starts[@]} -gt 0 && ${#ends[@]} -gt 0 ]]; then
        for ((i=0; i<${#starts[@]}; i++)); do
            start=${starts[$i]}
            end=${ends[$i]}
            # 判断start的位置和end的位置
            if [[ $start -lt $end ]]; then
                seqkit subseq -r $start:$end "$file" > "${output_dir}/${base}_${start}_${end}.fa"
            else
                seqkit subseq -r $((end-19)):$((start+19)) "$file" > "${output_dir}/${base}_${start}_${end}.fa"
            fi
        done
    else
        echo "No valid regions found in $file"
    fi
done

# 删除掉中间产物文件fai文件
rm -rf *.fai
