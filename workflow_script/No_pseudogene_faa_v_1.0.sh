#!/bin/bash
# **************************************************
# 共生菌基因预测以及去除假基因
# Date   : 2024-11-03
# Author : 孟庆瑶
# **************************************************

# 发生错误停止脚本
set -e

# 获取脚本名称并输出
name=$(basename "$0") 
echo 
echo "   脚本名称: $name"

# 设置参数以及处理
OUTPREFIX="Bacterial_results"
ONLY_PROKKA=false
MAX_JOBS=3  # 默认最大并行数量为 3

while getopts "hi:g:o:j:vp" opt; do
    case "$opt" in 
        h)
            echo -e "\n脚本说明: 基因预测、假基因预测和获得去除假基因的CDs序列，(*￣︶￣)，东西不多，一字一字看！"
            echo -e "使用说明: bash $name -i genome_file_name -g genome_fasta -o out_prefix -j max_jobs -p"
            echo -e "\t-i: 输入文本文件、包含基因文件名称、不带后缀、单列文件，一列为一个基因组；"
            echo -e "\t-g: 基因文件夹名称，文件夹中可以包含一个或者多个基因组；"
            echo -e "\t-o: 文件前缀, 默认:Bacterial_results；"
            echo -e "\t-j: 最大并行任务数量，默认:3个任务，9个线程；"
            echo -e "\t-v: 显示版本信息；"
            echo -e "\t-p: 仅进行 Prokka 预测；"
            echo "如有疑问请联系: <15877464851@163.com>\n"
            echo
            echo "`date '+Date: %D %T'`"
            exit 0
            ;;
        i)
            FILENAME="$OPTARG"
            ;;
        g)
            FOLDERNAME="$OPTARG"
            ;;
        o)
            OUTPREFIX="$OPTARG"
            ;;
        j)
            MAX_JOBS="$OPTARG"  # 用户指定的最大并行任务数量
            ;;
        v)
            echo -e "\n版本信息: v1.0\n"
            exit 0
            ;;
        p)
            ONLY_PROKKA=true
            ;;
        *)
            echo "未知参数: $opt" >&2
            exit 1
            ;;
    esac
done

# 检查输入文件名
if [ -z "$FILENAME" ]; then
    echo -e "\n请输入-h，查看帮助文档！"
    exit 1
fi

# 创建输出目录
if [ -d "$OUTPREFIX" ]; then
    rm -rf "$OUTPREFIX" 
fi
mkdir "$OUTPREFIX" && cd "$OUTPREFIX" || { echo "无法创建输出目录 $OUTPREFIX"; exit 1; }

# 处理输入文件中的基因组名称
cd .. || { echo "无法返回上级目录"; exit 1; }
while IFS= read -r genome; do
    mkdir -p "$FOLDERNAME/$genome" || { echo "无法创建目录 $FOLDERNAME/$genome"; exit 1; }
    mv "$FOLDERNAME/$genome."* "$FOLDERNAME/$genome" 2>/dev/null || echo "没有文件需要移动到 $FOLDERNAME/$genome"
done < "$FILENAME"

# Prokka
source /home/mengqingyao/miniconda3/bin/activate prokka || { echo "无法激活 prokka 环境"; exit 1; }
mkdir -p prokka && cd prokka || { echo "无法进入 prokka 目录"; exit 1; }

run_prokka() {
    prokka "$FOLDERNAME/$1/$1."* --cpus "$MAX_JOBS" --force --outdir "${1}_prokka" --prefix "$1" --addmrna --compliant 200 || echo "Prokka 执行失败: $1"
}

while IFS= read -r genome; do 
    while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
    echo "Prokka $genome"
    run_prokka "$genome" & 
done < "../$FILENAME"

wait
echo "Prokka 完成！"

# 仅进行 Prokka 预测则结束
if [ "$ONLY_PROKKA" = true ]; then
    echo "只进行了 Prokka 预测，脚本结束。"
    exit 0
fi

# Pseudofinder
conda activate pseudofinder || { echo "无法激活 pseudofinder 环境"; exit 1; }
mkdir -p pseudogene && cd pseudogene || { echo "无法进入 pseudogene 目录"; exit 1; }

run_pseudofinder() {
    pseudofinder.py annotate --genome "../prokka/${1}_prokka/${1}.gb*" --outprefix "$1" \
    -di --database /data/mengqy/database/diamond_nr/diamond_makedb.dmnd --threads "$MAX_JOBS" -e 1e-15 || echo "Pseudofinder 执行失败: $1"
}

while IFS= read -r genome; do 
    while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
    echo "pseudofinder $genome"
    run_pseudofinder "$genome" & 
done < "../$FILENAME"

wait
echo "pseudofinder 完成！"

# 去除假基因
cd .. || { echo "无法返回上级目录"; exit 1; }
cd "$OUTPREFIX" || { echo "无法进入输出目录 $OUTPREFIX"; exit 1; }
mkdir -p remove_pseudogene && cd remove_pseudogene || { echo "无法进入 remove_pseudogene 目录"; exit 1; }

run_remove_pseudogene() {
    awk '/old_locus_tag=/{print substr($0, RSTART+RLENGTH)}' "../pseudogene/${1}_pseudos.gff" | tr ',' '\t' | awk '{print $1}' | sed '/^$/d' | sed 's/_1_ign//g' | sort | uniq > "${1}_pseudogene.ids" || { echo "提取假基因ID失败: $1"; return 1; }
    
    comm -23 <(grep "^>" "../prokka/${1}_prokka/${1}.ffn" | cut -d ' ' -f 1 | cut -c 2- | sort) <(sort "${1}_pseudogene.ids") | seqkit grep -f - "../prokka/${1}_prokka/${1}.ffn" -o "${1}_remove.ffn" || { echo "提取非假基因序列失败: $1"; return 1; }
    comm -23 <(grep "^>" "../prokka/${1}_prokka/${1}.faa" | cut -d ' ' -f 1 | cut -c 2- | sort) <(sort "${1}_pseudogene.ids") | seqkit grep -f - "../prokka/${1}_prokka/${1}.faa" -o "${1}_remove.faa" || { echo "提取非假基因序列失败: $1"; return 1; }
}

while IFS= read -r genome; do 
    while [ "$(jobs -r | wc -l)" -ge "$MAX_JOBS" ]; do
        sleep 1
    done
    echo "获得去除假基因的序 $genome"
    run_remove_pseudogene "$genome" & 
done < "../$FILENAME"

wait
echo "获得去除假基因的序 完成！"

# 输出结果
echo "输出结果文件为: $(pwd)"
echo "输出结果中包含:"
echo "               基因预测文件;"
echo "               假基因预测文件;"
echo "               去除假基因序列文件。"
echo ""
echo "   运行结束！"
echo "   感谢使用本脚本！"
echo "   版本信息: v1.0"
echo "   日期: $(date '+%Y-%m-%d %H:%M:%S')"
echo "   作者: 孟庆瑶"
echo "   如有疑问请联系: <15877464851@163.com>"
