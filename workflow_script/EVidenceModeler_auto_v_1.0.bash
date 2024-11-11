#!/bin/bash
# **************************************************
# 整合各种预测证据
# Date   : 2024-11-03
# Author : 孟庆瑶
# Version: 1.0
# **************************************************

# 发生错误停止脚本
set -e

# 获取脚本名称并输出
name=$(basename $0) 
echo 
echo "   脚本名称: $name"

# 设置参数以及处理
OUTPREFIX="Euplotes"
CODE="TAA,TAG"
SPE="euplotes"
TAB="10"
NUM_THREADS=5  # 默认支持5个线程

# 设置封装
while getopts "hg:p:a:r:w:s:c:t:o:n:" opt; do
    case "$opt" in 
    h)
        echo
        echo -e "   脚本说明: [真核生物基因预测，将多个预测结果整合] \n\n             [不运行PASA，自动生成augustus和miniprot预测证据。]"
        echo -e "使用说明: bash $name -g genome.fasta -a assemblies.fasta -p PASA_result.gff3 -c TAA,TAG -r related_species.faa \n                                         -w weights.txt -s euplotes -t 10 -o out_prefix -n 5"
        echo -e "\t-g: 输入预测物种基因组 [.fasta]"
        echo -e "\t-p: PASA整合的结果 [.gff3]"
        echo -e "\t-a: PASA整合的结果 [.fasta]"
        echo -e "\t-r: 近缘物种蛋白序列 [.faa]"
        echo -e "\t-w: 权重文件 [.txt]"
        echo -e "\t-s: 物种 默认:[euplotes]"
        echo -e "\t-c: 终止密码子类型 默认:[TAA,TAG]"
        echo -e "\t-t: 翻译密码表 默认:[10]"
        echo -e "\t-o: 输出文件名称 默认:[Euplotes]"
        echo -e "\t-n: 并行计算线程数 默认:[5]"
        echo -e "\t-v: 显示版本信息；"
        echo
        echo "`date '+Date: %D %T'`"
        exit 0
        ;;
    a) PAS=$OPTARG ;;
    g) GEN=$OPTARG ;;
    p) PASA=$OPTARG ;;
    r) PROT=$OPTARG ;;
    w) WEI=$OPTARG ;;
    s) SPE=$OPTARG ;;
    c) CODE=$OPTARG ;;
    t) TAB=$OPTARG ;;
    o) OUTPREFIX=$OPTARG ;;
    n) NUM_THREADS=$OPTARG ;;
    v) echo -e "\n版本信息: v1.0\n"; exit 0 ;;
    *) echo "未知参数: $opt"; exit 1 ;;
    esac
done

# 检查必要参数
if [ -z $PROT ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $GEN ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $PASA ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $PAS ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    exit 1
fi

if [ -z $WEI ]
then
    echo -e "       [请输入-h，查看帮助文档！]       "
    echo -e "[除了具有默认的参数，其余参数都必需设置]"
    echo
    echo -e "  你真是个天才  "
    exit 1
fi

# 获取当前路径
path=$(pwd)

# 创建miniprot目录并运行miniprot
echo -e "\n [运行miniprot] \n"
mkdir -p "$path/miniprot"
source /home/mengqingyao/miniconda3/bin/activate miniprot

miniprot -t "$NUM_THREADS" --gff "$GEN" "$PROT" > "$path/miniprot/miniprot_$OUTPREFIX.gff"
grep -v "#" "$path/miniprot/miniprot_$OUTPREFIX.gff" > "$path/miniprot/miniprot_$OUTPREFIX\_mod.gff"
python3 /home/mengqingyao/miniconda3/envs/evidencemodeler/bin/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py \
    "$path/miniprot/miniprot_$OUTPREFIX\_mod.gff" > "$path/miniprot/miniprot_$OUTPREFIX\_mod_evm.gff3"

conda deactivate
echo -e "\n [蛋白证据准备完成] \n"

# 创建augustus目录并运行Augusuts
echo -e "\n [运行Augusuts] \n"
mkdir -p "$path/augusuts"
source /home/mengqingyao/miniconda3/bin/activate blat

blat -noHead "$GEN" "$PAS" "$path/augusuts/Augustus_est.psl"

data/wangruanlin/software/augustus-3.0.2/scripts/filterPSL.pl --best "$path/augusuts/Augustus_est.psl" > "$path/augusuts/Augustus_est.f.psl"
data/wangruanlin/software/augustus-3.0.2/scripts/blat2hints.pl --nomult --in="$path/augusuts/Augustus_est.f.psl" --out="$path/augusuts/Augustus_hints.est.gff"

data/wangruanlin/software/augustus-3.0.2/bin/augustus --gff3=on --species="$SPE" --protein=on --codingseq=on \
    --outfile="$path/augusuts/augustus_$OUTPREFIX.gff3" "$GEN" --translation_table="$TAB" --hintsfile="$path/augusuts/Augustus_hints.est.gff" \
    --extrinsicCfgFile=/data/wangruanlin/software/augustus-3.0.2/config/extrinsic/extrinsic.RM.E.W.cfg

awk '$3=="gene" || $3=="CDS" || $3=="transcript" {print}' "$path/augusuts/augustus_$OUTPREFIX.gff3" > "$path/augusuts/augustus_$OUTPREFIX\_mod.gff3"
python3 /home/mengqingyao/miniconda3/envs/evidencemodeler/bin/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl "$path/augusuts/augustus_$OUTPREFIX\_mod.gff3" \
    > "$path/augusuts/augustus_$OUTPREFIX\_mod_evm.gff3"

echo -e "\n [从头预测结束] \n"

# 创建EVidenceModeler目录并运行EVidenceModeler
echo -e "\n [运行EVidenceModeler] \n"
mkdir -p "$path/EVidenceModeler" && cd "$path/EVidenceModeler"
source /home/mengqingyao/miniconda3/bin/activate evidencemodeler

EVidenceModeler --sample_id "$OUTPREFIX" \
                   --genome "$path/$GEN" \
                   --weights "$path/$WEI" \
                   --gene_predictions "$path/augusuts/augustus_$OUTPREFIX\_mod_evm.gff3" \
                   --transcript_alignments "$path/$PASA" \
                   --protein_alignments "$path/miniprot/miniprot_$OUTPREFIX\_mod_evm.gff3" \
                   --segmentSize 100000 \
                   --overlapSize 10000  \
                   --stop_codons "$CODE" \
                   --min_intron_length 15 \
                   --CPU "$NUM_THREADS"

conda deactivate
cd ../
echo -e "\n [整合结束] \n"

# 获取蛋白序列
echo -e "\n [获得蛋白序列中......] \n"
source /home/mengqingyao/miniconda3/bin/activate gffread
gffread "$path/EVidenceModeler/$OUTPREFIX.EVM.gff3" -g "$GEN" -y "$path/EVidenceModeler/$OUTPREFIX.EVM.gff3.faa"
conda deactivate

python3 /home/mengqingyao/biosoftware/Script_mqy/stop_codon_replace.py "$path/EVidenceModeler/$OUTPREFIX.EVM.gff3.faa" > "$OUTPREFIX.faa"

echo "   运行结束！"
echo "   输出结果文件为: $OUTPREFIX.faa"
echo "   感谢使用本脚本！"
echo "   版本信息: v1.0"
echo "   日期: `date '+%Y-%m-%d %H:%M:%S'`"
echo "   作者: 孟庆瑶"
echo "   如有疑问请联系: <15877464851@163.com>"
