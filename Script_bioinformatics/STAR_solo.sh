#!/bin/bash

# STARsolo 传参脚本

# 基本用法提示
brief_usage() {
    echo "Usage: $0 -s <sample> -c <cDNA> -v <v2|v3> -i <in-path> -o <out-path> -w <whitelist>"
    echo "Use '-h' for detailed help"
    exit 1
}

# 详细帮助信息
detailed_help() {
    echo "STARsolo Processing Script"
    echo "Usage: $0 [options]"
    echo ""
    echo "Required Options:"
    echo "  -s <sampleID>         Sample ID"
    echo "  -c <cDNA_read>        cDNA reads filename"
    echo "  -v <v2|v3>            Barcode version (v2 or v3)"
    echo "  -i <path>             Input directory path"
    echo "  -o <path>             Output directory path"
    echo "  -w <file>             Whitelist file"
    echo ""
    echo "Other Options:"
    echo "  -h                    Show this help message"
    exit 0
}

# 如果没有提供任何参数，显示简要用法
if [ $# -eq 0 ]; then
    brief_usage
fi

# 使用getopt解析参数
PARSED=$(getopt -o s:c:v:i:o:w:h --long sample:,cDNA:,version:,in-path:,out-path:,whitelist:,help -n "$0" -- "$@")

if [ $? != 0 ]; then
    brief_usage
fi

eval set -- "$PARSED"

# 初始化变量
sampleID=""
cDNA_read=""
barcode_version=""
inPath=""
outPath=""
whiteList=""

# 解析参数
while true; do
    case "$1" in
        -s|--sample)
            sampleID="$2"
            shift 2
            ;;
        -c|--cDNA)
            cDNA_read="$2"
            shift 2
            ;;
        -v|--version)
            barcode_version="$2"
            shift 2
            ;;
        -i|--in-path)
            inPath="$2"
            shift 2
            ;;
        -o|--out-path)
            outPath="$2"
            shift 2
            ;;
        -w|--whitelist)
            whiteList="$2"
            shift 2
            ;;
        -h|--help)
            detailed_help
            ;;
        --)
            shift
            break
            ;;
        *)
            brief_usage
            ;;
    esac
done

# 检查必需参数
if [[ -z "$sampleID" || -z "$cDNA_read" || -z "$barcode_version" || -z "$inPath" || -z "$outPath" || -z "$whiteList" ]]; then
    echo "Error: Missing required parameters"
    brief_usage
fi

# 剩余脚本保持不变...
## STAR config
refIndex="~/upstream/refIndex/GRCh38-2020-A"
CPU=16

## 根据barcode版本选择设置
if [ "$barcode_version" == "v2" ]; then
    barcode_setting="--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloBarcodeReadLength 0"
elif [ "$barcode_version" == "v3" ]; then
    barcode_setting="--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloBarcodeReadLength 0"
else
    echo "Error: Invalid barcode version. Please specify 'v2' or 'v3'."
    exit 1
fi

## input config
inFASTQ_cDNA="$inPath/$sampleID/$cDNA_read"

## output config
soloFeatures="Gene Velocyto"
outSAMSettings="--outSAMtype None"
outPrefix="$outPath/$sampleID"

##==== cmds ====
echo "Running STARsolo with following parameters:"
echo "  Sample ID: $sampleID"
echo "  cDNA reads: $inFASTQ_cDNA"
echo "  Barcode version: $barcode_version"
echo "  Input path: $inPath"
echo "  Output path: $outPath"
echo "  Whitelist file: $whiteList"

# 检查输入文件是否存在
if [ ! -f "$inFASTQ_cDNA" ]; then
    echo "Error: cDNA file $inFASTQ_cDNA not found!"
    exit 1
fi

# 检查白名单文件是否存在
if [ ! -f "$whiteList" ]; then
    echo "Error: whitelist file $whiteList not found!"
    exit 1
fi

# 创建输出目录
mkdir -p "$outPrefix"

STAR --genomeDir $refIndex \
--runThreadN $CPU \
$outSAMSettings \
--outFileNamePrefix $outPrefix/ \
--readFilesIn $inFASTQ_cDNA \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
$barcode_setting \
--soloCBwhitelist $whiteList \
--soloCellFilter EmptyDrops_CR \
--soloFeatures $soloFeatures

if [ $? -eq 0 ]; then
    echo "STARsolo processing completed successfully for sample $sampleID"
else
    echo "Error: STARsolo processing failed for sample $sampleID"
    exit 1
fi