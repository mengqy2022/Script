#!/bin/bash

# 帮助文档
usage() {
    echo "Usage: $0 --option-file IDS_FILE --download-path DOWNLOAD_PATH"
    echo "  --option-file IDS_FILE: 包含SRR IDs的文件路径"
    echo "  --download-path DOWNLOAD_PATH: 下载的SRA文件保存的路径"
    exit 1
}

# 参数解析
if [[ $# -ne 4 ]]; then
    usage
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        --option-file)
            IDS_FILE="$2"
            shift 2
            ;;
        --download-path)
            DOWNLOAD_PATH="$2"
            shift 2
            ;;
        *)
            usage
            ;;
    esac
done

# 检查参数
if [[ ! -f "$IDS_FILE" ]]; then
    echo "Error: 文件 $IDS_FILE 不存在"
    exit 1
fi

# 使用prefetch下载SRA文件
prefetch --option-file "$IDS_FILE" -O "$DOWNLOAD_PATH"

# 遍历下载的SRA文件并使用fasterq-dump转换为fastq
for SRA_FILE in "$DOWNLOAD_PATH"/*.sra; do
    BASENAME=$(basename "$SRA_FILE" .sra)
    cd "$DOWNLOAD_PATH/$BASENAME"
    time fasterq-dump "$SRA_FILE" --split-3 -e 20 -o "$BASENAME" -p
done

echo "所有SRA文件已成功转换为fastq文件"
