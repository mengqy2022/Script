#!/bin/bash
# SCENIC (Single-Cell rEgulatory Network Inference and Clustering) 分析流程
# 日期: $(date +%Y-%m-%d)

# --------------------------
# 参数解析和初始化
# --------------------------

# 帮助文档函数
usage() {
    echo ""
    echo "──────────────────────────────────────────────────────────────────────────"
    echo " SCENIC (Single-Cell rEgulatory Network Inference and Clustering) 分析流程"
    echo "──────────────────────────────────────────────────────────────────────────"
    echo ""
    echo "用法:"
    echo "  bash run_pyscenic.sh [选项]"
    echo ""
    echo "选项:"
    echo "  -m <0.0-1.0>   设置内存使用阈值 (默认: 0.8，即80%内存使用时终止worker)"
    echo "  -d <目录路径>   设置工作目录路径 (包含输入文件)"
    echo "  -o <目录路径>   设置输出目录路径 (默认: 当前目录)"
    echo "  -h              显示此帮助信息"
    echo ""
    echo "示例:"
    echo "  bash run_pyscenic.sh -m 0.85 -d /data/input -o /data/output"
    echo ""
    echo "输入文件要求:"
    echo "  工作目录(-d)必须包含以下文件:"
    echo "  - hs_hgnc_tfs.txt"
    echo "  - hg38__refseq-r80__*.mc9nr.feather"
    echo "  - motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
    echo "  - scobj.loom"
    echo ""
    echo "输出文件:"
    echo "  - adj.sample.tsv (基因调控网络)"
    echo "  - reg.csv (调控模块)"
    echo "  - out_SCENIC_<日期>.loom (最终结果)"
    echo "  - scenic_analysis_<日期>.log (日志文件)"
    echo "──────────────────────────────────────────────────────────────────────────"
    exit 0
}

# 如果没有提供任何参数，显示帮助信息并退出
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# 默认参数
MEMORY_LIMIT=0.8
WORK_DIR="/data/nas1/mengqingyao_OD/project/Program-179/18_scRNA_SCENIC/cisTarget_databases"
OUTPUT_DIR=$(pwd)
DATE=$(date +%Y%m%d)
LOG_FILE="scenic_analysis_${DATE}.log"

# 使用getopts解析参数
while getopts ":m:d:o:h" opt; do
    case $opt in
        m)
            MEMORY_LIMIT=$OPTARG
            if ! awk -v n="$MEMORY_LIMIT" 'BEGIN{exit (n<0 || n>1)}'; then
                echo -e "\033[31m错误：内存阈值必须在0-1之间\033[0m" >&2
                echo "请使用 -h 查看帮助信息" >&2
                exit 1
            fi
            ;;
        d)
            WORK_DIR=$OPTARG
            if [ ! -d "$WORK_DIR" ]; then
                echo -e "\033[31m错误：工作目录不存在 - $WORK_DIR\033[0m" >&2
                exit 1
            fi
            ;;
        o)
            OUTPUT_DIR=$OPTARG
            if ! mkdir -p "$OUTPUT_DIR" 2>/dev/null; then
                echo -e "\033[31m错误：无法创建输出目录 - $OUTPUT_DIR\033[0m" >&2
                exit 1
            fi
            ;;
        h)
            usage
            ;;
        \?)
            echo -e "\033[31m错误：无效选项 -$OPTARG\033[0m" >&2
            echo "请使用 -h 查看帮助信息" >&2
            exit 1
            ;;
        :)
            echo -e "\033[31m错误：选项 -$OPTARG 需要参数\033[0m" >&2
            echo "请使用 -h 查看帮助信息" >&2
            exit 1
            ;;
    esac
done

# 记录开始时间
START_TIME=$(date +%s)
echo "====================================" >> "$LOG_FILE"
echo "SCENIC分析开始于: $(date)" >> "$LOG_FILE"
echo "参数设置:" >> "$LOG_FILE"
echo "内存阈值: $MEMORY_LIMIT" >> "$LOG_FILE"
echo "工作目录: $WORK_DIR" >> "$LOG_FILE"
echo "输出目录: $OUTPUT_DIR" >> "$LOG_FILE"
echo "====================================" >> "$LOG_FILE"

# 设置Dask内存阈值
export DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE=$MEMORY_LIMIT

# 内存监控函数
monitor_memory() {
    while true; do
        MEM_USAGE=$(free -m | awk '/Mem:/ {print $3/$2}')
        if (( $(echo "$MEM_USAGE > $MEMORY_LIMIT" | bc -l) )); then
            echo "警告: 内存使用超过阈值 ($(date))" >> "$LOG_FILE"
            echo "当前内存使用率: $(echo "$MEM_USAGE*100" | bc | cut -d. -f1)%" >> "$LOG_FILE"
            echo "系统将尝试终止最耗内存的进程..." >> "$LOG_FILE"
            
            # 找出最耗内存的进程
            PS_OUTPUT=$(ps -eo pid,%mem,cmd --sort=-%mem | head -n 2 | tail -n 1)
            PID=$(echo "$PS_OUTPUT" | awk '{print $1}')
            MEM=$(echo "$PS_OUTPUT" | awk '{print $2}')
            CMD=$(echo "$PS_OUTPUT" | awk '{for(i=3;i<=NF;i++) printf $i" "; print ""}')
            
            echo "终止进程: PID=$PID, 内存占用=$MEM%, 命令=$CMD" >> "$LOG_FILE"
            kill -9 "$PID" 2>/dev/null
            
            if [ $? -eq 0 ]; then
                echo "成功终止进程 $PID" >> "$LOG_FILE"
            else
                echo "无法终止进程 $PID" >> "$LOG_FILE"
            fi
        fi
        sleep 30
    done
}

# 启动内存监控(后台运行)
monitor_memory &
MONITOR_PID=$!

# 清理函数
cleanup() {
    echo "正在清理..." >> "$LOG_FILE"
    kill $MONITOR_PID 2>/dev/null
    pkill -f "dask-scheduler" 2>/dev/null
    pkill -f "dask-worker" 2>/dev/null
    exit
}

# 捕获退出信号
trap cleanup EXIT INT TERM

# --------------------------
# 文件路径设置
# --------------------------

tfs="$WORK_DIR/hs_hgnc_tfs.txt"
f_db_500bp="$WORK_DIR/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
f_db_10kb="$WORK_DIR/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
tbl="$WORK_DIR/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
input_loom="$WORK_DIR/scobj.loom"

# 清理残留的Dask进程
echo "清理残留Dask进程..." >> "$LOG_FILE"
pkill -f "dask-scheduler" 2>/dev/null
pkill -f "dask-worker" 2>/dev/null

# 检查所有输入文件是否存在
echo "检查输入文件..." >> "$LOG_FILE"
missing_files=0
for file in "$tfs" "$f_db_500bp" "$f_db_10kb" "$tbl" "$input_loom"; do
    if [ ! -f "$file" ]; then
        echo "错误: 文件不存在 - $file" >> "$LOG_FILE"
        missing_files=1
    fi
done

if [ $missing_files -eq 1 ]; then
    echo "错误: 必需的输入文件缺失!" >> "$LOG_FILE"
    exit 1
fi

# --------------------------
# 分析流程(其余部分保持不变，但重定向输出到日志文件)
# --------------------------

source /data/nas2/software/miniconda3/bin/activate pyscenic
echo -e "\tpyscenic环境已经激活！"

# 第一步：GRN推断
echo "开始第一步：GRN推断..." >> "$LOG_FILE"
pyscenic grn \
    --seed 777 \
    --num_workers 10 \
    --output "$OUTPUT_DIR/adj.sample.tsv" \
    --method grnboost2 \
    "$input_loom" "$tfs" >> "$LOG_FILE" 2>&1

if [ ! -s "$OUTPUT_DIR/adj.sample.tsv" ]; then
    echo "错误：adj.sample.tsv 为空或未生成!" >> "$LOG_FILE"
    exit 1
fi

# 第二步：Cistarget分析
echo "开始第二步：Cistarget分析..." >> "$LOG_FILE"
pyscenic ctx \
    "$OUTPUT_DIR/adj.sample.tsv" \
    "$f_db_500bp" "$f_db_10kb" \
    --annotations_fname "$tbl" \
    --expression_mtx_fname "$input_loom" \
    --mode "dask_multiprocessing" \
    --output "$OUTPUT_DIR/reg.csv" \
    --num_workers 15 \
    --mask_dropouts >> "$LOG_FILE" 2>&1

if [[ ! -s "$OUTPUT_DIR/reg.csv" ]]; then
    echo "致命错误：reg.csv生成失败!" >> "$LOG_FILE"
    # 诊断信息
    echo "收集诊断信息..." >> "$LOG_FILE"
    ls -lh /tmp/dask-worker-*.log >> "$LOG_FILE"
    grep -A 50 'Traceback' /tmp/dask-worker-*.log >> "$LOG_FILE"
    exit 1
fi

# 第三步：AUCell分析
echo "开始第三步：AUCell分析..." >> "$LOG_FILE"
pyscenic aucell \
    "$input_loom" \
    "$OUTPUT_DIR/reg.csv" \
    --output "$OUTPUT_DIR/out_SCENIC_${DATE}.loom" \
    --num_workers 15 \
    --seed 777 >> "$LOG_FILE" 2>&1

if [[ ! -s "$OUTPUT_DIR/out_SCENIC_${DATE}.loom" ]]; then
    echo "错误：out_SCENIC.loom生成失败!" >> "$LOG_FILE"
    exit 1
fi

# 计算运行时间
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED_TIME / 60))
ELAPSED_SEC=$((ELAPSED_TIME % 60))

echo "====================================" >> "$LOG_FILE"
echo "SCENIC分析成功完成!" >> "$LOG_FILE"
echo "总运行时间: ${ELAPSED_MIN}分${ELAPSED_SEC}秒" >> "$LOG_FILE"
echo "最终结果文件: $OUTPUT_DIR/out_SCENIC_${DATE}.loom" >> "$LOG_FILE"
echo "日志文件: $LOG_FILE" >> "$LOG_FILE"
echo "====================================" >> "$LOG_FILE"

# 将日志也输出到控制台
cat "$LOG_FILE"
