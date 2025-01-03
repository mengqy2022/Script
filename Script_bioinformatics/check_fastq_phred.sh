#!/usr/bin/bash

## ---------- 定义脚本的使用方法 ----------------------------------
usage() {
cat << EOF
Usage:bash $0 -i <fastq.gz>
 bash $0 -i <fastq.gz> -n3000
 bash $0 -i <fastq.gz> --read_nums=3000
 ls *.fastq.gz |xargs -i bash $0 -i {}

Options:
 -i, --input Input file.
 -o, --ouput Output file.
 -n, --read_nums A positive integer. Number of reads to test [default: 2500].
 -f, --force force tu rerun this step. 
 -h, --help Print the tips.
EOF
    exit 1
}

## ---------- 脚本选项设置 ----------------------------------
#  $@：表示执行脚本传入参数的所有个数（不包括$0）
#  $#：表示执行脚本传入参数的个数（包括$0）
OPTIONS=$(getopt -o hfi:o::n:: --long help,force,input:,output::,read_nums:: -n "ERROR:$0" -- "$@")

#if [ $? != 0]; then 
# echo "Terminating..." >&2 ; 
# exit 1; 
#fi

#  $? 表示上一个命令的返回值，0 表示执行成功，非0 表示失败。
#  -ne 表示不等于，即如果 $? 不等于 0，则执行 if后面的语句。

if [ $? -ne 0 ]; then 
 usage && exit 1
fi

# Note the quotes around `$OPTIONS': they are essential!
#  eval可读取一连串的参数，而后再依参数自己的特性来执行相应的操作。
eval set -- "$OPTIONS"

## ---------- 设置默认值 ----------------------------------
FILE1=""
FILE2=""
READ_NUMS=2500
FORCE=false

## ---------- 处理解析后的命令行参数 ----------------------------------
while true; do
  case "$1" in
    -h|--help)
     usage
     ;;
    -f|--force)
        FORCE=true
        shift
        ;;
    -i|--input)
     FILE1="$2" # 赋值
     shift 2
     ;;
 -o|--output)
     FILE2="$2" # 赋值
     shift 2
     ;;
 -n|--read_nums)
  READ_NUMS=$2
  shift 2
     ;;
    --)
     shift
     break
     ;;
    *)
     echo "Cannot recognize: $1"
     usage
     ;;
  esac
done

## ---------- 参数的处理 ----------------------------------
if [[  -z $FILE1 ]]; then # 字符串长度是否为0
 echo "ERROR: Input file is required."
 usage && exit 1
fi

if [[ $READ_NUMS -le 0 ]];then # 小于等于则为真
 echo "ERROR: optinon '-n' must be a positive integer!"
 usage && exit 1
else
 num_lines=$(( READ_NUMS * 4 ))
fi

echo "$FILE1"
echo "$READ_NUMS"

## ---------- 流程接续用文件 ----------------------------------
<<EOF
this_scirpt=$(readlink -f "$0")
#echo "$this_scirpt" # 测试脚本所在位置。
end_log="${this_scirpt}.end_log"
EOF

## ---------- 定义函数 ----------------------------------
### 处理压缩文件
function get_fq_lines(){ # $1 fastq.gz file; $2 number of lines of fastq

 if [[ "$1" =~ \.gz$ ]]; then
  zcat "$1" |head -n "$2"
 else
  cat "$1" |head -n "$2"
 fi
}

#get_fq_lines "$FILE1" "$num_lines" # 简单测试下这个函数

### 检查fastq文件质量体系
function check_fq (){ # $1 fastq.gz file; $2 number of lines of fastq
 file=$(basename $1) # print filename only
 #zcat  $1 |\
 #head -n 10000 | \
 get_fq_lines $1 $2 | \
  awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | \
  awk -v file=$file 'BEGIN{min=100;max=0;} \
      {for(i=1;i<=NF;i++) \
          {if($i>max) max=$i; \
               if($i<min) min=$i;}}END \
          {if(max<=74 && min<59) \
                     print file"\tPhred+33"; \
           else \
           if(max>73 && min>=64) \
                     print file"\tPhred+64"; \
           else \
           if(min>=59 && min<64 && max>73) \
                     print file"\tSolexa+64"; else print file"\tUnknown score encoding!";}'
}

## ---------- 检查与运行 ----------------------------------
## 检查输入文件是否存在
if [ ! -f "$FILE1" ]; then
    echo "ERROR: $i not exits."
 usage
    exit 1
fi

echo "---------- step 1 start: `date`----------------------------------"
echo $FILE1 $num_lines  
check_fq $FILE1 $num_lines  # check phred score
echo "---------- step 1 done: `date`----------------------------------"

## ---------- 流程接续 ----------------------------------
<<EOF
#[[ -f $end_log ]] || ( echo "$FILE1" && touch ${end_log} )

if [[ $FORCE ]]; then rm -rf $end_log ;fi
if [[ ! -f $end_log ]] ;then
 for i in $(ls $FILE1); do
  check_fq $i $num_lines
 done         # 主体命令
fi

if [ $? -ne 0 ]; then
  usage && exit 1
 else
  touch ${end_log}
fi
EOF

echo "---------- DONE: `date`----------------------------------"