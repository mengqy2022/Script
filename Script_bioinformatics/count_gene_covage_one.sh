#!/bin/bash/
# **************************************************
# 先利用sedkit计算基因组和预测基因总长度，在计算比率，输出列表格式文件
# Date   : 2024-03-11
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    filename=$(basename $0)
    echo -e "\n             >>>> ${filename} <<<< \n " 
    echo -e "usage: bash ${filename} [Genomics FOLDER] [Protein coding sequence FOLDER]\n"
    echo -e "Genomics FOLDER: the folder contains the genomics data, including the fasta file.\n"
    echo -e "Protein coding sequence FOLDER: the folder contains the protein coding sequence data.\n"
    echo -e "example: bash ${filename} /data/database/foldername /data/protein_coding_sequence_foldername"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

#  find找到文件并复制到另一个文件夹
# find ./prokka -name "*\.ffn" -exec cp {} ./cds_databas \;

pwd=$(pwd)
echo -e "\nCurrent directory: $pwd\n"

if [ ! -d $1 ] ;then
    echo "`date '+[error] %D %T'`"
    echo -e "[Genomics FOLDER: $1 does not exist.]\n"
    echo -e "Please check the folder name and try again. \n"
fi

if [ ! -d $2 ] ;then
    echo "`date '+[error] %D %T'`"
    echo -e "[Protein coding sequence FOLDER: $2 does not exist.]\n"
    echo -e "Please check the folder name and try again. \n"
fi

if [ -f $pwd/genomics.stats ] ;then
    rm -rf   $pwd/genomics.stats
fi

if [ -f $pwd/protein.stats ] ;then
    rm -rf   $pwd/protein.stats
fi

cd $pwd/$1

pwd_1=$(pwd)
echo -e "Current directory: $pwd_1\n"

for genomics in $(ls $pwd/$1);do
    seqkit stats -T $genomics | csvtk cut -t -f 1,5 | csvtk del-header >> $pwd/genomics.stats
done

cd $pwd/$2

pwd_2=$(pwd)
echo -e "Current directory: $pwd_2\n"

for protein in $(ls $pwd/$2);do
    seqkit stats -T $protein | csvtk cut -t -f 1,5 | csvtk del-header >> $pwd/protein.stats
done

cat $pwd/genomics.stats | sed 's/\.fna//g' | sed 's/\.fasta//g' | sed 's/\.fa//g' | sed 's/\.fnd//g' | sed 's/\.faa//g' | sed 's/\.ffn//g' | sed 's/\./_/g' | awk -F"\t" '{print $1"\t"$2}' > 123 && mv 123 $pwd/genomics.stats

cat $pwd/protein.stats | sed 's/\.fna//g' | sed 's/\.fasta//g' | sed 's/\.fa//g' | sed 's/\.fnd//g' | sed 's/\.faa//g' | sed 's/\.ffn//g' | sed 's/\./_/g' | awk -F"\t" '{print $1"\t"$2}' > 123 && mv 123 $pwd/protein.stats

if [ -f $pwd/cds_coverage.txt ] ;then
    rm -rf $pwd/cds_coverage.txt
fi

paste $pwd/genomics.stats $pwd/protein.stats | cut -f 1-2,4 | awk -v OFS="\t" '{print $1,$3/$2*100}' > $pwd/cds_coverage.txt

cd $pwd

awk 'BEGIN{print "name\tratio"; while (getline < "cds_coverage.txt") print}' > 123 && mv 123 cds_coverage.txt

rm -rf $pwd/genomics.stats $pwd/protein.stats

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are stored in $pwd \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"