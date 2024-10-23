#!/bin/bash/
# **************************************************
# 提取基因组中的16s序列，细菌
# Date   : 2024-03-11
# Author : 孟庆瑶
# **************************************************

set -e

if [ -z $1 ] ;then
    echo "`date '+[error] %D %T'`"
    filename=$(basename $0)
    echo -e "\nCharacterization and extraction of ribosomal genomes. [16S rRNA] [Bacteria]\n" 
    echo -e "usage: bash ${filename} [DATA BASE FOLDER] [OUTPUT FOLDER]"
    echo -e "example: bash ${filename} /data/database/foldername"
    echo -e "\n             >>>> mqy <<<< \n "
    exit 1
fi

pwd=$(pwd)
echo "Current directory: $pwd"

if [ -a $pwd/$2 ] ;then
    rm -rf  $pwd/$2
fi

mkdir $pwd/$2

for filename in $(ls $1 | sed 's/\.[^.]*$//'); do
  barrnap $pwd/$1/${filename}\.* > $pwd/$1/${filename}_16s.gff
  cat $pwd/$1/${filename}_16s.gff | grep "16S" | head -n 1 > temp && mv temp $pwd/$1/${filename}_16s.gff
  gff2bed --do-not-sort < $pwd/$1/${filename}_16s.gff > $pwd/$1/${filename}_16s.bed
  seqkit subseq --bed $pwd/$1/${filename}_16s.bed $pwd/$1/${filename}\.* >$pwd/$2/${filename}_16s.fasta
  rm $pwd/$1/${filename}_16s.gff $pwd/$1/${filename}_16s.bed $pwd/$1/*.fai
done

echo -e " \n >>>> Done <<<< \n "
echo -e "Results are stored in $pwd/$2 \n"
echo -e "Thank you for using our software! \n"
echo -e "Please contact us if you have any questions. \n"
echo -e "Email: <15877464851@163.com> \n"
echo -e "Sincerely, \n"
echo -e "MQY \n"