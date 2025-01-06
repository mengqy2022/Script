#!/usr/bin/env python
# # -*- coding: utf-8 -*-
import sys
import argparse
import logging
import gffutils
from collections import defaultdict

def rename(args):
    seqid2name = {}
    # 读取变更映射文件，并捕获错误
    try:
        with open(args.change, 'r') as change_file:
            for line in change_file:
                tem = line.strip().split()
                if len(tem) >= 2:
                    seqid2name[tem[0]] = tem[1]
    except Exception as e:
        logging.error(f"读取变更文件 {args.change} 失败: {e}")
        sys.exit(1)

    # 创建GFF数据库
    try:
        db = gffutils.create_db(args.gff, ':memory:', force=True, keep_order=True, merge_strategy="create_unique", sort_attribute_values=True)
    except Exception as e:
        logging.error(f"创建GFF数据库失败: {e}")
        sys.exit(1)

    mRNA_children = ("exon", "CDS")
    idmap = {"CDS": "CDS", "exon": "exon"}

    # 使用with语句管理文件打开，确保正确关闭文件
    with open(f'{args.prefix}.rename.gff3', 'w') as f_out:
        genenum = 10  # 初始化基因编号
        for gene in db.features_of_type("gene", order_by=('seqid', 'start', 'end')):
            seqid = gene.seqid
            genename = f'{seqid2name.get(seqid, "unknown")}g{genenum:06}'
            f_out.write(f'{gene.seqid}\t{gene.source}\t{gene.featuretype}\t{gene.start}\t{gene.end}\t{gene.score}\t{gene.strand}\t{gene.frame}\tID={genename}\n')

            for t, mRNA in enumerate(db.children(gene, featuretype="mRNA", order_by=('seqid', 'start', 'end'))):
                mrna_num = t + 1
                mrnaid = f'{seqid2name.get(seqid, "unknown")}t{genenum:06}.{mrna_num}'
                f_out.write(f'{mRNA.seqid}\t{mRNA.source}\t{mRNA.featuretype}\t{mRNA.start}\t{mRNA.end}\t{mRNA.score}\t{mRNA.strand}\t{mRNA.frame}\tID={mrnaid};Parent={genename}\n')

                numdict = defaultdict(int)
                for child in db.children(mRNA, featuretype=mRNA_children, order_by=("start", 'end')):
                    numdict[child.featuretype] += 1
                    child_id = f'{mrnaid}.{idmap[child.featuretype]}{numdict[child.featuretype]}'
                    f_out.write(f'{child.seqid}\t{child.source}\t{child.featuretype}\t{child.start}\t{child.end}\t{child.score}\t{child.strand}\t{child.frame}\tID={child_id};Parent={mrnaid}\n')
                
                genenum += 10  # 为下一个基因递增编号

def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="rename gff3 file"
    )
    parser.add_argument('-g', '--gff', required=True, help='gff3 file')
    parser.add_argument('-c', '--change', required=True, help='a file, correspondence between sequence name and gene name prefix')
    parser.add_argument('-a', '--addnum', type=int, default=1, help='diff in gene number, such as if addnum = 10, xx1G000010, xx1G000020')
    parser.add_argument('-p', '--prefix', default='result', help='prefix of output')
    args = parser.parse_args()

    rename(args)

if __name__ == "__main__":
    main()
