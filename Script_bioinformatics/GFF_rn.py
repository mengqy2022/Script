#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script is used to rename genes and mRNAs in a GFF3 file based on a given mapping and to add a numerical 
increment to the gene names. It supports the renaming of sequence identifiers to a specified prefix and ensures 
the uniqueness of identifiers by adding a numerical suffix.

Parameters:
    -g, --gff:         The input GFF3 file.
    -c, --change:      A file mapping original sequence names to the desired gene name prefixes.
    -a, --addnum:      The increment added to gene numbers to ensure uniqueness.
    -p, --prefix:      The prefix for the output GFF3 file.
"""
import re
import sys
import argparse
import logging
import gffutils
from collections import defaultdict

def rename(args):
    """
    Renames genes and mRNAs in the GFF3 file based on the provided mapping and adds a numerical suffix to gene names.

    Parameters:
        args:           Namespace object containing command line arguments.

    Returns:
        None:           This function does not return any value. It writes the renamed GFF3 data to a new file.
    """
    # Load the mapping from original sequence IDs to new gene name prefixes
    seqid2name = dict()
    for line in open(args.change, 'r'):
        tem = line.strip().split()
        seqid2name[tem[0]] = tem[1]

    # Create an in-memory database from the input GFF3 file
    db = gffutils.create_db(args.gff, ':memory:', force=True, keep_order=True, merge_strategy="create_unique", sort_attribute_values=True)
    # Define feature types that are children of mRNAs
    mRNA_children=("five_prime_UTR","three_prime_UTR","CDS","exon")
    # Map feature types to specific identifiers
    idmap = {
        "CDS" : "cds",
        "exon" : "exon",
        "five_prime_UTR":"5utrp",
        "three_prime_UTR":"3utrp"
    }

    # Open the output GFF3 file for writing
    f_out = open('%s.rename.gff3' % args.prefix, 'w')
    seqid = None
    # Iterate over genes, renaming and writing them to the output file
    for gene in db.features_of_type("gene", order_by=('seqid','start','end')):
        if gene.seqid != seqid:
            genenum = 0
        genenum += args.addnum
        seqid = gene.seqid

        # Generate the new gene name
        genename = '{0}G{1:06}'.format(seqid2name[seqid], genenum)
        # Write the gene record to the output file
        f_out.write('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={geneid}\n'.format(seqid=gene.seqid, source=gene.source, featuretype=gene.featuretype, start=gene.start, end=gene.end, score=gene.score, strand=gene.strand, frame=gene.frame, geneid=genename))

        # Iterate over mRNAs belonging to the gene, rename, and write them
        for t, mRNA in enumerate(db.children(gene, featuretype="mRNA", order_by=('seqid','start','end'))):
            mrna_num = t + 1
            mrnaid = '{genename}.mRNA{num}'.format(genename=genename, num=mrna_num)
            # Write the mRNA record to the output file
            f_out.write('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={mrnaid};Parent={geneid}\n'.format(seqid=mRNA.seqid, source=mRNA.source, featuretype=mRNA.featuretype, start=mRNA.start, end=mRNA.end, score=mRNA.score, strand=mRNA.strand, frame=mRNA.frame, mrnaid=mrnaid, geneid=genename))

            # Manage child features of mRNAs (e.g., CDS, UTRs), rename, and write them
            numdict = defaultdict(int)
            for child in db.children(mRNA, featuretype=mRNA_children, order_by=("start",'end')):
                numdict[child.featuretype] += 1
                child_id = '{genename}.{childid}{num}'.format(genename=genename, childid=idmap[child.featuretype], num=numdict[child.featuretype])
                # Write the child feature record to the output file
                f_out.write('{seqid}\t{source}\t{featuretype}\t{start}\t{end}\t{score}\t{strand}\t{frame}\tID={child_id};Parent={mrnaid}\n'.format(seqid=child.seqid, source=child.source, featuretype=child.featuretype, start=child.start, end=child.end, score=child.score, strand=child.strand, frame=child.frame, child_id=child_id, mrnaid=mrnaid))

    # Close the output file
    f_out.close()

def main():
    """
    Sets up logging and parses command line arguments before calling the rename function.

    Parameters:
        None:           This function does not take any parameters.

    Returns:
        None:           This function does not return any value.
    """
    # Configure logging
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    # Create the argument parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="rename gff3 file"
    )

    # Add command line arguments
    parser.add_argument('-g', '--gff', required=True, help='gff3 file')
    parser.add_argument('-c', '--change', required=True, help='a file, correspondence between sequence name and gene name prefix')
    parser.add_argument('-a', '--addnum', type=int, default=1, help='diff in gene number, such as if addnum = 10, xx1G000010, xx1G000020')
    parser.add_argument('-p', '--prefix', default='result', help='prefix of output')
    args = parser.parse_args()

    # Call the renaming function with the parsed arguments
    rename(args)


if __name__ == "__main__":
    main()