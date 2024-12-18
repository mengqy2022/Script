#!/usr/bin/perl -w

#�ýű���������С��fasta�ļ�֮�����е��໥blast
#ͬ��д�ģ����������ٻ�����������д�����Ƶ�ת���ű�����һֱ��������

use strict;
use Bio::SearchIO;

if(@ARGV != 2)  {
    print "perl blast_parser.pl INPUT  Ref \n";
	exit;
}

my($input, $ref) = @ARGV;
my $db = $ref;
system("makeblastdb -in $ref -dbtype nucl");
my $blastout = $input . ".blast";
my $output = $blastout . ".result";
system("blastn -query $input -db $ref -out $blastout -evalue 1e-5");

open(OUTPUT, ">$output") || die "could not open.\n";
print OUTPUT "Query_name\tQuery_description\tQuery_length\tQuery_start\tQuery_end\tHit_name\tHit_description\tHit_length\tHit_start\tHit_end\tAln_length\tIdentity\tHit_Strand\tE_value\n";
my $searchio = Bio::SearchIO->new(-format=>"blast", -file => "$blastout");

while(my $result = $searchio->next_result)  {
    my $query_name = $result->query_name;
    my $query_description = $result->query_description;
    my $query_length = $result->query_length;
    while(my $hit = $result->next_hit)  {
        my $hit_name = $hit->name();
        my $hit_desc = $hit->description();
        my $hit_length = $hit->length;
        while(my $hsp = $hit->next_hsp)  {
            my $evalue = $hsp->evalue;
            my $strand = $hsp->strand('hit');
            my $query_start = $hsp->start('query');
            my $query_end = $hsp->end('query');
            my $hit_start = $hsp->start('hit');
            my $hit_end = $hsp->end('hit');
            my $aln_length = $hsp->length('total');
            my $identity = $hsp->percent_identity;
            print OUTPUT "$query_name\t$query_description\t$query_length\t$query_start\t$query_end\t$hit_name\t$hit_desc\t$hit_length\t$hit_start\t$hit_end\t$aln_length\t$identity\t$strand\t$evalue\n";            
        }
    }
}
close OUTPUT;
exit;
