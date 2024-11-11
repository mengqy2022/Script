#!/usr/bin/perl -w

#该脚本用于转换Blast比对结果的格式
#同事写的，个人懒得再换个其它语言写个类似的转换脚本，就一直在用他的

use strict;
use Bio::SearchIO;

if(@ARGV != 2)  {
    print "perl blast_parser.pl INPUT  OUTPUT \n";
	exit;
}

my($input, $output) = @ARGV;

open(OUTPUT, ">$output") || die "could not open.\n";
print OUTPUT "Query_name\tQuery_description\tQuery_length\tQuery_start\tQuery_end\tHit_name\tHit_description\tHit_length\tHit_start\tHit_end\tAln_length\tIdentity\Evalue\n";
my $searchio = Bio::SearchIO->new(-format=>"blast", -file => "$input");

while(my $result = $searchio->next_result)  {
    my $query_name = $result->query_name;
    my $query_description = $result->query_description;
    my $query_length = $result->query_length;
    if($result->num_hits) {
    my $hit = $result->next_hit;
        my $hit_name = $hit->name();
        my $hit_desc = $hit->description();
        my $hit_length = $hit->length;
        my $hsp = $hit->next_hsp;
            my $evalue = $hsp->evalue;
            my $strand = $hsp->strand('hit');
            my $query_start = $hsp->start('query');
            my $query_end = $hsp->end('query');
            my $hit_start = $hsp->start('hit');
            my $hit_end = $hsp->end('hit');
            my $aln_length = $hsp->length('total');
            my $identity = $hsp->percent_identity;
            print OUTPUT $query_name,"\t",$query_description,"\t",$query_length,"\t",$query_start,"\t",$query_end,"\t",$hit_name,"\t",$hit_desc,"\t",$hit_length,"\t",$hit_start,"\t",$hit_end,"\t",$aln_length,"\t",$identity,"\t",$evalue,"\n";            
    }
}
close OUTPUT;
exit;
