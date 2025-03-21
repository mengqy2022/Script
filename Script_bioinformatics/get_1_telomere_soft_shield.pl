#!/usr/bin/perl -w 

use strict;

my $fasta=$ARGV[0];

open (FASTA,$fasta);
my %hash;
my $tag;
while (<FASTA>){
        $_=~s/[\r\n]+//g;

        if ($_=~/^\>/){
                ($tag) =($_ =~ /^\>(\S+)/);
        }else{
                 $hash{$tag}.=$_;
                }

}
foreach my $tag(keys %hash) {
		if($hash{$tag}=~/[ca]*ccccaaaacccc/ or $hash{$tag}=~/ggggttttgggg[gt]*/)
	{print "\>$tag\n$hash{$tag}\n"}

}
