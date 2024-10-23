#!/usr/bin/perl
use strict;

print STDERR "This script calculate each seq's GC content, and make a distribution\n\n\n";


######################################hashing fasta############
print STDERR "Reading the fasta file\n";
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
##########################
my $max=0;
my $min=1000000000000000;
my $tolseqs;
my $k50;
my $k5;
my $tollen;
foreach my $key (keys %hash){
	my $len=length($hash{$key});
	$tollen+=$len;
	if ($len>$max){
		$max=$len;
	}
	
	if ($len<$min){
		$min=$len;
	}
	
	if ($len>=50000){
		$k50++;
	}
	
	if ($len<=5000){
		$k5++;
	}
	
	$tolseqs++;
	
}

my $out=$ARGV[0] . "\.stats";
open OUT, ">$out";
print OUT "$tolseqs\t$max\t$min\t$k50\t$k5\t$tollen\n";


my $out2=$ARGV[0] . "\.perseqs";
open OUT2, ">$out2";
#foreach seqs file, calculate the GC content
my %gchis;
foreach my $key (keys %hash){
	$key=~s/[\r\n]+//g;
	$hash{$key}=uc($hash{$key});
	my $tmp=$hash{$key};
	$tmp=~s/N//g;
	my $ln1=length($hash{$key});
	my $ln2=length($tmp);
	
	(my $gc)=($tmp=~tr/GC//);
	my $gc_content=sprintf( "%.2f", $gc/length($tmp));
	$gchis{$gc_content}+=1;
	print OUT2 "$key\t$ln1\t$ln2\t$gc_content\n";	
}

#############################################

my $out1=$ARGV[0] . "\.GCs";
open OUT1, ">$out1";
#######################################
############print files
foreach my $key (sort {$a<=>$b} keys %gchis){
	print OUT1 "$key\t$gchis{$key}\n";
	
}
#######################################




