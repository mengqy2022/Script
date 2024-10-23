#!/usr/bin/perl -w 

use strict;
#########################################################################################
#USAEGE:using cmd (current directory), command: "perl getseq.pl fastafile YOUR-ID-file" #
#########################################################################################

my $fasta=$ARGV[0];                                                                         
my $ids=$ARGV[1];

open (ID,$ids);
my @id=<ID>;

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

my $out=$ids . "\.fasta";
open O,">>$out";

my $c=0;
foreach my $uid(@id){
	$uid=~s/[\r\n \t]+//g;
	foreach my $key(keys %hash){
		if ($uid eq $key){
				print O "\>$uid\n$hash{$key}\n";
				$c++;
		}		
	}	
}

my $in_id_no=@id;

if ($in_id_no==$c){
	print "\n\nYOUR INPUT ID NUMBER IS: $in_id_no\; HITS NO IS : $c\;\n\n"; 
}else{
	print "\nYOUR INPUT ID NUMBER IS: $in_id_no\;\n\nHITS NUMBER IS : $c\;\n\n!!!WARNNING!!!:YOUR INPUT IDS NUMBER NOT EQUAL TO HITS, CHECK YOUR IDS!!!\n";	
}
