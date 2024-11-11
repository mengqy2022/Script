#!/usr/bin/perl
use strict;
my %hash;
my $tag;
open (TEXT,$ARGV[0]);
while(<TEXT>)
{
	chomp;
	my @ar=split(/[\t\s]/,$_);
	$tag=$ar[0];
	$hash{$tag}.="$_\(";
}
my %hs;
foreach my $key (keys %hash)
{
	my @ar=split(/\(/,$hash{$key});
	$hs{$key}.=$ar[0];
}
foreach my $id (keys %hs)
{
	print "$hs{$id}\n";
}

