#!/usr/bin/perl

##############################################################################
#
#                                 cmpfastq
#
# Concept: Stephen Newhouse (stephen.newhouse@kcl.ac.uk)
# Author:  David To (david.to@kcl.ac.uk)
# Written: 30th June 2010
# 
# DESCRIPTION:
# This script is designed to compare two fastq files.
# It produces 4 output files
# - 2 files that is a list of common sequences of both files
# - 2 files that are unique to each file
#
# Known issues: This script appears to use a lot of memory for many
#               sequences (still much less than reading all sequences into
#               memory).
#
# Copyright 2010 NIHR Biomedical Research Centre for Mental Health
#                South London and Maudsley NHS Foundation Trust &
#                Institute of Psychiatry
#                Kings College London
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see:
# <http://www.gnu.org/licenses/>.
##############################################################################

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

## CONSTANTS ##
my $TRUE              = 1;
my $FALSE             = 0;
my $DEBUG             = $FALSE;
my $EXITSTATUS        = 0;

# Default umask
umask 027;

# Define variables

# Get the options passed at the command prompt
GetOptions();

##############################
#   M A I N   P R O G R A M  #
##############################

# Check to see we received two files in the arguments
if(scalar(@ARGV) != 2)
{
    print STDERR "Incorrect number of arguments\n";
    Usage();
    exit(1);
}

my $fail = $FALSE;

# Check to see if the files exist
foreach my $file (@ARGV)
{
    if(!-e $file)
    {
        print STDERR "File $file didn't exist\n";
        $fail = $TRUE;
    }
}

# If any of the files didn't exist， let's kill it
if($fail)
{
    exit(1);
}

# Read the file names in from the command line
my $file1 = shift(@ARGV);
my $file2 = shift(@ARGV);

# Index the first file.
my %fastqIndex1 = %{IndexFastq($file1)};

# Compare the two files
CompareFastq($file1, $file2, \%fastqIndex1);

exit($EXITSTATUS);

# Subroutines
sub Usage
{
    my $base = basename($0);
    print "Usage: $base [dh] file1 file2\n";
    print "\td:\tDebug mode on (default off)\n";
    print "\th:\tPrint this usage\n";
    print "\nDESCRIPTION:\n";
    print "\tThis script compares two FASTQ files and produces output files that list:\n";
    print "\t- Common sequences found in both files\n";
    print "\t- Unique sequences found in each file\n";
    print "\nEXAMPLES:\n";
    print "\t$base file1.fastq file2.fastq\n";
    print "\t$base -d file1.fastq file2.fastq  # Enable debug mode\n";
}

sub GetOptions
{
    # Get the options passed at the command prompt
    my %options=();
    getopts("dh", \%options);

    if(defined($options{'d'}))
    {
        $DEBUG = $TRUE;
    }

    if(defined($options{'h'}))
    {
        Usage();
        exit($EXITSTATUS);
    }
}

sub IndexFastq
{
    my $file = shift;
    my %fastqIndex;

    open(IN, $file) or die("Could not open $file\n");
    my $pos = tell(IN);
    my $lineCounter = 1;
    while(my $line = <IN>)
    {
        chomp($line);

        # Each block is going to be of 4 lines
        # Let's get the seq ID from the sequence name
        if($line =~ m/^@(.*)#.*/)
        {
            $fastqIndex{$1} = $pos;
            # Skip the next 3 lines
            for(my $i=0; $i<3; $i++)
            {
                <IN>;
                $lineCounter++;
            }
        }
        elsif($line =~ m/^#/)
        {
            print STDERR "File: $file\[$lineCounter]: Skipping comment line: $line\n" if($DEBUG);
        }
        elsif($line =~ m/^$/)
        {
            print STDERR "File: $file\[$lineCounter]: Skipping empty line: $line\n" if($DEBUG);
        }
        else
        {
            print STDERR "File: $file\[$lineCounter]: Could not match the sequence ID from the name: $line\n" if($DEBUG);
        }
        $pos = tell(IN);
        $lineCounter++;
    }
    close(IN);

    return \%fastqIndex;
}

sub CompareFastq
{
    my $file1          = shift;
    my $file2          = shift;
    my $fastqIndex1Ref = shift;
    my %fastqIndex1    = %{$fastqIndex1Ref};
    my %found1;

    # We don't want to have to open/close file handles each time, so let's open them here
    open(F1COUT, ">$file1-common.fq") or die("Could not write to file: $file1-common.fq\n");
    open(F2COUT, ">$file2-common.fq") or die("Could not write to file: $file2-common.fq\n");
    open(F1UOUT, ">$file1-unique.fq") or die("Could not write to file: $file1-unique.fq\n");
    open(F2UOUT, ">$file2-unique.fq") or die("Could not write to file: $file2-unique.fq\n");


    open(F1IN, $file1) or die("Could not open $file1\n");
    open(F2IN, $file2) or die("Could not open $file2\n");
    while(my $line = <F2IN>)
    {
        chomp($line);

        # Skip empty lines or comments
        if($line =~ m/^$/g or $line =~ m/^\s*#/)
        {
            next;
        }

        # Each block is going to be of 4 lines
        # Let's get the seq ID from the sequence name
        if($line =~ m/^@(.*)#.*/)
        {
            my $seqId = $1;

            if(defined($fastqIndex1{$seqId}))
            {
                $found1{$seqId} = $TRUE;

                # Print out from file1
                seek(F1IN, $fastqIndex1{$seqId}, 0);
                for(my $i=0;$i<4;$i++)
                {
                    my $tmpLine = <F1IN>;
                    print F1COUT $tmpLine;
                }
                
                # Print out from file 2
                print F2COUT $line . "\n";
                for(my $i=0; $i<3; $i++)
                {
                    my $tmpLine = <F2IN>;
                    print F2COUT $tmpLine;
                }
            }
            else
            {
                # Print out from file 2
                print F2UOUT $line . "\n";
                for(my $i=0; $i<3; $i++)
                {
                    my $tmpLine = <F2IN>;
                    print F2UOUT $tmpLine;
                }
            }
        }
        else
        {
            print STDERR "Could not match the sequence ID from the name: $line\n";;
            next;
        }
    }
    close(F1COUT);
    close(F2COUT);
    close(F2UOUT);
    close(F2IN);
    # Now let's worry about the sequences that weren't common in file 1

    # File 1
    if(keys(%fastqIndex1) != keys(%found1))
    {
        foreach my $seqId (keys %fastqIndex1)
        {
            if(!defined($found1{$seqId}))
            {
                seek(F1IN, $fastqIndex1{$seqId}, 0);
                for(my $i=0;$i<4;$i++)
                {
                    my $tmpLine = <F1IN>;
                    print F1UOUT $tmpLine;
                }
            }
        }
    }
    close(F1UOUT);
    close(F1IN);
}
