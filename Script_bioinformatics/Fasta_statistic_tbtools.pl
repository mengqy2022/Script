#!/usr/bin/perl

use strict;
use Getopt::Long;

my ($input_file,$output_file,$help);
GetOptions(
    "i:s"=>\$input_file,
    "o:s"=>\$output_file,
    "h|help"=>\$help,
    );
&usage if (defined $help || !$input_file || !$output_file);
sub usage {
    die "
        Usage: 
          此脚本用于统计核酸或氨基酸序列文件，统计时去除序列中可能存在的*，匹配GAP(N)、ATCG以及氨基酸序列时均不区分大小写;
          统计包括文件的序列条数，N50，每条序列和整体序列的长度、GC含量、GAP(N)含量，最长和最短序列长度及所有对应的序列；
          如果是氨基酸序列则GAP和GC含量可忽略；
          如果是核酸序列，则：单条序列长度、总序列长度、N50长度、平均长度、最长最短序列长度计算均包括GAP(N)，但单条序列GC含量、整体核酸序列GC含量为GC碱基占总ATCG碱基的比例（去除GAP）
          perl fasta_statistic.pl <options>
            
        Version 1.01 Feb 02,2023
        Version 2.01 Jan 16,2024
        
        Options:
        -i      <str> 输入文件名
        -o      <str> 输出文件名
        -h|help      <str> Pod this help
        
        eg:perl fasta_statistic.pl -i genome.fasta -o genome_stats \n";
}

our @temp;
our $total_num=0;
our $total_len=0;
our $total_gc_number=0;
our $total_gap_number=0;
our $total_atcg_number=0;
our $last;
my @len_array=();
my %seq_length;
my $possible;
my $cal50=0;
my $N50;

our $output = "ID\tGeneLength(bp)\tGAPLength(bp)\tGAP_ratio(%)\tGClength(bp)\tGCcontent(%)(Exclude N)\n";
open IN,"<$input_file",or die "input_file not exists!\n";
open OUT, ">$output_file",or die "output_file not exists!\n";
$/=">";<IN>;
while (<IN>) {
    #去除分隔符 >
    chomp;
    
    #判断fasta文件中序列是否为一行输出，如果不是则改为一行输出
    #将序列按照换行符储存在数组，在第一个序列ID和最后一行序列末尾添加一个换行符
    @temp=(split /\n/,$_);
    $last = $#temp;
    for (0..$last)
    {
       if ($_==0 or $_==$#temp)
        {
                $temp[$_] .="\n";}
        else{next;}
    }
    my $new = join("",@temp);
    
    my ($id,$seq)=(split /\n/,$new,2)[0,1];
    @temp = ();
    
    #删除最后换行符和可能存在的*
    $seq=~ s/\*//g;
    $seq=~ s/\n//g;
    #大小长度
    my $len = length ($seq);
    $total_len+=$len;
    $seq_length{$id} = $len;
    push @len_array,$len;
    
    #小gap数目比例，大gap数目
    my $gap_number =($seq=~s/N/N/ig) || 0;
    my $gap_ratio=($gap_number/$len)*100;
    $total_gap_number +=$gap_number;
   #小gc数目比例，大gc数目
   #小atcg数目比例，大atcg数目
    my $gc_number =($seq=~s/G/G/ig+$seq=~s/C/C/ig);
    my $gc_ratio=($gc_number/($len-$gap_number))*100;
    $total_gc_number +=$gc_number;
    my $atcg_number = ($seq=~s/G/G/ig+$seq=~s/C/C/ig+$seq=~s/A/A/ig+$seq=~s/T/T/ig);
    $total_atcg_number +=$atcg_number;
    $output .= "$id\t$len\t$gap_number\t";
    $output .= sprintf "%.2f\t","$gap_ratio";
    $output .= "$gc_number\t";
    $output .= sprintf "%.2f\n","$gc_ratio";
    $total_num++;
}

#总长和总GC
my $avg_len=($total_len/$total_num);
my $total_GC=($total_gc_number/$total_atcg_number)*100;

#最长和最短
my ($max_seq, $max_length);
my @max_seqs;
my $min_seq;
my @min_seqs;

#获取所有长度
my @lengths = values %seq_length;
#找到最大和最小长度
my $max_length = (sort { $b <=> $a } @lengths)[0];
my $min_length = (sort { $a <=> $b } @lengths)[0];

#找到所有最大长度对应的序列
foreach my $max_seq_ID (keys %seq_length) {
    if ($seq_length{$max_seq_ID} == $max_length) {
        push @max_seqs, $max_seq_ID;
    }
}
#找到所有最小长度对应的序列
foreach my $min_seq_ID (keys %seq_length) {
    if ($seq_length{$min_seq_ID} == $min_length) {
        push @min_seqs, $min_seq_ID;
    }
}

#N50
@len_array=sort {$a<=>$b} @len_array;
foreach $possible(@len_array) 
{
  $cal50 += $possible;
  if ($cal50 >= $total_len/2) {
    $N50=$possible;
    last;
  }
}
#找出 N50 对应的所有序列名
my @N50_seqs;
foreach my $seq_name (keys %seq_length) {
    if ($seq_length{$seq_name} == $N50) {
        push @N50_seqs, $seq_name;
    }
}

$output .= "\n========================================================\n";
$output .= "Total State\n";
$output .= "Total sequence number:$total_num\n";
$output .= "Total Length:$total_len bp\n";
$output .= "N50:$N50 bp\n";
$output .= "Sequences with length equal to N50: " . join(', ', @N50_seqs) . "\n";
$output .= sprintf "Average Length:%.2f\n","$avg_len bp";
$output .= "Gap(N):$total_gap_number bp\n";
$output .= sprintf "GC content (Exclude N):%.2f%\n","$total_GC";
$output .= "Maximum Length:$max_length bp\n";
$output .= "Minimum Length:$min_length bp\n";
$output .= "Sequences with Max Length: " . join(', ', @max_seqs). "\n";
$output .= "Sequences with Min Length: " . join(', ', @min_seqs). "\n";

print OUT $output;