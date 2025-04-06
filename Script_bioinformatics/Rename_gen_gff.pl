#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

=pod

=head1 NAME

rename_contigs.pl - 重命名FASTA序列ID并同步更新GFF3文件中的对应ID

=head1 SYNOPSIS

rename_contigs.pl [选项]

 必选参数:
  -f, --fasta         输入的FASTA文件
  -g, --gff3          输入的GFF3文件
  -o, --output_fasta  输出的重命名后FASTA文件
  -p, --output_gff3   输出的更新后GFF3文件

 其他选项:
  -h, --help          显示此帮助信息
  -v, --version       显示版本信息

=head1 DESCRIPTION

本脚本执行以下操作:
1. 将FASTA文件中的序列ID重命名为Contig1, Contig2...格式
2. 同步更新GFF3文件中对应的旧ID
3. 确保GFF3文件中所有相关引用都被正确更新

=head1 EXAMPLES

基本用法:
  rename_contigs.pl -f input.fasta -g input.gff3 -o renamed.fasta -p updated.gff3

显示帮助:
  rename_contigs.pl -h

=cut

# 版本信息
our $VERSION = '1.0.0';

# 重命名fasta文件中的序列ID
sub rename_fasta {
    my ($input_fasta, $output_fasta) = @_;
    my %id_mapping;
    my $seq_in = Bio::SeqIO->new(-file => $input_fasta, -format => 'fasta')
        or die "无法打开FASTA输入文件 $input_fasta: $!";
    my $seq_out = Bio::SeqIO->new(-file => ">$output_fasta", -format => 'fasta')
        or die "无法创建FASTA输出文件 $output_fasta: $!";
    
    my $i = 1;
    while (my $seq = $seq_in->next_seq) {
        my $old_id = $seq->id;
        my $new_id = "Contig$i";
        $id_mapping{$old_id} = $new_id;
        
        $seq->id($new_id);
        $seq->desc('');  # 清空描述信息
        $seq_out->write_seq($seq);
        $i++;
    }
    
    return \%id_mapping;
}

# 更新gff3文件
sub update_gff3 {
    my ($input_gff3, $output_gff3, $id_mapping) = @_;
    open(my $in, '<', $input_gff3) or die "无法打开GFF3输入文件 $input_gff3: $!";
    open(my $out, '>', $output_gff3) or die "无法创建GFF3输出文件 $output_gff3: $!";
    
    while (my $line = <$in>) {
        if ($line =~ /^#/) {
            print $out $line;
            next;
        }
        
        chomp $line;
        my @fields = split(/\t/, $line);
        if (scalar @fields < 9) {
            print $out $line, "\n";
            next;
        }
        
        # 替换第一列的ID
        my $old_seqid = $fields[0];
        if (exists $id_mapping->{$old_seqid}) {
            $fields[0] = $id_mapping->{$old_seqid};
        }
        
        # 替换第九列属性中的ID
        my $attributes = $fields[8];
        foreach my $old_id (keys %$id_mapping) {
            my $new_id = $id_mapping->{$old_id};
            $attributes =~ s/\b\Q$old_id\E\b/$new_id/g;
        }
        $fields[8] = $attributes;
        
        print $out join("\t", @fields), "\n";
    }
    
    close $in;
    close $out;
}

# 主程序
sub main {
    my %args;
    GetOptions(
        'f|fasta=s'         => \$args{fasta},
        'g|gff3=s'          => \$args{gff3},
        'o|output_fasta=s'  => \$args{output_fasta},
        'p|output_gff3=s'   => \$args{output_gff3},
        'h|help'            => sub { pod2usage(1) },
        'v|version'         => sub { print "rename_contigs.pl 版本 $VERSION\n"; exit },
    ) or pod2usage(2);
    
    # 检查必需参数
    for my $req (qw(fasta gff3 output_fasta output_gff3)) {
        pod2usage("错误: 缺少必需参数 --$req") unless defined $args{$req};
    }
    
    print "开始处理FASTA文件...\n";
    my $id_mapping = rename_fasta($args{fasta}, $args{output_fasta});
    my $count = scalar keys %$id_mapping;
    print "成功重命名 $count 条序列\n";
    
    print "开始更新GFF3文件...\n";
    update_gff3($args{gff3}, $args{output_gff3}, $id_mapping);
    print "GFF3文件更新完成\n";
    
    print "所有操作完成！\n";
}

main();
