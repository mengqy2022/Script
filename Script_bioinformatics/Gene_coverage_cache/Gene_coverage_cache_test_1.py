#!/usr/bin/env python
import os
import pytest
from Gene_coverage_cache import read_fasta_lengths, calculate_cds_proportions, write_results_to_file, process_file_pair, process_file_pairs, main

# 为read_fasta_lengths函数编写单元测试
def test_read_fasta_lengths():
    # 假设有一个合法的fasta文件，其中包含两个序列
    test_fasta = "test.fasta"
    with open(test_fasta, 'w') as f:
        f.write(">seq1\nATCG\n>seq2\nCGTA")
    
    # 调用函数
    lengths = read_fasta_lengths(test_fasta)
    
    # 断言结果是正确的
    assert lengths == {'seq1': 4, 'seq2': 4}
    
    # 清理创建的测试文件
    os.remove(test_fasta)

# 为calculate_cds_proportions函数编写单元测试
def test_calculate_cds_proportions():
    # 构造输入数据
    genome_lengths = {'genome1': 1000, 'genome2': 2000}
    cds_lengths = {'cds1': 400, 'cds2': 800}
    
    # 调用函数
    proportions = calculate_cds_proportions(genome_lengths, cds_lengths)
    
    # 断言结果是正确的
    assert proportions == {'cds1': 0.4, 'cds2': 0.4}

# 为write_results_to_file函数编写单元测试
def test_write_results_to_file(tmp_path):
    # 构造输入数据
    output_file = tmp_path / "output.txt"
    results = {'file1': {'cds1': 0.5}, 'file2': {'cds2': 0.25}}
    
    # 调用函数
    write_results_to_file(str(output_file), results)
    
    # 检查文件是否已正确写入
    with open(output_file) as f:
        lines = f.readlines()
        assert len(lines) == 3  # 包含表头的文件应该有3行
        assert lines[1].strip() == "file1\t0.5000\n"
        assert lines[2].strip() == "file2\t0.2500\n"

# 为process_file_pair函数编写单元测试
def test_process_file_pair(mocker):
    # 使用mock来模拟文件内容
    mocker.patch('Gene_coverage_cache.read_fasta_lengths', return_value={'genome1': 1000, 'cds1': 400})
    
    # 调用函数
    result = process_file_pair('genome.fasta', 'cds.fasta')
    
    # 断言结果是正确的
    assert result == ('genome.fasta', {'cds1': 0.4})

# 为process_file_pairs函数编写单元测试
def test_process_file_pairs(mocker, tmp_path):
    # 使用mock来模拟文件内容和文件对
    mocker.patch('Gene_coverage_cache.read_fasta_lengths', side_effect=[{'genome1': 1000}, {'cds1': 400}])
    mocker.patch('Gene_coverage_cache.process_file_pair', return_value=('genome.fasta', {'cds1': 0.4}))
    
    # 设置测试目录和输出文件
    test_dir = tmp_path / "test_dir"
    test_dir.mkdir()
    genome_file = test_dir / "genome.fasta"
    cds_file = test_dir / "cds.fasta"
    genome_file.touch()
    cds_file.touch()
    
    output_file = tmp_path / "output.txt"
    
    # 调用函数
    process_file_pairs(str(test_dir), str(output_file))
    
    # 检查输出文件是否已正确写入
    with open(output_file) as f:
        lines = f.readlines()
        assert len(lines) == 3  # 包含表头的文件应该有3行
        assert lines[1].strip() == "genome.fasta\t0.4000\n"

# 为main函数编写单元测试
def test_main(mocker):
    # 使用mock来模拟命令行参数和文件存在性
    mocker.patch('argparse.ArgumentParser.parse_args', return_value=argparse.Namespace(directory='test_dir', output_file='output.txt'))
    mocker.patch('os.path.exists', return_value=True)
    mocker.patch('Gene_coverage_cache.process_file_pairs')
    
    # 调用函数
    main()
    
    # 检查是否调用了process_file_pairs函数
    Gene_coverage_cache.process_file_pairs.assert_called_once_with('test_dir', 'output.txt')