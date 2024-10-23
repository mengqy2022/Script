# 导入必要的模块和函数
import os
import unittest
from unittest.mock import patch, mock_open
from Gene_coverage_cache import read_fasta_lengths, calculate_cds_proportions, write_results_to_file, process_file_pair, process_file_pairs, main

# 定义一个辅助函数，用于生成模拟的FASTA文件内容
# 定义测试文件内容的辅助函数
def mock_fasta_content(file_name, lengths):
    """
    根据给定的ID和长度生成模拟的FASTA文件内容。

    :param file_name: 文件名，用于描述生成内容的来源
    :param lengths: 一个字典，包含ID和对应的序列长度
    :return: 字符串，表示生成的FASTA文件内容
    """
    content = ""
    for id, length in lengths.items():
        content += f">{id}\n{'A' * length}\n"
    return content

# 定义测试类，继承自unittest.TestCase
# 单元测试类
class TestGeneCoverageCache(unittest.TestCase):
    
    # 在每个测试之前运行，用于设置测试环境
    def setUp(self):
        """
        在执行每个测试之前调用，用于创建临时的FASTA文件。
        """
        # 创建临时的基因组FASTA文件和CDS FASTA文件
        # 创建临时的FASTA文件
        self.genome_file = "test_genome.fasta"
        self.cds_file = "test_cds.fasta"
        with open(self.genome_file, "w") as f:
            f.write(mock_fasta_content("genome", {"genome1": 1000, "genome2": 2000}))
        with open(self.cds_file, "w") as f:
            f.write(mock_fasta_content("cds", {"cds1": 500, "cds2": 1000, "cds3": 1500}))
        
        # 定义预期的结果，用于与实际测试结果进行比较
        self.expected_lengths = {"genome1": 1000, "genome2": 2000}
        self.expected_proportions = {"cds1": 0.5, "cds2": 0.5, "cds3": 0.75}

    # 在每个测试之后运行，用于清理测试环境
    def tearDown(self):
        """
        在执行每个测试之后调用，用于删除临时的FASTA文件。
        """
        # 删除之前创建的临时文件
        # 清理测试文件
        os.remove(self.genome_file)
        os.remove(self.cds_file)

    # 测试read_fasta_lengths函数
    def test_read_fasta_lengths(self):
        """
        测试read_fasta_lengths函数是否能正确读取并解析FASTA文件中的序列长度。
        """
        # 调用函数并获取结果
        lengths = read_fasta_lengths(self.genome_file)
        # 与预期结果进行比较
        self.assertEqual(lengths, self.expected_lengths)

    # 测试calculate_cds_proportions函数
    def test_calculate_cds_proportions(self):
        """
        测试calculate_cds_proportions函数是否能正确计算CDS的比例。
        """
        # 使用预期的基因组长度和实际的CDS长度进行计算
        genome_lengths = self.expected_lengths
        cds_lengths = {"cds1": 500, "cds2": 1000, "cds3": 1500}
        proportions = calculate_cds_proportions(genome_lengths, cds_lengths)
        # 与预期结果进行比较
        self.assertEqual(proportions, self.expected_proportions)

    # 测试write_results_to_file函数
    @patch('builtins.open', new_callable=mock_open)
    def test_write_results_to_file(self, mock_file):
        """
        测试write_results_to_file函数是否能正确将结果写入文件。
        """
        # 准备要写入的结果
        results = {"file1": self.expected_proportions}
        # 调用函数
        write_results_to_file("test_output.txt", results)
        # 检查是否调用了文件打开和写入操作
        # 检查是否正确写入文件
        mock_file.assert_called_once_with("test_output.txt", 'w')
        handle = mock_file()
        handle.write.assert_called_with("File Name\tProportion\nfile1\t0.5000\nfile1\t0.5000\nfile1\t0.7500\n")

    # 测试process_file_pair函数
    def test_process_file_pair(self):
        """
        测试process_file_pair函数是否能正确处理单个文件对并返回预期结果。
        """
        # 使用创建的临时文件进行处理
        genome_file = self.genome_file
        cds_file = self.cds_file
        file_name, proportions = process_file_pair(genome_file, cds_file)
        # 与预期结果进行比较
        self.assertEqual(proportions, self.expected_proportions)
        self.assertEqual(file_name, os.path.basename(genome_file))

    # 测试process_file_pairs函数
    @patch('os.listdir')
    @patch('os.path.exists')
    @patch('concurrent.futures.ThreadPoolExecutor')
    def test_process_file_pairs(self, mock_executor, mock_exists, mock_listdir):
        """
        测试process_file_pairs函数是否能正确处理目录中的多个文件对。
        """
        # 模拟文件列表和文件存在性检查
        # 模拟os.listdir和os.path.exists的返回值
        mock_listdir.return_value = ["genome1.fasta", "genome2.fasta"]
        mock_exists.side_effect = lambda x: True

        # 模拟处理函数的返回值
        mock_process_file_pair = mock_executor.return_value.map.return_value.__iter__.return_value
        mock_process_file_pair.__next__.side_effect = lambda: ("file1", self.expected_proportions)

        # 调用函数
        process_file_pairs("test_dir", "test_output.txt")

        # 检查模拟的处理函数是否被调用
        mock_process_file_pair.__next__.assert_called()
        
# 如果直接运行这个文件，则执行单元测试
if __name__ == '__main__':
    unittest.main()