#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com
# @Time      : 2024/05/23 8:57
# @File      : Mod_files_reads_name.py

import os
import unittest
from Mod_files_reads_name import main  # 此处假设你的脚本名为 Mod_files_reads_name.py

class TestMainFunction(unittest.TestCase):

    def setUp(self):
        # 创建临时的测试目录和文件
        self.input_dir = "test_input"
        self.output_dir = "test_output"
        self.prefix = "test_prefix"
        os.makedirs(self.input_dir, exist_ok=True)
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 创建测试的FASTA文件
        with open(os.path.join(self.input_dir, "test1.fa"), 'w') as handle:
            handle.write(">seq1\nATCG\n>seq2\nCGTA")
        
        # 创建不应该被处理的其他文件
        with open(os.path.join(self.input_dir, "not_a_fasta.txt"), 'w') as handle:
            handle.write("This is not a fasta file.")

    def tearDown(self):
        # 清理测试环境
        for filename in os.listdir(self.input_dir):
            os.remove(os.path.join(self.input_dir, filename))
        for filename in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, filename))
        os.rmdir(self.input_dir)
        os.rmdir(self.output_dir)

    def test_main_function(self):
        # 调用待测函数
        main(self.input_dir, self.prefix, self.output_dir)
        
        # 验证输出文件是否存在
        expected_files = ["test1_modified.fasta"]
        for expected_file in expected_files:
            self.assertTrue(os.path.exists(os.path.join(self.output_dir, expected_file)),
                            msg=f"Output file {expected_file} does not exist.")

        # 验证输出文件内容是否正确
        with open(os.path.join(self.output_dir, "test1_modified.fasta"), 'r') as output_handle:
            output_content = output_handle.read()
            self.assertIn(">test_prefix_seq1\nATCG\n", output_content)
            self.assertIn(">test_prefix_seq2\nCGTA\n", output_content)
            self.assertNotIn("not_a_fasta.txt", output_content)

if __name__ == '__main__':
    unittest.main()