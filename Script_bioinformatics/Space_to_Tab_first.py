import argparse
import os

def process_file(input_file_path, output_file_path):
    """
    读取输入文件，将每行的第一个空格替换为制表符，其他空格保持不变，
    然后将结果写入输出文件。
    """
    with open(input_file_path, 'r') as input_file, \
         open(output_file_path, 'w') as output_file:
        for line in input_file:
            # 使用字符串的replace方法，只替换第一个空格为制表符
            processed_line = line.replace(' ', '\t', 1)
            output_file.write(processed_line)

def main():
    # 创建解析器
    parser = argparse.ArgumentParser(description="脚本用于将文本文件中每行的第一个空格转换为制表符。")
    
    # 添加参数
    parser.add_argument("-i","--input_file", required=True, help="要处理的输入文本文件路径。")
    parser.add_argument("-o","--output_file", required=True, help="处理后输出的文本文件路径。")
    
    # 解析参数
    args = parser.parse_args()
    
    # 检查输入文件是否存在
    if not os.path.exists(args.input_file):
        parser.error(f"输入文件 {args.input_file} 未找到，请检查路径是否正确。")
    
    # 处理文件并写入新文件
    process_file(args.input_file, args.output_file)

    print('\nRunning successfully!')
    print(f"\nProcessing complete. Results saved to {args.output_file}\n")
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()