import argparse
import re

def process_and_write(input_file, output_file):
    """
    读取文件，将多个空格转换为单个制表符，并去除每行末尾的制表符，
    然后将处理后的内容写入新文件。
    """
    with open(input_file, 'r') as file_read, open(output_file, 'w') as file_write:
        for line in file_read:
            # 处理每行，替换空格为制表符并去除行尾制表符
            processed_line = re.sub(r'\s+', '\t', line.rstrip())
            # 将处理后的行写入输出文件
            file_write.write(processed_line + '\n')

def main():
    # 创建解析器并添加参数
    parser = argparse.ArgumentParser(description="Convert multiple spaces to tabs and write to a new file.")
    parser.add_argument("-i", "--input_file", help="The input file to be processed.", required=True)
    parser.add_argument("-o", "--output_file", help="The output file to write the processed content.", required=True)
    
    # 解析参数
    args = parser.parse_args()
    
    # 调用函数处理并写入文件
    process_and_write(args.input_file, args.output_file)

    print('\nRunning successfully!')
    print(f"\nProcessing complete. Results saved to {args.output_file}\n")
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()