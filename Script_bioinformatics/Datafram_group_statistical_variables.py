#!/bin/bash/
# **************************************************
# 分组一个变量，提取另一变量中的观测值
# Date   : 2024-06-25
# Author : 孟庆瑶
# **************************************************

import argparse
import pandas as pd
def read_and_filter_comments(file_path):
    """
    读取文件并过滤掉以 '#' 开头的注释行。
    
    :param file_path: 文件路径
    :return: 去除注释后的Pandas DataFrame
    """
    # 使用列表推导式和pandas的read_csv函数结合，去除以'#'开头的行
    # 这里假设文件的第一行如果不是注释，则为列名
    data = pd.read_csv(file_path, comment='#', header=0, sep="\t")
    return data

def handle_file(df, group_column, agg_column, output_file_path):

    # 按group_column分组，agg_column列的观测值数量统计
    # 结果以字典形式返回，键为group_column的类别，值为agg_column的描述列表
    description_by_category = df.groupby(group_column)[agg_column].apply(list).to_dict()
    
    # 打开文件准备写入
    with open(output_file_path, 'w', encoding='utf-8') as file:
        # 遍历字典，将类别和描述列表转换为字符串并写入文件
        for category, descriptions in description_by_category.items():
            # 将类别与描述列表中的每个元素用制表符连接
            # 注意：此处将描述列表转换为一个由制表符分隔的字符串
            descriptions_str = ",".join(map(str, descriptions))  # 确保所有元素都被转换为字符串
            file.write(f"{category}\t{descriptions_str}\n")  # 每个类别及其描述后跟换行符

def main():
    parser = argparse.ArgumentParser(description="Process data file and generate statistics.")
    parser.add_argument("-i", "--input_file", help="Path to the input data file.", required=True)
    parser.add_argument("-o", "--output_file", help="Path to the output txt file.", required=True)
    parser.add_argument("-g", "--group_column", help="Column to group by.", required=True)
    parser.add_argument("-a", "--agg_column", help="Column to aggregate.", required=True)
    args = parser.parse_args()

    # 读取并清理文件内容
    cleaned_data = read_and_filter_comments(args.input_file)

    # 假设我们要统计'column1'中每个'group_column'的观测值数量，并将结果保存在'new_count'列
    group_column = args.group_column
    agg_column = args.agg_column
    output_file_path = args.output_file
    handle_file(cleaned_data, group_column, agg_column, output_file_path)

    print('\nRunning successfully!')
    print(f"\nProcessing complete. Results saved to {args.output_file}\n")
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

if __name__ == "__main__":
    main()