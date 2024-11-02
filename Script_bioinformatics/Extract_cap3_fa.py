import argparse
import re

def parse_file(file_path):
    contig_sequences = {}
    current_contig_id = None

    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            # 匹配 Contig ID 行
            contig_match = re.search(r'\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*(.*?)\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*', line)
            if contig_match:
                current_contig_id = contig_match.group(1).strip()  # 提取 Contig ID
            
            # 匹配包含 'consensus' 的行
            if 'consensus' in line and current_contig_id:
                consensus_sequence = line.split()[-1]  # 获取最后一个元素，即共识序列
                
                # 如果当前 Contig ID 不在字典中，初始化一个列表
                if current_contig_id not in contig_sequences:
                    contig_sequences[current_contig_id] = []
                
                # 将共识序列添加到当前 Contig ID 的列表中
                contig_sequences[current_contig_id].append(consensus_sequence)

    return contig_sequences

def main():
    parser = argparse.ArgumentParser(description='提取文件中的共识序列')
    parser.add_argument('file', type=str, help='输入的文件路径')
    parser.add_argument('-o', '--output', type=str, help='输出文件路径', default='output.txt')

    args = parser.parse_args()

    contig_sequences = parse_file(args.file)

    # 输出到指定文件，按照FASTA格式输出
    with open(args.output, 'w') as out_file:
        for contig_id, sequences in contig_sequences.items():
            # 处理 Contig ID，替换 '*' 和 空格，并去掉开头的下划线
            formatted_contig_id = contig_id.replace('*', '').replace(' ', '_').lstrip('_')
            out_file.write(f">{formatted_contig_id}\n")  # 输出处理后的 Contig ID，前面加 ">"
            for sequence in sequences:  # 遍历所有共识序列
                out_file.write(f"{sequence}\n")  # 输出对应的共识序列
    
    print(f"结果已写入到 {args.output}")

if __name__ == "__main__":
    main()
