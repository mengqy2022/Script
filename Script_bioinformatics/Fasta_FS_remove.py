import argparse
from Bio import SeqIO
import pandas as pd
import os
import sys

def read_input_file(file_path):
    """读取输入文件（Excel或TSV）"""
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()
    
    try:
        if ext in ('.xlsx', '.xls'):
            return pd.read_excel(file_path, dtype={'FS_type': str, 'FS': int})
        elif ext in ('.tsv', '.txt', '.tab'):
            return pd.read_csv(file_path, sep='\t', dtype={'FS_type': str, 'FS': int})
        else:
            raise ValueError(f"不支持的文件格式: {ext}. 请使用.xlsx, .xls, .tsv或.txt")
    except Exception as e:
        print(f"读取输入文件时出错: {str(e)}")
        sys.exit(1)

def apply_frameshift(sequence, operation, position):
    """应用单个frameshift操作到序列上"""
    if operation == '+1':
        return sequence[:position-1] + sequence[position:]
    elif operation == '+2':
        return sequence[:position-2] + sequence[position:]
    elif operation == '-1':
        return sequence[:position] + sequence[position-1:]
    elif operation == '-2':
        return sequence[:position] + sequence[position-2:]
    else:
        raise ValueError(f"未知操作类型: {operation}")

def adjust_position(current_position, operation, fs_position, fs_operation):
    """
    根据前一个操作调整当前操作的位置
    """
    if fs_operation == '+1' and current_position > fs_position:
        return current_position - 1
    elif fs_operation == '+2' and current_position > fs_position:
        return current_position - 2
    elif fs_operation == '-1' and current_position > fs_position:
        return current_position - 1
    elif fs_operation == '-2' and current_position > fs_position:
        return current_position - 2
    return current_position

def process_sequences(fasta_file, input_file, output_file):
    # 读取FASTA文件
    try:
        sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    except Exception as e:
        print(f"读取FASTA文件时出错: {str(e)}")
        sys.exit(1)
    
    # 读取输入文件并按Gene_id分组
    try:
        df = read_input_file(input_file)
        required_columns = ['Gene_id', 'FS_type', 'FS']
        if not all(col in df.columns for col in required_columns):
            raise ValueError(f"输入文件必须包含以下列: {', '.join(required_columns)}")
        
        # 按Gene_id分组并排序（按FS位置升序）
        grouped = df.groupby('Gene_id', group_keys=False).apply(
            lambda x: x.sort_values('FS')
        )
    except Exception as e:
        print(f"处理输入文件时出错: {str(e)}")
        sys.exit(1)
    
    # 处理每个基因的所有修改
    results = []
    processed_genes = set()
    
    for gene_id, group in grouped.groupby('Gene_id'):
        if gene_id not in sequences:
            print(f"警告: 序列ID {gene_id} 在FASTA文件中未找到")
            continue
            
        seq_record = sequences[gene_id]
        current_seq = str(seq_record.seq)
        operations = []
        previous_operations = [] 
        
        # 按顺序应用所有修改
        for _, row in group.iterrows():
            operation = row['FS_type']
            original_position = row['FS']
            
            # 调整当前位置，考虑之前的操作
            adjusted_position = original_position
            for prev_op, prev_pos in previous_operations:
                adjusted_position = adjust_position(adjusted_position, operation, prev_pos, prev_op)
            
            try:
                current_seq = apply_frameshift(current_seq, operation, adjusted_position)
                operations.append(f"{operation}_{original_position}(adjusted_to_{adjusted_position})")
                previous_operations.append((operation, original_position))
            except IndexError:
                print(f"警告: 基因{gene_id}的位置{adjusted_position}(原位置{original_position})超出序列范围(长度{len(current_seq)})")
                continue
            except ValueError as e:
                print(f"警告: {str(e)} 对于基因{gene_id}")
                continue
        
        # 保存结果
        if operations:
            ops_str = ",".join(operations)
            results.append(f">{gene_id}_modified_{ops_str}\n{current_seq}\n")
            processed_genes.add(gene_id)
    
    # 写入输出文件
    try:
        with open(output_file, 'w') as f:
            f.writelines(results)
        print(f"成功处理了 {len(processed_genes)} 个基因的序列")
    except Exception as e:
        print(f"写入输出文件时出错: {str(e)}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='根据输入文件中的指令处理FASTA序列',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''输入文件要求:
  支持格式: Excel(.xlsx, .xls) 或 制表符分隔文件(.tsv, .txt)
  必须包含以下列:
    - Gene_id: 与FASTA文件中的序列ID匹配
    - FS_type: 操作类型 ("+1", "+2", "-1", "-2")
    - FS: 位置坐标 (整数)
  
示例文件结构:
  Gene_id  FS_type  FS  Sequence
  gene1    +1       10  ATGC...
  gene1    -2       20  CGTA...
  gene2    +1       15  GCTA...'''
    )
    
    parser.add_argument('-f', '--fasta', required=True, help='输入FASTA文件路径')
    parser.add_argument('-i', '--input', required=True, 
                       help='输入指令文件路径(支持.xlsx, .xls, .tsv, .txt)')
    parser.add_argument('-o', '--output', required=True, help='输出FASTA文件路径')
    
    args = parser.parse_args()
    
    process_sequences(args.fasta, args.input, args.output)
    print(f"处理完成，结果已保存到 {args.output}")

if __name__ == '__main__':
    main()