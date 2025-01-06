import pandas as pd
import matplotlib.pyplot as plt
import argparse

def format_dbcan2list(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t', header=None)
    df_filtered = df[[0, 1]]  # 选择需要的字段
    df_filtered.to_csv(output_file, sep='\t', index=False, header=False)

def calculate_abundance(count_file, gene_list_file, column_index, separator, output_file):
    counts = pd.read_csv(count_file, sep=separator, header=None)
    gene_list = pd.read_csv(gene_list_file, sep='\t', header=None)

    gene_counts = {gene: counts.iloc[:, column_index].sum() for gene in gene_list[0]}
    
    total_counts = sum(gene_counts.values())
    tpm = {gene: (count / total_counts) * 1e6 for gene, count in gene_counts.items()}

    abundance_df = pd.DataFrame(list(tpm.items()), columns=['Gene', 'TPM'])
    abundance_df.to_csv(output_file, sep='\t', index=False)

def plot_abundance(tpm_file):
    tpm_data = pd.read_csv(tpm_file, sep='\t')
    plt.bar(tpm_data['Gene'], tpm_data['TPM'])
    plt.xlabel('基因')
    plt.ylabel('丰度 (TPM)')
    plt.title('基因丰度分布')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='处理Diamond比对结果并计算基因丰度')
    
    parser.add_argument('-d', '--diamond_output', required=True, help='Diamond比对结果文件 (.f6)')
    parser.add_argument('-c', '--count_file', required=True, help='基因计数文件')
    parser.add_argument('-g', '--gene_list_output', required=True, help='输出的基因列表文件')
    parser.add_argument('-o', '--tpm_output', required=True, help='输出的TPM文件')
    parser.add_argument('-ci', '--column_index', type=int, required=True, help='计数文件中的列索引 (1-based)')
    parser.add_argument('-s', '--separator', default=',', help='计数文件的分隔符，默认为逗号')

    args = parser.parse_args()

    # 1. 预处理Diamond比对结果
    format_dbcan2list(args.diamond_output, args.gene_list_output)

    # 2. 计算丰度
    calculate_abundance(args.count_file, args.gene_list_output, args.column_index - 1, args.separator, args.tpm_output)

    # 3. 绘制丰度分布图
    plot_abundance(args.tpm_output)
