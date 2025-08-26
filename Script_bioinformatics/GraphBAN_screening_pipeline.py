import argparse
import torch
import numpy as np
import pandas as pd
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from transformers import AutoTokenizer, AutoModel
from tqdm import tqdm
import os

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='GraphBAN 多序列生物标志物药物筛选工具')
    
    parser.add_argument('--fasta_file', type=str, required=True,
                      help='包含生物标志物序列的FASTA文件路径')
    parser.add_argument('--compound_lib', type=str, required=True,
                      help='化合物SMILES数据集路径(ZINC-250K)')
    parser.add_argument('--model_path', type=str, required=True,
                      help='预训练GraphBAN模型路径')
    parser.add_argument('--threshold', type=float, default=0.5,
                      help='相互作用概率阈值 (默认: 0.5)')
    parser.add_argument('--batch_size', type=int, default=1000,
                      help='化合物批处理大小 (默认: 1000)')
    parser.add_argument('--output_dir', type=str, default='results',
                      help='输出目录 (默认: results)')
    
    return parser.parse_args()

def load_fasta_sequences(fasta_path):
    """加载FASTA文件中的所有序列"""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # 取header第一部分作为ID
                current_seq = []
            else:
                current_seq.append(line)
    
        if current_header is not None:
            sequences[current_header] = ''.join(current_seq)
    
    if not sequences:
        raise ValueError("FASTA文件中未找到有效序列")
    
    return sequences

def load_compound_library(smiles_path):
    """加载化合物库"""
    if not os.path.exists(smiles_path):
        raise FileNotFoundError(f"化合物库文件不存在: {smiles_path}")
    
    return pd.read_csv(smiles_path)

def extract_protein_features(sequence, tokenizer, model):
    """提取蛋白质序列特征"""
    inputs = tokenizer(sequence, return_tensors="pt", 
                      padding=True, truncation=True, 
                      max_length=1024)
    with torch.no_grad():
        outputs = model(**inputs)
    return outputs.last_hidden_state.mean(dim=1)

def extract_compound_features(smiles_list):
    """提取化合物特征（简化示例）"""
    features = []
    for smiles in smiles_list:
        # 实际应用中应使用更复杂的分子特征
        features.append([
            len(smiles),                    # SMILES长度
            len(set(smiles)),               # 独特字符数
            smiles.count('C'),              # 碳原子数
            smiles.count('N'),              # 氮原子数
            smiles.count('O'),              # 氧原子数
            smiles.count('='),              # 双键数
            smiles.count('#'),              # 三键数
            smiles.count('@') + smiles.count('/') + smiles.count('\\')  # 手性中心
        ])
    return torch.tensor(features, dtype=torch.float)

def create_graph_data(protein_feat, compound_feats):
    """创建图数据"""
    num_compounds = compound_feats.shape[0]
    
    # 蛋白质节点特征（复制与化合物数量相同次数）
    protein_nodes = protein_feat.repeat(num_compounds, 1)
    
    # 创建全连接的蛋白质-化合物图
    edge_index = []
    for i in range(num_compounds):
        edge_index.append([0, i+1])  # 蛋白质(0)连接到每个化合物
        edge_index.append([i+1, 0])  # 双向连接
    
    # 合并所有节点特征
    x = torch.cat([protein_feat, compound_feats], dim=0)
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    
    return Data(x=x, edge_index=edge_index)

def predict_interactions(model, data_loader, device):
    """批量预测相互作用"""
    predictions = []
    with torch.no_grad():
        for data in data_loader:
            data = data.to(device)
            out = model(data)
            preds = out.sigmoid().cpu().numpy()
            predictions.extend(preds)
    return np.array(predictions).flatten()

def process_sequence(header, sequence, compound_df, model, tokenizer, protein_model, args, device):
    """处理单个蛋白序列"""
    print(f"\n处理序列: {header} (长度: {len(sequence)} aa)")
    
    # 1. 提取蛋白质特征
    protein_feat = extract_protein_features(sequence, tokenizer, protein_model)
    
    # 2. 分批处理化合物库
    results = []
    for i in tqdm(range(0, len(compound_df), args.batch_size), desc="处理化合物"):
        batch = compound_df.iloc[i:i+args.batch_size]
        
        # 3. 提取化合物特征
        compound_feats = extract_compound_features(batch['smiles'].tolist())
        
        # 4. 构建图数据
        graph_data = create_graph_data(protein_feat, compound_feats)
        data_loader = DataLoader([graph_data], batch_size=1)
        
        # 5. 预测相互作用
        preds = predict_interactions(model, data_loader, device)
        
        # 6. 保存结果
        batch_results = batch.copy()
        batch_results['interaction_prob'] = preds
        batch_results['protein_id'] = header
        results.append(batch_results[preds > args.threshold])
    
    if not results:
        print(f"未找到相互作用概率 > {args.threshold}的化合物")
        return None
    
    return pd.concat(results).sort_values('interaction_prob', ascending=False)

def main():
    args = parse_arguments()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 1. 加载数据
    print(f"\n加载FASTA文件: {args.fasta_file}")
    protein_sequences = load_fasta_sequences(args.fasta_file)
    
    print(f"加载化合物库: {args.compound_lib}")
    compound_df = load_compound_library(args.compound_lib)
    
    # 2. 初始化特征提取器
    print("初始化蛋白质特征提取器...")
    protein_tokenizer = AutoTokenizer.from_pretrained("Rostlab/prot_bert")
    protein_model = AutoModel.from_pretrained("Rostlab/prot_bert")
    
    # 3. 加载GraphBAN模型
    print(f"加载GraphBAN模型: {args.model_path}")
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = torch.load(args.model_path, map_location=device)
    model.eval()
    
    # 4. 处理每条序列
    all_results = []
    for header, sequence in protein_sequences.items():
        result = process_sequence(
            header, sequence, compound_df, 
            model, protein_tokenizer, protein_model, 
            args, device
        )
        
        if result is not None:
            # 保存单个蛋白的结果
            output_file = os.path.join(args.output_dir, f"{header}_results.csv")
            result.to_csv(output_file, index=False)
            print(f"结果已保存到: {output_file}")
            all_results.append(result)
    
    # 5. 保存合并结果
    if all_results:
        combined_results = pd.concat(all_results)
        combined_file = os.path.join(args.output_dir, "combined_results.csv")
        combined_results.to_csv(combined_file, index=False)
        print(f"\n所有结果已合并保存到: {combined_file}")
    
    print("\n药物筛选完成!")

if __name__ == "__main__":
    main()
