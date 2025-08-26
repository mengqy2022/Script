import pandas as pd
import warnings
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Lipinski, QED
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds, CalcNumHBD, CalcNumHBA
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import os

# 抑制RDKit的警告和特定的boost警告
RDLogger.DisableLog('rdApp.*')
warnings.filterwarnings("ignore", 
    category=RuntimeWarning, 
    message="to-Python converter for class boost::shared_ptr<class RDKit::FilterHierarchyMatcher> already registered; second conversion method ignored.")

def standardize_smiles(smiles):
    """标准化SMILES（保留立体化学信息）"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)  # 保持立体化学信息
    return None

def calculate_molecular_properties(smiles):
    """计算分子属性和ADMET相关特征，包括QED评分"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # 基本分子属性
    props = {
        'MolecularWeight': Descriptors.MolWt(mol), # 分子量
        'LogP': MolLogP(mol), # 脂水分配系数
        'NumHAcceptors': CalcNumHBA(mol), # 氢键受体数
        'NumHDonors': CalcNumHBD(mol), # 氢键供体数
        'PolarSurfaceArea': Descriptors.TPSA(mol), # 极性表面积
        'QED': QED.qed(mol),  # QED 定量估计类药性
        'RotatableBonds': CalcNumRotatableBonds(mol), # 可旋转键数
        'AromaticRings': Lipinski.NumAromaticRings(mol), # 芳香环数
        'HeavyAtoms': Lipinski.HeavyAtomCount(mol), # 重原子数
        'FractionCSP3': Lipinski.FractionCSP3(mol), # sp³杂化碳比例
        'SMILES': standardize_smiles(smiles)
    }
    
    # 计算Lipinski规则违反次数
    violations = 0
    if props['MolecularWeight'] > 500: violations += 1
    if props['LogP'] > 5: violations += 1
    if props['NumHAcceptors'] > 10: violations += 1
    if props['NumHDonors'] > 5: violations += 1
    props['LipinskiViolations'] = violations
    
    return props

def filter_pains(smiles_list):
    """PAINS过滤器"""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    
    # 使用tqdm显示进度
    filtered = []
    for smiles in tqdm(smiles_list, desc="Filtering PAINS", unit="mol"):
        mol = Chem.MolFromSmiles(smiles)
        if mol and not catalog.HasMatch(mol):
            filtered.append(smiles)
    return filtered

def process_smiles_chunk(smiles_chunk):
    """处理SMILES块的辅助函数，用于多进程"""
    return [calculate_molecular_properties(smiles) for smiles in smiles_chunk]

def process_file(input_file, top_n=10):
    print(f"\nProcessing {input_file}...")
    
    # 读取数据并去重（基于标准化后的SMILES）
    df = pd.read_csv(input_file)
    df['StandardizedSMILES'] = df['SMILES'].apply(standardize_smiles)
    df = df.drop_duplicates(subset=['StandardizedSMILES'])
    print(f"Initial compounds: {len(df)} (after deduplication)")
    
    # 使用多核计算分子属性
    smiles_list = df['SMILES'].tolist()
    chunk_size = max(1, len(smiles_list) // (cpu_count() * 2))  # 确保至少为1
    
    # 使用tqdm显示进度
    with Pool() as pool:
        results = list(tqdm(
            pool.imap(process_smiles_chunk, 
                     [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]),
            total=len(smiles_list)//chunk_size + 1,
            desc="Calculating properties",
            unit="chunk"
        ))
    
    # 合并结果并再次去重（双重保险）
    properties = [prop for chunk in results for prop in chunk if prop is not None]
    props_df = pd.DataFrame(properties).drop_duplicates(subset=['SMILES'])
    print(f"\nUnique valid SMILES parsed: {len(props_df)}")
    
    # PAINS过滤
    pains_filtered = filter_pains(props_df['SMILES'].tolist())
    pains_df = props_df[props_df['SMILES'].isin(pains_filtered)]
    print(f"After PAINS filter: {len(pains_df)}")
    
    # 按QED评分排序并取前N个
    top_compounds = pains_df.sort_values('QED', ascending=False).head(top_n)
    
    # 保存结果
    output_file = input_file.replace('.csv', '_QED_top{}.csv'.format(top_n))
    top_compounds.to_csv(output_file, index=False)
    
    print(f"\nTop {top_n} compounds by QED score:")
    print(top_compounds[['SMILES', 'QED', 'MolecularWeight', 'LogP', 
                        'NumHAcceptors', 'NumHDonors', 'LipinskiViolations']])
    print(f"\nResults saved to {output_file}")
    
    return top_compounds

if __name__ == '__main__':
    # 处理所有文件
    files = [
        #"./TRPM4.csv",
        "./SLC12A1.csv",
        #"./ITPR1.csv"
        "./GRIA3.csv"
    ]
    
    for file in files:
        if not os.path.exists(file):
            print(f"Warning: File {file} not found, skipping...")
            continue
        top_10 = process_file(file)
    
    print("\nAll files processed successfully!")
