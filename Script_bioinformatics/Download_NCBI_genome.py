import requests
import os
import time
import argparse

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='从NCBI下载基因组序列')
    
    # 添加参数
    parser.add_argument('-i', '--input', required=True, 
                       help='包含GCA编号的输入文件路径')
    parser.add_argument('-o', '--output', required=True,
                       help='输出文件的目录路径')
    parser.add_argument('-d', '--delay', type=float, default=1.0,
                       help='下载间隔时间（秒），默认为1秒')
    
    # 解析参数
    args = parser.parse_args()
    
    # 从文件中读取GCA编号
    try:
        with open(args.input, 'r') as gca_file:
            gca_numbers = [line.strip() for line in gca_file if line.strip()]
    except FileNotFoundError:
        print(f"错误：找不到输入文件 {args.input}")
        return
    except Exception as e:
        print(f"读取输入文件时出错：{e}")
        return
    
    # 创建输出目录
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    def download_genome(gca_number):
        try:
            # Using NCBI Datasets API
            base_url = f"https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{gca_number}/download"
            params = {"include_annotation_type": "GENOME_FASTA"}
            
            print(f"正在下载 {gca_number}...")
            response = requests.get(base_url, params=params)
            
            if response.status_code == 200:
                output_path = os.path.join(args.output, f"{gca_number}.zip")
                with open(output_path, "wb") as f:
                    f.write(response.content)
                print(f"成功下载 {gca_number}")
                return True
            else:
                print(f"下载 {gca_number} 失败: HTTP {response.status_code}")
                return False
                
        except Exception as e:
            print(f"下载 {gca_number} 时出错: {e}")
            return False
    
    # 下载每个基因组
    success_count = 0
    total_count = len(gca_numbers)
    
    for i, gca_number in enumerate(gca_numbers, 1):
        print(f"\n进度: {i}/{total_count}")
        if download_genome(gca_number):
            success_count += 1
        
        if i < total_count:  # 最后一个不需要等待
            time.sleep(args.delay)
    
    print(f"\n下载完成！成功: {success_count}/{total_count}")

if __name__ == "__main__":
    main()
