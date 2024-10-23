from Bio import SeqIO
import os
 
# 定义原始序列文件所在目录路径
input_dir = "E:/Euplotes/8_HGT/euk"
 
# 获取该目录下所有的序列文件
file_list = [f for f in os.listdir(input_dir) if f.endswith(".fa")]
 
for file_name in file_list:
    # 构建输入文件的完整路径
    input_file = os.path.join(input_dir, file_name)
    
    # 读取序列记录
    records = list(SeqIO.parse(input_file, "fasta"))
    
    # 遍历每条序列记录并修改其ID（也就是序列名称）
    for record in records:
        new_id = "euk_" + record.description   # 这里将新的序列名称设置为“new_”加上原始序列名称
        
        # 修改序列记录的ID属性
        record.id = new_id
        record.description = ""
        
    # 保存修改后的序列到新的文件中
    output_file = os.path.splitext(os.path.basename(input_file))[0] + "_modified.fasta"
    with open(output_file, 'w') as handle:
        SeqIO.write(records, handle, "fasta")