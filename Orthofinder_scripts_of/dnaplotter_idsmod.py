import argparse
from Bio import SeqIO

def replace_cds_with_unsure(gbk_file, ids_file, output_file):
    # 读取需要转换的ID
    with open(ids_file, 'r') as f:
        ids_to_replace = set(line.strip() for line in f)

    # 读取GBK文件并进行修改
    with open(output_file, 'w') as output_handle:
        for record in SeqIO.parse(gbk_file, "genbank"):
            new_features = []
            for feature in record.features:
                if feature.type == 'CDS' and feature.qualifiers.get('locus_tag', [''])[0] in ids_to_replace:
                    feature.type = 'unsure'
                new_features.append(feature)
            record.features = new_features
            SeqIO.write(record, output_handle, "genbank")

def main():
    parser = argparse.ArgumentParser(description='将GBK文件中指定ID的CDS替换为unsure')
    parser.add_argument('-g','--gbk_file', type=str, required=True, help='输入GBK文件路径')
    parser.add_argument('-i','--ids_file', type=str, help='包含需要替换ID的文件路径')
    parser.add_argument('-o','--output_file', type=str, help='输出GBK文件路径')

    args = parser.parse_args()
    replace_cds_with_unsure(args.gbk_file, args.ids_file, args.output_file)

if __name__ == '__main__':
    main()
