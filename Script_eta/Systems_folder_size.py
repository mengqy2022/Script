import os
import argparse

def get_folder_size(folder):
    """
    计算指定文件夹的总大小。
    
    :param folder: str, 文件夹路径
    :return: int, 文件夹大小（以字节为单位）
    """
    total_size = 0
    large_files = []  # 存储大文件的列表
    for dirpath, dirnames, filenames in os.walk(folder):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            try:
                file_size = os.path.getsize(file_path)
                total_size += file_size
                if file_size > 5 * 1024 * 1024 * 1024:  # 大于 5GB
                    large_files.append(file_path)  # 加入大文件列表
            except FileNotFoundError:
                print(f'文件未找到: {file_path}')
    return total_size, large_files  # 返回文件夹大小和大文件列表

def main():
    # 设置argparse
    parser = argparse.ArgumentParser(description='查看指定文件夹下每个子文件夹的大小，并孵出大于10GB的文件.')
    parser.add_argument('folder_path', type=str, help='要检查的文件夹路径')

    args = parser.parse_args()

    folder_path = args.folder_path

    if not os.path.isdir(folder_path):
        print("输入的路径不是一个有效的文件夹。")
        return

    # 遍历该文件夹下的每个项目
    for item in os.listdir(folder_path):
        item_path = os.path.join(folder_path, item)

        if os.path.isdir(item_path):  # 如果是子文件夹
            folder_size, large_files = get_folder_size(item_path)
            print(f'文件夹: {item} - 大小: {folder_size / (1024 * 1024):.2f} MB')  # 转换为MB并保留两位小数
            
            # 打印大于10GB的文件
            for large_file in large_files:
                print(f'大于5GB的文件: {large_file} - 大小: {os.path.getsize(large_file) / (1024 * 1024 * 1024):.2f} GB')
        else:
            print(f'跳过非文件夹项目: {item}')

if __name__ == '__main__':
    main()
