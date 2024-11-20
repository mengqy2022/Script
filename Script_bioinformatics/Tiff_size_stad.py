#!/usr/bin/env python3
import os
import argparse
from PIL import Image

def resize_images(input_folder, output_folder, target_size):
    """
    调整指定文件夹中所有 TIFF 格式图像的大小，并保存到输出文件夹。
    
    :param input_folder: 输入文件夹路径
    :param output_folder: 输出文件夹路径
    :param target_size: 目标图像大小（宽度, 高度）
    """
    # 确保输出文件夹存在
    os.makedirs(output_folder, exist_ok=True)

    # 遍历输入文件夹中的所有文件
    for filename in os.listdir(input_folder):
        # 确保文件是 TIFF 图像
        if filename.lower().endswith(('.tif', '.tiff')):
            input_path = os.path.join(input_folder, filename)
            output_path = os.path.join(output_folder, filename)

            try:
                # 打开图像并调整大小
                with Image.open(input_path) as img:
                    img_resized = img.resize(target_size, Image.LANCZOS)
                    img_resized.save(output_path, format='TIFF')
                    print(f"成功处理文件: {filename}")
            except FileNotFoundError:
                print(f"文件未找到: {input_path}")
            except IsADirectoryError:
                print(f"遇到目录而非文件: {input_path}")
            except Exception as e:
                print(f"处理文件 {filename} 时出错: {e}")

    print("所有文件处理完成。")

def main():
    parser = argparse.ArgumentParser(description='调整图像大小的脚本')
    parser.add_argument('-i', '--input_folder', required=True, type=str, help='输入文件夹路径')
    parser.add_argument('-o', '--output_folder', required=True, type=str, help='输出文件夹路径')
    parser.add_argument('-wi', '--width', required=True, type=int, help='目标图像宽度')
    parser.add_argument('-he', '--height', required=True, type=int, help='目标图像高度')

    args = parser.parse_args()
    
    # 确保宽度和高度有效
    if args.width <= 0 or args.height <= 0:
        print("宽度和高度必须为正数。")
        return

    target_size = (args.width, args.height)
    resize_images(args.input_folder, args.output_folder, target_size)

if __name__ == '__main__':
    main()
