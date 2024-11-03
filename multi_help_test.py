# MixingTool 帮助文档模板
import argparse
import sys

def display_intro():
    print("欢迎使用 MixingTool！")
    print("MixingTool 是一个用于各种配液和浓度计算的工具。")
    print("您可以使用以下功能模块：\n")

def display_basic_mixing_help(args):
    print("=== 基础配液 ===")
    print("功能描述：用于计算混合两种或多种溶液的最终浓度。")
    print("命令格式：")
    print(f"   python mixing_tool.py basic {args.concentration1} {args.volume1} {args.concentration2} {args.volume2}")
    print("示例：")
    print("   python mixing_tool.py basic 1.0 100 2.0 50")
    print("这个命令将混合1.0 M的100 mL溶液和2.0 M的50 mL溶液。\n")

def display_dilution_help(args):
    print("=== 稀释配液 ===")
    print("功能描述：用于计算将浓缩溶液稀释到目标浓度所需的稀释剂体积。")
    print("命令格式：")
    print(f"   python mixing_tool.py dilution {args.target_concentration} {args.target_volume} {args.original_concentration}")
    print("示例：")
    print("   python mixing_tool.py dilution 0.5 200 2.0")
    print("这个命令将计算如何将2.0 M的浓缩溶液稀释到0.5 M，最终体积为200 mL。\n")

def display_concentration_calculation_help(args):
    print("=== 浓度计算 ===")
    print("功能描述：用于根据配制的溶液和体积计算其浓度。")
    print("命令格式：")
    print(f"   python mixing_tool.py concentration {args.solute_mass} {args.solution_volume}")
    print("示例：")
    print("   python mixing_tool.py concentration 10 500")
    print("这个命令将计算10克溶质在500 mL溶液中的浓度。\n")

def display_contact_info():
    print("联系方式:")
    print("如果您有任何问题或反馈，请联系：support@example.com")

def main():
    parser = argparse.ArgumentParser(description='MixingTool 帮助文档')
    subparsers = parser.add_subparsers(dest='command', required=True)

    # 基础配液
    basic_parser = subparsers.add_parser('basic', help='基础配液帮助')
    basic_parser.add_argument('concentration1', type=float, help='溶液1的浓度（单位：M）')
    basic_parser.add_argument('volume1', type=float, help='溶液1的体积（单位：mL）')
    basic_parser.add_argument('concentration2', type=float, help='溶液2的浓度（单位：M）')
    basic_parser.add_argument('volume2', type=float, help='溶液2的体积（单位：mL）')

    # 稀释配液
    dilution_parser = subparsers.add_parser('dilution', help='稀释配液帮助')
    dilution_parser.add_argument('target_concentration', type=float, help='目标浓度（单位：M）')
    dilution_parser.add_argument('target_volume', type=float, help='目标体积（单位：mL）')
    dilution_parser.add_argument('original_concentration', type=float, help='原始浓度（单位：M）')

    # 浓度计算
    concentration_parser = subparsers.add_parser('concentration', help='浓度计算帮助')
    concentration_parser.add_argument('solute_mass', type=float, help='溶质质量（单位：克）')
    concentration_parser.add_argument('solution_volume', type=float, help='溶液体积（单位：mL）')

    args = parser.parse_args()

    # 根据命令选择调用相应的函数
    if args.command == 'basic':
        display_basic_mixing_help(args)
    elif args.command == 'dilution':
        display_dilution_help(args)
    elif args.command == 'concentration':
        display_concentration_calculation_help(args)

    display_contact_info()

if __name__ == "__main__":
    main()
