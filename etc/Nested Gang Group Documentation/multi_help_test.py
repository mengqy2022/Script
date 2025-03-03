# MixingTool 帮助文档模板
import argparse
import sys

def display_intro():
    print(f"\n欢迎使用 MixingTool！\n")
    print(f"MixingTool 是一个用于各种配液和浓度计算的工具。\n")
    print(f"添加\"-h\"查看可以使用的功能模块：\n")

def display_help(command):
    help_text = {
        "basic": (
            "=== 基础配液 ===\n"
            "功能描述：用于计算混合两种或多种溶液的最终浓度。\n"
            "命令格式：python mixing_tool.py basic <浓度1> <体积1> <浓度2> <体积2>\n"
            "示例：python mixing_tool.py basic 1.0 100 2.0 50\n"
        ),
        "dilution": (
            "=== 稀释配液 ===\n"
            "功能描述：用于计算将浓缩溶液稀释到目标浓度所需的稀释剂体积。\n"
            "命令格式：python mixing_tool.py dilution <目标浓度> <目标体积> <原始浓度>\n"
            "示例：python mixing_tool.py dilution 0.5 200 2.0\n"
        ),
        "concentration": (
            "=== 浓度计算 ===\n"
            "功能描述：用于根据配制的溶液和体积计算其浓度。\n"
            "命令格式：python mixing_tool.py concentration <溶质质量> <溶液体积>\n"
            "示例：python mixing_tool.py concentration 10 500\n"
        )
    }
    print(help_text.get(command, "无效的帮助信息。"))

def display_contact_info():
    print("\n联系方式: 15877464851\n")
    print("如果您有任何问题或反馈，请联系：15877464851t@163.com\n")

def main():
    parser = argparse.ArgumentParser(description='MixingTool 帮助文档')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-i', '--intro', action='store_true', help='显示简介')
    parser.add_argument('-c', '--contact', action='store_true', help='显示联系方式')
    
    subparsers = parser.add_subparsers(dest='command')

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

    # 在解析参数之前添加一个默认值，确保至少可以执行一个命令
    try:
        args = parser.parse_args()

        if args.intro:
            display_intro()
        elif args.contact:
            display_contact_info()
        elif args.command:
            display_help(args.command)
        else:
            # 默认显示简介
            display_intro()

    except Exception as e:
        print(f"错误：{str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
