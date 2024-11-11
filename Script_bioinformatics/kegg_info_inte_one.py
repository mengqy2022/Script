#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com
# @Time      : 2024/11/09

import requests
import json
import sys
import re
import argparse

def print_colored(text, color):
    # 模拟颜色输出，这里可以根据需要实现不同的颜色样式
    color_codes = {
        'purple': '\033[95m',
        'green': '\033[92m',
        'red': '\033[91m',
        'reset': '\033[0m'
    }
    print(f"{color_codes.get(color, '')}{text}{color_codes['reset']}")

def fetch_and_process_kegg_module(output_file):
    url = "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00002&format=json&filedir="
    response = requests.get(url)

    if response.status_code == 200:
        module = response.json()
    else:
        print(f"请求失败，状态码：{response.status_code}")
        return

    with open(output_file, 'w') as outFile:
        outFile.write('moduleID\tdescription\tpathway\tlevel0\tlevel1\tlevel2\n')
        for level0 in module['children']:
            level0_name = level0['name']
            for level1 in level0['children']:
                level1_name = level1['name']
                for level2 in level1['children']:
                    level2_name = level2['name']
                    for level3 in level2['children']:
                        module_info = level3['name']

                        ml = module_info.strip().split('  ')
                        moduleID = ml[0]
                        if re.search(r'\[PATH:', module_info):
                            module_name = ml[1].split('[')[0]
                            relatePath = re.findall(r'\[PATH:.*?]', module_info)[0].split(":")[1].strip(']')
                        else:
                            module_name = ml[1]
                            relatePath = 'NA'
                        
                        outFile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (moduleID, module_name, relatePath, level0_name, level1_name, level2_name))

def fetch_and_process_kegg_pathway(output_file):
    url = "https://www.kegg.jp/kegg-bin/download_htext?htext=br08901&format=json&filedir="
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        pathway = response.json()
    except requests.RequestException as e:
        print(f"请求失败，错误信息：{e}")
        return
    except json.JSONDecodeError:
        print("JSON解析失败")
        return

    with open(output_file, 'w') as outFile:
        outFile.write('mapID\tdescription\tlevel1\tlevel2\n')
        for level0 in pathway.get('children', []):
            level0_name = level0.get('name', 'NA')
            for level1 in level0.get('children', []):
                level1_name = level1.get('name', 'NA')
                for level2 in level1.get('children', []):
                    level2_name = level2.get('name', '')
                    ll = level2_name.strip().split('  ')
                    if len(ll) >= 2:
                        outFile.write('map%s\t%s\t%s\t%s\n' % (ll[0], ll[1], level0_name, level1_name))
                    else:
                        print(f"警告：未能正确解析level2名称：{level2_name}")

def fetch_and_process_kegg_compounds(output_file):
    url = "https://www.genome.jp/kegg-bin/download_htext?htext=br08001&format=json&filedir="

    try:
        response = requests.get(url)
        response.raise_for_status()
        compounds = response.json()
    except requests.RequestException as e:
        print(f"请求失败，错误信息：{e}")
        return
    except json.JSONDecodeError:
        print("JSON解析失败")
        return

    with open(output_file, 'w') as outFile:
        outFile.write('CompoundID\tdescription\tlevel1\tlevel2\tlevel3\n')
        for level0 in compounds.get('children', []):
            level0_name = level0.get('name', 'NA')
            for level1 in level0.get('children', []):
                level1_name = level1.get('name', 'NA')
                for level2 in level1.get('children', []):
                    level2_name = level2.get('name', 'NA')
                    for level3 in level2.get('children', []):
                        cl = level3.get('name', '').strip().split('  ')
                        if len(cl) >= 2:
                            outFile.write('%s\t%s\t%s\t%s\t%s\n' % (cl[0], cl[1], level0_name, level1_name, level2_name))
                        else:
                            print(f"警告：未能正确解析level3名称：{level3.get('name', '')}")

if __name__ == "__main__":

    print_colored("\n                由于kegg数据库更新的速度较快，因此我们需要不断更新其信息。\n", 'purple')
    print_colored("                     此脚本用于获取并输出KEGG模块、通路和化合物信息。\n", 'purple')
    print_colored("                             >>> 注意：五种参数独立使用！ <<<", 'red')
    print_colored("                           >>> 注意：五种参数必须选择一种！ <<<\n", 'red')

    parser = argparse.ArgumentParser(description=print_colored('                              >>> KEGG信息获取<<<', 'green'), epilog=print_colored('usage example：python kegg_info_inte_on.py -o kegg_info.txt [--modules] [--pathways] [--compounds]\n', 'green'))
    parser.add_argument('-o', '--output', required=True, help='指定输出文件的名称')
    parser.add_argument('--modules', action='store_true', help = '获取模块信息')
    parser.add_argument('--pathways', action='store_true', help = '获取通路信息')
    parser.add_argument('--compounds', action='store_true', help = '获取化合物信息')

    args = parser.parse_args()
    
    if not (args.modules or args.pathways or args.compounds):
        print("请至少选择一个数据类型：--modules, --pathways, --compounds")
        sys.exit(1)

    if args.modules:
        fetch_and_process_kegg_module(args.output)
    if args.pathways:
        fetch_and_process_kegg_pathway(args.output)
    if args.compounds:
        fetch_and_process_kegg_compounds(args.output)
