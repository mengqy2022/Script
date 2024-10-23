#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/03/28
# @Description : TCDB转运蛋白注释

import re
import pandas as pd
import requests
from bs4 import BeautifulSoup
import html
import argparse

def use(diamond_out, listSuperfamilies, output_file):
    data = requests.get(r'https://tcdb.org/browse.php')
    htm = data.content.decode('utf-8')
    ID_name={}
    for term in re.findall(r'&nbsp;([0-9A-Z\.]+):&nbsp;(.*?)</div>',htm):
        ID = term[0]
        name = term[1]
        name = html.unescape(name)
        name = name.replace(r'<sup>+</sup>', '+')
        name = name.replace(r'<sup>2+</sup>', '2+')
        name = name.replace(r'<sup>3+</sup>', '3+')
        ID_name[ID] = name

    alig = pd.read_csv(diamond_out,sep='\t',header=None)
    alig = alig.iloc[:,[0,1]]
    alig.columns = ['GENE','TCID']
    alig['TCID'] = alig['TCID'].map(lambda x : x.split('|')[-1])

    alig_hier = pd.read_csv(listSuperfamilies,sep='\t')
    alig_hier=alig_hier.loc[:,['Family','Fam_abbreviation','Superfamily']]
    alig_hier.drop_duplicates(subset=['Family','Fam_abbreviation','Superfamily'],keep='first',inplace=True)
    alig['class'] = alig['TCID'].map(lambda x : '.'.join(x.split('.')[:1]))
    alig['class_term'] = alig['class'].map(lambda x : ID_name[x])
    alig['subclass'] = alig['TCID'].map(lambda x : '.'.join(x.split('.')[:2]))
    alig['subclass_term'] = alig['subclass'].map(lambda x:ID_name[x])
    alig['family'] = alig['TCID'].map(lambda x : '.'.join(x.split('.')[:3]))
    alig['family_term'] = alig['family'].map(lambda x:ID_name[x])
    alig = alig.rename(columns={'family':'Family'})
    TCDB_anno = pd.merge(alig,alig_hier,on='Family',how='left')
    TCDB_anno.to_csv(output_file,sep='\t',index=0)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='TCDB注释 >>>> mqy <<<<')
    parser.add_argument('--version','-v',action='version',version='%(prog)s 1.0')
    parser.add_argument('--diamond','-d',type=str, help='diamond：默认比对结果文件',required=True)
    parser.add_argument('--listSuperfamilies','-l',type=str, help='listSuperfamilies：TCDB家族列表',required=True)
    parser.add_argument('--output','-o',type=str, default='TCDS_result.csv', help='output：输出文件名',required=False)
    args = parser.parse_args()

    if args.output == None:
        args.output = 'TCDS_result.csv'
    else:
        args.output = args.output
    
    use(diamond_out=args.diamond,listSuperfamilies=args.listSuperfamilies,output_file=args.output)
    print('\nTCDB注释完成，输出文件：', args.output)
    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')

try:
    use()
except IndexError:
    print('\n一般是由于第一次在TCDB下载listSuperfamilies文件时，里面有一个官方小错误，所有文件应该都只有5列，但是有6列的数据，因此机会输出报错，应该更改6列的数据为5列\n')   