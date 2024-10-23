#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/03/08
# @Description : 批量重命名文件

import os
import sys

def usage():
    print('\nDescription: Batch Rename Files\n')
    print('Usage: python3 File_rename.py [Absolute Path]\n')
    print('Example: python3 File_rename.py /Users/username/Desktop/files\n')
    print('警告：请确保输入的路径是正确的，否则可能导致文件丢失或损坏！\n')
    print('警告：脚本简化文件名称，只保留登入号\n')
    print('    >>>> mqy <<<<\n')
    sys.exit()


def changename(orignname):
      picture=os.listdir(orignname)
      for filename in picture:
          # filename1 = filename.split(".")[0]
          # filename2=re.findall(r"\d+\.?\d*", filename1)[0]+".png"
          # srcpath = os.path.join(orignname,filename)
          # allpath = os.path.join(orignname,filename2)
          # os.rename(srcpath,allpath)
 
         #split("_",2)[1]    “_”表示分隔符 ; 2表示分割次数 ; [1]表示选取第 i 个片段
         filename1=filename.split("_")[0] + "_" + filename.split("_")[1] + ".fna"
         
         #设置旧文件名（就是路径+文件名）
         srcpath=os.path.join(orignname,filename)
         #设置新文件名
         allpath= os.path.join(orignname,filename1)
         os.rename(srcpath, allpath)
 
 
def main():
    if __name__ == '__main__':
        if len(sys.argv) != 2 or sys.argv[1] == "-h" or sys.argv[1] == "--help" or sys.argv[1] == "-help":
            usage()
        orignname=str(sys.argv[1])
        changename(orignname)
        print('Done!\n')
        print('End of run, thanks for using it!')
        print('\t>>> mqy <<<<\n')
        print('\t如果有任何问题请及时联系\n')
        print('\t邮箱：<15877464851@163.com>\n')
        sys.exit()

try:
    main()
except IndexError:
    usage()

