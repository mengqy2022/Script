#! /usr/bin/env python
# coding: utf-8

import pyperclip, re, sys, argparse

def read_file(filename):
    with open(filename, 'r') as f:
        return f.read()

def main(text,output_file):
    phoneRegex = re.compile(r'''(
                            (\d{3}|\(\d{3}\))?  # area code
                            (\s|-|\.)?          # separatorS 电话号码分割字符可以是空格（\s）、 短横（-）或句点（.）所以这些部分也应该用管道连接
                            \d{3}               # first 3 digits
                            (\s|-|\.)           # separator
                            \d{4}               # last 4 digits
                            (\s*(ext|x|ext.)\s*\d{2,5})?  # extension
                            )''', re.VERBOSE)

    #  Create email regex pattern

    emailRegex = re.compile(r'''(
                            [a-zA-Z0-9._%+-]+   # username
                            @                   # @ symbol
                            [a-zA-Z0-9.-]+      # domain name
                            (\.[a-zA-Z]{2,4})   # dot-something
                            )''', re.VERBOSE)

    #  Find matches in clipboard text

    matches = []

    for groups in phoneRegex.findall(text):
        phoneNum = '-'.join([groups[1], groups[3], groups[5]])
        if groups[8] != '':
            phoneNum += 'x' + groups[8]
        matches.append(phoneNum)

    for groups in emailRegex.findall(text):
        matches.append(groups[0])

    #  Copy all matches to the clipboard

    if len(matches) > 0:
        output = open( output_file , 'w')
        output.write('\n'.join(matches))
        output.close()
    else:
        print('No phone numbers or email addresses found.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract phone numbers and email addresses from clipboard text or a file.')
    parser.add_argument("--file","-f",required=True, help="The file to extract phone numbers and email addresses from.")
    parser.add_argument("--output","-o",required=False, default='output.txt',help="The output file to write the extracted phone numbers and email addresses to. Default is 'output.txt'.")
    argparse_args = parser.parse_args()

    if argparse_args.file and argparse_args.output:
        main(read_file(argparse_args.file),argparse_args.output)
    else:
        main(read_file(argparse_args.file))

    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')