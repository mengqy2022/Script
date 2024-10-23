#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mqy
# @Time    : 2024/03/28
# @Description : 合并多个PDF文件

import argparse
import PyPDF2
from pdf2docx import Converter

#  合并pdf文件
def merge_pdf(input_pdf, output_merge_pdf):
    merger = PyPDF2.PdfMerger()
    for path in input_pdf:
        merger.append(path)
    merger.write(output_merge_pdf)
    merger.close()

#  去除PDF文件的空白页
def remove_blank_pages(input_pdf, output_pdf):
    reader = PyPDF2.PdfReader(input_pdf[0])
    writer = PyPDF2.PdfWriter()
 
    # 检查是否存在空白页
    for page in reader.pages:
        text = page.extract_text()
        if text:
            writer.add_page(page)

    with open(output_pdf, 'wb') as file:
        writer.write(file)

#  pdf转换为docx
def pdf_to_docx(input_pdf, output_docx):
    converter = Converter(input_pdf[0])
    converter.convert(output_docx)
    converter.close()

#  主函数
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""[Merge multiple PDF files into one.]
                                                [Transform PDF to DOCX and remove blank pages.]
                                                >>>> mqy <<<<
                                                """)
    parser.add_argument('--version', '-v', action='version', version='%(prog)s 1.0')
    parser.add_argument('--input_files', '-i', nargs='+', required=True, help='Input PDF files to merge')
    parser.add_argument('--ouput_merge_file', '-m', help='Output PDF file after merge')
    parser.add_argument('--output_docx', '-d', help='Output DOCX file after convert PDF to DOCX')
    parser.add_argument('--ouptut_remove_blanks', '-r', help='Output PDF file after remove blank pages')
    args = parser.parse_args()

    if args.input_files and args.ouput_merge_file:
        print('\nMerging PDF files...')
        merge_pdf(input_pdf=args.input_files, output_merge_pdf=args.ouput_merge_file)
        print('\nMerged PDF file saved to', args.ouput_merge_file)

    if args.input_files and args.output_docx:
        print('\nConverting PDF to DOCX...')
        pdf_to_docx(input_pdf=args.input_files, output_docx=args.output_docx)
        print('\nConverted DOCX file saved to', args.output_docx)
    
    if args.ouptut_remove_blanks and args.input_files:
        print('\nRemoving blank pages from the PDF file...')
        remove_blank_pages(input_pdf=args.input_files, output_pdf=args.ouptut_remove_blanks)
        print('\nRemoved blank pages PDF file saved to', args.ouptut_remove_blanks) 

    print('\t>>> mqy <<<<\n')
    print('\t如果有任何问题请及时联系\n')
    print('\t邮箱：<15877464851@163.com>\n')
