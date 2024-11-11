#!/usr/bin/env python 
# -*- coding: utf-8 -*-
# @Author  : mengqy
# @Time    : 2023/10/08

import argparse

def parse_go_file(input_file, output_file):
    with open(input_file, "r") as file:
        lib = {}
        for line in file:
            line = line.strip()
            col_name = line.split(":")[0]
            if col_name == "id":
                id = line.split(" ", maxsplit=1)[1]
                lib[id] = ""
            if col_name == "name":
                name = line.split(" ", maxsplit=1)[1]
                lib[id] = lib[id] + "@" + name
            if col_name == "namespace":
                namespace = line.split(" ", maxsplit=1)[1]
                lib[id] = lib[id] + "@" + namespace
                
    with open(output_file, "a+") as out:
        out.write("Class" + "\t" + "GO_IDs" + "\t" + "Description" + "\n")
        for key in lib.keys():
            go_id = key
            go_name = lib[key].split("@")[1]
            go_namespace = lib[key].split("@")[2]
            if go_namespace == "molecular_function":
                go_namespace = "MF"
                out.write(go_namespace + "\t" + go_id + "\t" + go_name + "\n")
            if go_namespace == "biological_process":
                go_namespace = "BP"
                out.write(go_namespace + "\t" + go_id + "\t" + go_name + "\n")
            if go_namespace == "cellular_component":
                go_namespace = "CC"
                out.write(go_namespace + "\t" + go_id + "\t" + go_name + "\n")

def main():
    parser = argparse.ArgumentParser(description='Process GO files.',
                                    epilog='更详细的信息请访问: https://mengqy2022.github.io/genomics/GO/')
    parser.add_argument('-i', '--input_file', type=str, help='The input GO file to be processed', required=True)
    parser.add_argument('-o', '--output_file', type=str, help='The output file to save results', required=True)

    args = parser.parse_args()
    
    parse_go_file(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
