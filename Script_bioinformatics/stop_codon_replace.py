#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys

def replace(infile):
    with open(infile) as fi:
        for line in fi:
            line = line.strip()
            if line.startswith(">"):
                print(line)
            else:   
                if "." in line:
                    line = line.replace(".","C")
                    print(line)
                else:
                    print(line)


if __name__ == "__main__":
    replace(sys.argv[1])
