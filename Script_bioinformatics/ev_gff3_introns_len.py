import re
import sys

def usage():
    print('Description: 从EVM整合的.gff3文件中得到内含子的数量和长度')
    print('Usage: ev_gff3_introns_len [EVM.gff3] > [outfile.txt]')

class gene():  #  面向对象gene
    #  第一种形式：def __init__(self)：
    #  这种形式在__init__方法中，只有一个self，指的是实例的本身，但是在方法的类部，包含三个属性:feature_index、threshold和alpha。
    #  它允许定义一个空的结构，当新数据来时，可以直接添加。实例化时，需要实例化之后，再进行赋值。

    #   __init__函数的作用是定义一个函数，该函数会创建所需类的对象。
    #  作用：当我们创建好一个实例对象之后，会自动调用这个方法，来初始化这个对象
    def __init__(self): 
        self.id = None
        self.exonnum = 0
        self.exon = []
        self.start = 0
        self.end = 0
        self.intron = []

    def calcintron(self):  #  定义函数
        self.intron.append((self.start, self.exon[0][0]))  #  append：追加一个元素。   [0][0]访问数组
        for i in range(len(self.exon)-1): self.intron.append((self.exon[i][1], self.exon[i+1][0]))
        self.intron.append((self.exon[-1][1], self.end))
        return None

    def insert(self, tup, index = 0):
        while index < len(self.exon) and tup[0] > self.exon[index][0]: index += 1
        self.exon.insert(index, tup)
        return None

def read(key, line):
    #  re.split 返回一个列表，其中字符串在每次匹配时被拆分。
    lis = re.split('[\f\n\r\t\v]+', each) #  匹配多个分隔符，每一个都拆。
    if len(lis) != 10: 
        return False
    if type(key) == int:  #  如果等于正整数
        return lis[key-1]  #  因为python索引是重0开始。
    elif key == 'ID':
        ID = re.match('ID=[A-Za-z0-9]+', lis[-2])
        return ID.group()   #  返回匹配到的内容

def main():
    lis = []
    f = open(sys.argv[1], 'r')
    for each in f:
        # 读取type
        typ = read(3, each)
        # 读取ID
        ID = read('ID', each)
        if typ == 'gene':
            temp = gene()  #  调动类
            temp.id = ID
            temp.start = int(read(4, each))
            temp.end = int(read(5, each))
            lis.append(temp)  #  想lis列表中添加内容
        elif typ == 'exon':
            temp.exonnum += 1
            # 按序插入外显子序列
            temp.insert((int(read(4, each)), int(read(5, each))))   #  

    for each in lis:
        if each.exonnum >= 2:
            each.calcintron()
            leng = []
            for i in each.intron:
                leng.append(i[1]-i[0])
            print('intron_length{}'.format(leng))
            #print('{}的外显子数量为{}, 内含子长度为{}'.format(each.id, each.exonnum, leng))
        else:
            #pass
            print('%s的外显子数量为%d'%(each.id, each.exonnum))

try:
    main()
except IndexError:
    usage()
