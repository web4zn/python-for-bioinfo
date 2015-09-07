# coding=utf-8
import requests


# ACGT转数字
def symboltonumber(symbol):
    symtonum = {'A':0,'C':1,'G':2,'T':3}
    return symtonum[symbol]

def patterntonumber(pattern):
    if pattern == '':
        return 0
    symbol = pattern[-1]
    prefix = pattern[:-1]
    return 4 * patterntonumber(prefix) + symboltonumber(symbol)

# 数字转ATCG
def quotient(index):
    return index / 4


def remainder(index):
    return index - 4 * (index / 4)


def numbertosymbol(num):
    numtosym = {0:'A',1:'C',2:'G',3:'T'}
    return numtosym[num]

def numbertopattern(index, k):
    if k == 1:
        return numbertosymbol(index)
    prefixIndex = quotient(index)
    r = remainder(index)
    symbol = numbertosymbol(r)
    prefixPattern = numbertopattern(prefixIndex, k - 1)
    return prefixPattern + symbol


# 序列法
def computingFrequence(text, k):
    frequenceArray = range(4 ** k)
    for i in range(4 ** k):
        frequenceArray[i] = 0
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        j = patterntonumber(pattern)
        frequenceArray[j] += 1
    return frequenceArray


# 求互补链
def nucleotideComplement(nucleotide):
    """
    匹配核苷酸
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        return 'A'


def reverseComplement(dnaStr):
    # 迭代太占用资源
    # if dnaStr == '':
    #     return ''
    # nucleotide = dnaStr[-1]
    # prefixdnastr = dnaStr[:-1]
    # return nucleotideComplement(nucleotide) + reverseComplement(prefixdnastr)
    newdna = ''
    for i in range(-1, -len(dnaStr) - 1, -1):
        newdna += nucleotideComplement(dnaStr[i])
    return newdna


# 得出pattern出现的位置
def patterncount(dna, pattern):
    for i in range(len(dna) - len(pattern) + 1):
        if pattern == dna[i:len(pattern) + i]:
            print str(i) + ' ',


"""
#
# file = open('text.txt','r')
# dna = file.read()
# file.close()
# pattern = 'CTTGATCAT'
# # pattern = 'ATGATCAAG'
# patterncount(dna,pattern)
"""

# 初级clumpfinding
"""
def clumpfinding(genome, k, t, l):
    frequentpattern = []
    clump = range(4 ** k)
    for i in range(4 ** k):
        clump[i] = 0
    for i in range(len(genome) - l + 1):
        text = genome[i:i + l]
        frequentarray = computingFrequence(text, k)
        # print frequentarray
        for index in range(4 ** k):
            if frequentarray[index] >= t:
                clump[index] = 1

    for i in range(4 ** k):
        if clump[i] == 1:
            pattern = numbertopattern(i, k)
            frequentpattern.append(pattern)
    return frequentpattern
"""


# 进阶clumpfinding
"""
def betterclumpfinding(genome, k, t, l):
    frequentpattern = []  # 空的模式集
    clump = range(4 ** k)
    for i in range(4 ** k):
        clump[i] = 0  # 初始化所有clump的值
    text = genome[0:l]
    frequencyarray = computingFrequence(text, k)  # 计算一次frequencyarray
    for i in range(4 ** k):
        if frequencyarray[i] >= t:
            clump[i] = 1  # 统计clump，并把大于t的赋值为一

    for i in range(len(genome) - l + 1):
        frequencyarray[patterntonumber(genome[i:i + k])] -= 1  # 前k位的频率减一
        index = patterntonumber(genome[i + l - k + 1:i + l + 1])  # 移动窗口并得出后k位的序号
        frequencyarray[index] += 1  # 后k为的序列增加一次频率
        if frequencyarray[index] >= t:
            clump[index] = 1  # 统计clump

    for i in range(4 ** k):
        if clump[i] == 1:
            pattern = numbertopattern(i, k)
            frequentpattern.append(pattern)  # 把clump=1加入到模式集中

    return frequentpattern

file = open('E-coli.txt', 'r')
genome = file.read()
file.close()
result = betterclumpfinding(genome, 9, 3, 500)
"""


