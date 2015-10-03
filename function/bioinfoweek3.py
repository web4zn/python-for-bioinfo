# coding=utf-8
__author__ = 'web4z'


# 数字转ATCG
def quotient(index):
    return index / 4


def remainder(index):
    return index - 4 * (index / 4)


def numbertosymbol(num):
    numtosym = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    return numtosym[num]


def numbertopattern(index, k):
    if k == 1:
        return numbertosymbol(index)
    prefixIndex = quotient(index)
    r = remainder(index)
    symbol = numbertosymbol(r)
    prefixPattern = numbertopattern(prefixIndex, k - 1)
    return prefixPattern + symbol


def k_merBuild(dna, k):
    """
    得出长度为k的模式，并返回
    :param dna:
    :param k:
    :return:
    """
    k_mer = []
    for i in range(len(dna) - k + 1):
        k_mer.append(dna[i:i + k])
    return k_mer


def Suffix(pattern):
    """
    去掉第一位
    :param pattern:
    :return:
    """
    return pattern[1:]


def HammingDistance(Str1, Str2):
    """
    计算汉明距
    :param Str1:
    :param Str2:
    :return:
    """
    d = 0
    for i in range(len(Str1)):
        if Str1[i] != Str2[i]:
            d += 1
    return d


def neighbors(pattern, d):
    """
    得出汉明距小于d的相似序列
    :param pattern:
    :param d:
    :return:
    """
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ['A', 'T', 'C', 'G']
    Neighborhood = []
    SuffixNeighbors = neighbors(Suffix(pattern), d)
    for i in SuffixNeighbors:
        if HammingDistance(Suffix(pattern), i) < d:
            for j in ['A', 'T', 'C', 'G']:
                Neighborhood.append(j + i)
        else:
            Neighborhood.append(pattern[0] + i)
    return Neighborhood


def unequalHammingDistance(Str1, Str2, d):
    """
    不等长汉民距
    :param Str1:
    :param Str2:
    :param d:
    :return:
    """
    if len(Str1) > len(Str2):
        for i in range(len(Str1) - len(Str2) + 1):
            if HammingDistance(Str2, Str1[i:i + len(Str2)]) <= d:
                return True
    else:
        for i in range(len(Str2) - len(Str1) + 1):
            if HammingDistance(Str1, Str2[i:i + len(Str1)]) <= d:
                return True
    return False


def MOTIFENUMERATION(dna, k, d):
    pattern = []
    for i in dna:
        k_mers = k_merBuild(i, k)
    for i in k_mers:
        k_mer_neighbors = neighbors(i, d)
        for j in k_mer_neighbors:
            count = 0
            for k in dna:
                if unequalHammingDistance(j, k, d):
                    count += 1
            if count == len(dna):
                pattern.append(j)
    pattern = {}.fromkeys(pattern).keys()  # 排序
    return pattern


def DistanceBetweenPatternAndStrings(pattern, dna):
    k = len(pattern)
    distance = 0
    for eachstring in dna:
        hammingDistance = float('Inf')
        for i in range(len(eachstring) - k + 1):
            if hammingDistance > HammingDistance(pattern, eachstring[i:i + k]):
                hammingDistance = HammingDistance(pattern, eachstring[i:i + k])
        distance = distance + hammingDistance
    return distance


def readDnasegement(filename):
    """
    分行获取内容并赋值到一个list
    :param filename:
    :return:
    """
    temp = open(filename, 'r')
    fulldnasegement = temp.read()
    temp.close()
    DNA = []
    k = ''
    for i in fulldnasegement:
        if i == '\n':
            DNA.append(k)
            k = ''
        else:
            k += i
    DNA.append(k)
    return DNA


def MedianString(dna, k):
    """
    找到medianstring
    :param dna:
    :param k:
    :return:
    """
    distance = float('Inf')
    for i in range(4 ** k):
        pattern = numbertopattern(i, k)
        if distance > DistanceBetweenPatternAndStrings(pattern, dna):  # 找到小值
            distance = DistanceBetweenPatternAndStrings(pattern, dna)  # 赋值小值给distance
            median = pattern  # 最后得到最小值对应的pattern即medianstring
    return median


def findK_merInString(dna, k, profile):
    minInf = float('-Inf')
    for i in range(len(dna) - k + 1):
        probability = 1
        pattern = dna[i:i + k]
        for j in range(len(pattern)):
            probability *= profile[pattern[j]][j]
        if minInf < probability:
            minInf = probability
            probabilityString = pattern
    return probabilityString


# ----------------------------------------------------------------------------------------------------------------------
def motifFind(dna, k, profile):
    minInf = float('-Inf')
    nucletides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(dna) - k + 1):
        probability = 1
        pattern = dna[i:i + k]
        for j in range(len(pattern)):
            probability *= profile[j][nucletides[pattern[j]]] / float(sum(profile[j]))
        if minInf < probability:
            minInf = probability
            probabilityString = pattern
    return probabilityString


def findConsensus(motifs, k):
    nucletides = {'A':0,'C':1,'G':2,'T':3}
    profile = []
    consensus =''
    for i in range(k):
        profile.append([0, 0, 0, 0])
    for i in range(len(motifs)):
        motifs_i = motifs[i]
        for j in range(k):
            profile[j][nucletides[motifs_i[j]]]+=1
    for i in profile:
        position = i.index(max(i))
        consensus += 'ACGT'[position]

    return consensus


def Score(motifs, k):
    consensus = findConsensus(motifs,k)
    minscore = float('Inf')
    distance = 0
    for i in motifs:
        distance += HammingDistance(consensus, i)
    return distance


def addingProfileMatrix(profile, motif):
    nucletides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(motif)):
        profile[i][nucletides[motif[i]]] += 1
    return profile


def GREEDYMOTIFSEARCH(dna, k, t):
    bestMotifs = []
    initalString = dna[0]
    otherString = dna[1:]
    for i in dna:
        bestMotifs.append(i[0:k])
    for i in range(len(initalString) - k + 1):
        Motifs = [initalString[i:i + k]]
        firstmotif = Motifs[0]

        profileMatrix = []
        for i in range(k):
            profileMatrix.append([1, 1, 1, 1])

        profileMatrix = addingProfileMatrix(profileMatrix, firstmotif)

        for j in otherString:
            nextmotif = motifFind(j, k, profileMatrix)
            profileMatrix = addingProfileMatrix(profileMatrix, nextmotif)
            Motifs.append(nextmotif)
        if Score(Motifs, k) < Score(bestMotifs, k):
            bestMotifs = Motifs
    return bestMotifs


dna = readDnasegement('read.txt')
motifs = GREEDYMOTIFSEARCH(dna, 7, 3)

for i in motifs:
    print i
# dna = ['GGCGTTCAGGCA','AAGAATCAGTCA','CAAGGAGTTCGC','CACGTCAATCAC','CAATAATATTCG']
# motifs = GREEDYMOTIFSEARCH(dna,3,5)
#
# for i in motifs:
#     print i,
