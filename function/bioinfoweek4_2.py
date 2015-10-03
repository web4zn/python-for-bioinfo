# coding=utf-8
__author__ = 'web4z'
import random


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


def motifFind(dna, profile):
    k = len(profile)
    minInf = float('-Inf')
    nucletides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(dna) - k + 1):
        probability = 1
        pattern = dna[i:i + k]
        for j in range(len(pattern)):
            probability *= profile[j][nucletides[pattern[j]]] / float(sum(profile[j]))
        if minInf < probability:
            minInf = probability
            probabilitystring = pattern
    return probabilitystring

def profileMatrix(Motifs, k):
    profile = []
    nucletides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(k):
        profile.append([1, 1, 1, 1])
    for i in range(len(Motifs)):
        ss = Motifs[i]
        for j in range(len(ss)):
            profile[j][nucletides[ss[j]]] += 1
    return profile


def findConsensus(motifs, k):
    consensus = ''
    profile = profileMatrix(motifs, k)
    for i in profile:
        position = i.index(max(i))
        consensus += 'ACGT'[position]
    return consensus


def Score(motifs, k):
    consensus = findConsensus(motifs, k)
    minscore = float('Inf')
    distance = 0
    for i in motifs:
        distance += HammingDistance(consensus, i)
    return distance


def GIBBSSAMPLER(Dna, k, t, N):
    """
    presudocode：
    #
    GIBBSSAMPLER(Dna, k, t, N)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string
            from Dna
        BestMotifs ← Motifs
        for j ← 1 to N
            i ← Random(t)
            Profile ← profile matrix constructed from all strings in Motifs
                       except for Motifi
            Motifi ← Profile-randomly generated k-mer in the i-th sequence
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
        return BestMotifs
    :param Dna:
    :param k:
    :param t:
    :return:
    """
    motifs = []
    for i in Dna:
        randomInt = random.randint(0, len(i) - k)
        motifs.append(i[randomInt:randomInt + k])
    # 生成随机motifs
    bestMotifs = motifs[:]
    for j in range(N):
        i = random.randint(0, len(Dna) - 1)
        del motifs[i]
        profile = profileMatrix(motifs, k)
        newk_mer = motifFind(Dna[i],profile)
        motifs.insert(i, newk_mer)
        if Score(motifs, k) < Score(bestMotifs, k):
            print bestMotifs,Score(bestMotifs, k)
            print motifs and Score(motifs, k)
            bestMotifs = motifs[:]
    return bestMotifs


Dna = readDnasegement('read.txt')
k = 15
t = 20
N = 2000
best =  GIBBSSAMPLER(Dna, k, t,N)
for i in range(5):
    print i
    newmotif =  GIBBSSAMPLER(Dna, k, t,N)
    if Score(newmotif,k)<Score(best,k):
        best = newmotif
for i in newmotif:
    print i
