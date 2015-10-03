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
            probabilityString = pattern
    return probabilityString


def buildMotifs(profile, Dna):
    motifs = []
    for i in Dna:
        motifs.append(motifFind(i, profile))
    return motifs


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


def RANDOMIZEDMOTIFSEARCH(Dna, k, t):
    """
    presudocode：
    #
    RANDOMIZEDMOTIFSEARCH(Dna, k, t)
        randomly select k-mers Motifs = (Motif1, …, Motift) in each string
            from Dna
        BestMotifs ← Motifs
        while forever
            Profile ← Profile(Motifs)
            Motifs ← Motifs(Profile, Dna)
            if Score(Motifs) < Score(BestMotifs)
                BestMotifs ← Motifs
            else
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
    bestMotifs = motifs
    while True:
        profile = profileMatrix(motifs, k)
        motifs = buildMotifs(profile, Dna)
        if Score(motifs, k) < Score(bestMotifs, k):
            bestMotifs = motifs
        else:
            return bestMotifs


def thousandRun(Dna, k, t):
    bestmotifs = RANDOMIZEDMOTIFSEARCH(Dna, k, t)
    for i in range(300):
        print i
        motifs = RANDOMIZEDMOTIFSEARCH(Dna, k, t)
        if Score(motifs, k) < Score(bestmotifs, k):
            print motifs, Score(motifs, k)
            bestmotifs = motifs
    return bestmotifs


Dna = readDnasegement('read.txt')
k = 15
t = 20
motifs = ['CCA','CCT','CTT','TTG']
print buildMotifs(profileMatrix(motifs,3), Dna)