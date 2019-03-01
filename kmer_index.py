#!/usr/bin/env python

"""kmer_index.py: A k-mer index for indexing a text."""

__author__ = "Ben Langmead"

import bisect
from bm_preproc import BoyerMoore


class Index(object):
    """ Holds a substring index for a text T """

    def __init__(self, t, k):
        """ Create index from all substrings of t of length k """
        self.k = k  # k-mer length (k)
        self.index = []
        for i in range(len(t) - k + 1):  # for each k-mer
            self.index.append((t[i:i + k], i))  # add (k-mer, offset) pair
        self.index.sort()  # alphabetize by k-mer

    def query(self, p):
        """ Return index hits for first k-mer of p """
        kmer = p[:self.k]  # query with first k-mer
        i = bisect.bisect_left(self.index, (kmer, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


t = readGenome("chr1.GRCh38.excerpt.fasta")
index = Index(t, 8)


def queryIndex(p, t, index):
    k = index.k
    offsets = []

    for i in index.query(p):
        if p[k:] == t[i + k:i + len(p)]:
            offsets.append(i)
    return offsets


def pigeon_hole_index_matching(p, t, index, n):
    """ 
    n = maximum number of mismatches allowed
    
    """
    segment_length = int(len(p) / (n + 1))
    matches_set = set()

    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))

        occurrences = queryIndex(p[start:end], t, index)

        for o in occurrences:
            if o < start or o - start + len(p) > len(t):
                continue

            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[o - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[o - start + j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                matches_set.add(o - start)
        return len(list(matches_set))


def naive_2mm(p, t):
    occurences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        error = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                error += 1
                if error > 2:
                    match = False
                    break
        if match:
            occurences.append(i)
    return len(occurences)


p = 'GGCGCGGTGGCTCACGCCTGTAAT'

print(pigeon_hole_index_matching(p, t, index, n=2))
print(naive_2mm(p, t))
