from itertools import permutations
from time import perf_counter_ns


def readFastq(filename):
    sequences = []

    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
    return sequences


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. 
        Example: 
        >>> a = ""
        >>> b = ""
        >>> print(overlap(a, b))

        """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move just past previous match


def overlap_graph(reads, k):
    """

    """
    start = perf_counter_ns()
    k_reads = {}
    overlapgraph = {}
    len_reads = len(reads)
    for i in range(len_reads):
        for j in range((len(reads[i])) - k + 1):
            kmer = reads[i][j:j + k]
            if kmer not in k_reads:
                k_reads[kmer] = set()
            k_reads[kmer].add(reads[i])

    for kmer in k_reads:
        for a, b in permutations(k_reads[kmer], 2):
            olen = overlap(a, b, k)
            if olen > 0:
                overlapgraph[(a, b)] = olen

    count = 0
    seen = set()
    for entry in overlapgraph:
        if entry[0] not in seen:
            seen.add(entry[0])
        else:
            count += 1
    print("count is: ", count)
    # print(len(overlapgraph))
    end = perf_counter_ns() - start
    print(f"Time overlap_graph takes: {end/100} s")



# reads = ['CGTACG', 'TACGTA', 'GTACGT', 'ACGTAC', 'GTACGA', 'TACGAT']
reads = readFastq("ERR266411_1.for_asm.fastq")
# print(len(reads))
overlap_graph(reads, 30)
