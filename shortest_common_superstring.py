import itertools


def overlap(a, b, min_length=3):
    """ Return length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # move just past previous match


def scs(ss):
    """ 
    Returns shortest common superstring of given
    strings, which must be the same length 

    Example 1:
    >>> scs(["ABC", "BCA", "CAB"])
    'ABCAB'

    Example 2:
    >>> len(scs(["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]))
    11
    """
    shortest_sup = None
    all_substrings = set()
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss) - 1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i + 1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i + 1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring

    return shortest_sup  # return shortest


def scs_list(ss):
    """ 
    Returns shortest common superstring of given
    strings, which must be the same length 

    Example 1:
    >>> scs_list(['ABC', 'BCA', 'CAB'])
    3
    >>> scs_list(['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA'])
    10
    >>> scs_list(["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"])
    4
    """
    shortest_sup_1 = scs(ss)
    all_substrings = {shortest_sup_1}
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss) - 1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i + 1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i + 1][olen:]

        if len(sup) == len(shortest_sup_1):
            shortest_sup = sup  # found shorter superstring
            all_substrings.add(shortest_sup)

    return len(all_substrings)


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


reads = readFastq("ads1_week4_reads.fq")


def pick_max_overlap(reads, k):
    """ Returns 2 reads with maximum overlap"""
    reada, readb = None, None
    best_olen = 0
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > best_olen:
            best_olen = olen
            reada, readb = a, b
    return reada, readb, best_olen


def greedy_scs(reads, k):
    """
    Example 1:
    >>> greedy_scs(["ABC", "BCA", "CAB"], 2)
    'CABCA'
    """

    reada, readb, olen = pick_max_overlap(reads, k)
    while olen > 0:
        reads.remove(reada)
        reads.remove(readb)
        reads.append(reada + readb[olen:])
        reada, readb, olen = pick_max_overlap(reads, k)
    return "".join(reads)


print(greedy_scs(reads, 30))

if __name__ == "__main__":
    import doctest
    if doctest.testmod().failed == 0:
        print("Tests passed.")
