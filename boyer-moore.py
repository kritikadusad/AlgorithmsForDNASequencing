from bm_preproc import BoyerMoore


filename = "chr1.GRCh38.excerpt.fasta"


def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p) - 1, -1, -1):
            break
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences


def boyer_moore_with_counts(p, p_bm, t):
    """ 
    Returns occurences, # alignments done and # characters matched
    Example 1:
    >>> p = 'word'
    >>> t = 'there would have been a time for such a word'
    >>> lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
    >>> p_bm = BoyerMoore(p, lowercase_alphabet)
    >>> occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
    >>> print(occurrences, num_alignments, num_character_comparisons)
    [40] 12 15

    Example 2:
    >>> p = 'needle'
    >>> t = 'needle need noodle needle'
    >>> lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
    >>> p_bm = BoyerMoore(p, lowercase_alphabet)
    >>> occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
    >>> print(occurrences, num_alignments, num_character_comparisons)
    [0, 19] 5 18

    Example 3: 
    >>> p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    >>> t = readGenome("chr1.GRCh38.excerpt.fasta")
    >>> p_bm = BoyerMoore(p)
    >>> occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
    >>> print(num_alignments)
    127974
    >>> print(len(p))
    """

    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False

        for j in range(len(p) - 1, -1, -1):
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break

        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)

        i += shift
        num_alignments += 1
        num_character_comparisons += (len(p) - j)

    return occurrences, num_alignments, num_character_comparisons


def naive(p, t):
    """ Returns a list of indeces of chars in p matching in t.
    Example:
    >>> naive("ATGC", "AATGCTTTATGC")
    [1, 8]

    """
    occurences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i + j] != p[j]:
                match = False
                break
        if match:
            occurences.append(i)
    return occurences


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


def naive_with_counts(p, t):
    """
    Returns occurences, # alignments done and # characters matched
    Example 1:
    >>> p = 'word'
    >>> t = 'there would have been a time for such a word'
    >>> occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
    >>> print(occurrences, num_alignments, num_character_comparisons)
    [40] 41 46

    Example 2:
    >>> p = 'needle'
    >>> t = 'needle need noodle needle'
    >>> occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
    >>> print(occurrences, num_alignments, num_character_comparisons)
    [0, 19] 20 35

    Example 3:
    >>> p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
    >>> t = readGenome("chr1.GRCh38.excerpt.fasta")
    >>> occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
    >>> print(num_alignments, num_character_comparisons)
    799954 984143
    """
    occurences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i + j] != p[j]:
                match = False
                break
        if match:
            occurences.append(i)
        num_alignments += 1
        num_character_comparisons += (j + 1)
    return occurences, num_alignments, num_character_comparisons


if __name__ == "__main__":
    import doctest
    if doctest.testmod().failed == 0:
        print("Tests passed.")
