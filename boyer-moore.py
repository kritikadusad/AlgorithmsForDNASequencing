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
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p) - 1, -1, -1):
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                num_character_comparisons.append(i)
                break

        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)

        i += shift
        num_alignments += 1

    return occurrences, num_alignments, len(num_character_comparisons)


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


p = "word"
t = 'there would have been a time for such a word'
lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
p_bm = BoyerMoore(p, lowercase_alphabet)
occurrences, num_alignments, num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
print(occurrences, num_alignments, num_character_comparisons)
