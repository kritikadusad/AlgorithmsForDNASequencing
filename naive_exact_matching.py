def reverse_complement(strand):
    """returns complementary dna strand (string)
    Example:
    >>> reverse_complement("ATGC")
    'GCAT'
    """

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }
    complementary_strand = ""
    for base in strand:
        complementary_strand = complement[base] + complementary_strand
    return complementary_strand


def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip()  # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome


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


def naive_with_rc(p, t):
    """ Returns a list of indeces of chars in p and complement_p matching in t
    Example:
    >>> naive_with_rc("ATGC", "AATGCTACGTTATGC")
    [1, 11]

    """
    complement_p = reverse_complement(p)
    occurences = []
    if p != complement_p:
        for i in range(len(t) - len(p) + 1):
            match = True
            for j in range(len(p)):
                if t[i + j] != p[j]:
                    match = False
                    break
            if match:
                occurences.append(i)
        for i in range(len(t) - len(p) + 1):
            match = True
            for j in range(len(p)):
                if t[i + j] != complement_p[j]:
                    match = False
                    break
            if match:
                occurences.append(i)

    else:
        for i in range(len(t) - len(p) + 1):
            match = True
            for j in range(len(p)):
                if t[i + j] != p[j]:
                    match = False
                    break
            if match:
                occurences.append(i)
    return occurences


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
    return occurences


""" Quiz questions: """


filename = "lambda_virus.fa"
virus_genome = readGenome(filename)


# 1. How many times does AGGT or its reverse complement ACCT occur in the lambda virus genome? E.g. if AGGT occurs 10 times and ACCT occurs 12 times, you should report 22.
p = "AGGT"

print(len(naive(p, virus_genome)) + len(naive(reverse_complement(p), virus_genome)))
print("ans1. ", len(naive_with_rc(p, virus_genome)))
# 2. How many times does TTAA or its reverse complement occur in the lambda virus genome? Hint: TTAA and its reverse complement are equal, so remember not to double count.
p = "TTAA"
complement = reverse_complement(p)
print(len(naive(p, virus_genome)))
print("ans.2: ", len(naive_with_rc(p, virus_genome)))


# 3. What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome? E.g. if the leftmost occurrence of ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse complement ACTTAGT is at offset 29, then report 29.
p = "ACTAAGT"
reverse_p = reverse_complement(p)
leftmost_occur = min(naive_with_rc(p, virus_genome))
print("ans.3 ", leftmost_occur)


# 4. What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome?
p = "AGTCGA"
print(naive_with_rc(p, virus_genome)[0])
leftmost_occur = min(naive_with_rc(p, virus_genome))

print("ans.4: ", leftmost_occur)

# 5. As we will discuss, sometimes we would like to find approximate matches for P in T. That is, we want to find occurrences with one or more differences.
# For Questions 5 and 6, make a new version of the naive function called naive_2mm that allows up to 2 mismatches per occurrence. Unlike for the previous questions, do not consider the reverse complement here. We're looking for approximate matches for P itself, not its reverse complement.

p = "TTCAAGCC"
occurrences = naive_2mm(p, virus_genome)

print("ans.5: ", len(occurrences))


# 6. What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?

p = "AGGAGGTT"
leftmost_occur = naive_2mm(p, virus_genome)[0]
print("ans.5: ", leftmost_occur)

# 7. Finally, download and parse the provided FASTQ file containing real DNA sequencing reads derived from a human:
#  This dataset has something wrong with it; one of the sequencing cycles is poor quality.

# Report which sequencing cycle has the problem. Remember that a sequencing cycle corresponds to a particular offset in all the reads. For example, if the leftmost read position seems to have a problem consistently across reads, report 0. If the fourth position from the left has the problem, report 3. Do whatever analysis you think is needed to identify the bad cycle. It might help to review the "Analyzing reads by position" video.
filename = "ERR037900_1.first1000.fastq"
human_sequences, _ = readFastq(filename)


def find_N_by_pos(sequences):
    n = [0] * 100

    for sequence in sequences:
        for i in range(len(sequence)):
            if sequence[i] == "N":
                n[i] += 1

    maximum_freq_N = max(n)
    for i in range(100):
        if n[i] == maximum_freq_N:
            return i


print("ans.6: ", find_N_by_pos(human_sequences))

if __name__ == "__main__":
    import doctest
    if doctest.testmod().failed == 0:
        print("\nAll tests have passed.\n")
