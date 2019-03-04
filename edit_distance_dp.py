filename = "chr1.GRCh38.excerpt.fasta"


def readGenome(file):
    genome = ""
    with open(file, 'r') as f:
        for line in f:
            if not line[0] == ">":
                genome += line.rstrip()
    return genome



def edit_distance(p, t):
    """
    Returns the edit distance between pattern, p and 
    text, t using dynamic programming. 
    Example:
    >>> p = "shake spea"
    >>> t = "Shakespear"
    >>> edit_distance(p, t)
    3
    """
    matrix = []
    # Matrix will be a list of lists:
    # [[0][0], [0][1], ...]
    for i in range(len(p) + 1):
        matrix.append([0] * (len(t) + 1))

    # Initializing the first row and column:
    for i in range(len(p) + 1):
        matrix[i][0] = i

    for i in range(len(t) + 1):
        matrix[0][i] = i
    # Fill in the rest of the matrix:
    for i in range(1, len(p) + 1):
        for j in range(1, len(t) + 1):
            dist_hor = matrix[i][j - 1] + 1
            dist_vert = matrix[i - 1][j] + 1

            if p[i - 1] == t[j - 1]:
                dist_diag = matrix[i - 1][j - 1]
            else:
                dist_diag = matrix[i - 1][j - 1] + 1

            matrix[i][j] = min(dist_hor, dist_vert, dist_diag)

    # Return the edit distance between the two strings.
    # This would be the last element in the matrix.
    return matrix[-1][-1]


def edit_distance_approx(p, t):
    """
    Returns the edit distance between pattern, p and 
    text, t using dynamic programming. 
    Example 1:
    >>> p = "GCGTATGC"
    >>> t = "TATTGGCTATACGGTT"
    >>> edit_distance_approx(p, t)
    2

    Example 2:
    >>> p = "GCTGATCGATCGTACG"
    >>> human_chr1 = readGenome("chr1.GRCh38.excerpt.fasta")
    >>> edit_distance_approx(p, human_chr1)
    3

    Example 3:
    >>> p = "GATTTACCAGATTGAG"
    >>> human_chr1 = readGenome("chr1.GRCh38.excerpt.fasta")
    >>> edit_distance_approx(p, human_chr1)
    2

    """
    matrix = []
    # Matrix will be a list of lists:
    # [[0][0], [0][1], ...]
    for i in range(len(p) + 1):
        matrix.append([0] * (len(t) + 1))

    # Initializing the first row and column:
    for i in range(len(p) + 1):
        matrix[i][0] = i

    for i in range(len(t) + 1):
        matrix[0][i] = 0
    # Fill in the rest of the matrix:
    for i in range(1, len(p) + 1):
        for j in range(1, len(t) + 1):
            dist_hor = matrix[i][j - 1] + 1
            dist_vert = matrix[i - 1][j] + 1

            if p[i - 1] == t[j - 1]:
                dist_diag = matrix[i - 1][j - 1]
            else:
                dist_diag = matrix[i - 1][j - 1] + 1

            matrix[i][j] = min(dist_hor, dist_vert, dist_diag)

    # Return the edit distance between the two strings.
    # This would be the last element in the matrix.
    i = len(p)
    last_row = []
    for j in range(len(t) + 1):
        last_row.append(matrix[i][j])

    return min(last_row)


if __name__ == "__main__":
    import doctest
    if doctest.testmod().failed == 0:
        print("Tests passed.")
