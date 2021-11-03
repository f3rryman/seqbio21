import numpy as np

"""Task 01"""
def is_match(c1, c2, match=1, mismatch=-1):
    """Get either match or mismatch value after comparing two chars"""
    if c1 == c2:
        return match
    else:
        return mismatch


def nw_basic(x, y, match=1, mismatch=-1, gap=1):
    """Global alignment via Needleman-Wunsch-algorithm:
        Initialization: F(i,0) = -i * d for all i = 0, 1, 2, ...., n
                        F(j,0) = -j * d for all i = 0, 1, 2, ...., m
        Recursion:                  (F(i-1, j-1) + s(x_i, y_i)
                        F(i,j) = max(F(i-1, j) - d
                                    (F(i, j-1) - d
                        Score is alpha = F(n, m)
        Traceback: From lower right corner to upper left."""
    n = len(x)
    m = len(y)
    # Build score matrix
    f = np.zeros((n+1, m+1))
    # Initialization
    f[:, 0] = np.linspace(0, -n*gap, n+1)
    f[0, :] = np.linspace(0, -m*gap, m+1)
    # Recursion
    score = np.zeros(3)
    for i in range(n):
        for j in range(m):
            score[0] = f[i, j] + is_match(x[i], y[j], match, mismatch) # upper left -> match or mismatch
            score[1] = f[i, j+1] - gap # gap y
            score[2] = f[i+1, j] - gap # gap x
            maximum = np.max(score)
            f[i+1, j+1] = maximum

    # Traceback
    seqX = ""
    seqY = ""
    i = n
    j = m

    while i > 0 and j > 0:
        s = f[i][j]
        sDiag = f[i-1][j-1]
        sUp = f[i-1][j]
        sLeft = f[i][j-1]

        if s == sDiag + is_match(x[i-1], y[j-1], match, mismatch): # match or mismatch case from diagonal
            seqX += x[i-1] # append x[i] to x
            seqY += y[j-1] # append y[i] to y
            i -= 1 # move to diag
            j -= 1
        elif s == sUp - gap: # insert gap from upper
            seqX += x[i-1] # append x[i] to x
            seqY += '-' # append gap to y
            i -= 1 # move to upper cell
        elif s == sLeft - gap: # insert gap from left
            seqX += '-' # append gap to x
            seqY += y[j-1] # append y[j] to y
            j -= 1 # move to left cell
        else: # should never be reached
            print("Something is wrong.")
            break

    # in case of left chars after i or j == 0
    while i > 0:
        seqX += x[i-1]
        seqY += '-'
        i -= 1
    while j > 0:
        seqX += '-'
        seqY += y[j-1]
        j -= 1

    # reverse strings since backtracing is backwards
    seqX = seqX[::-1]
    seqY = seqY[::-1]

    return (seqX, seqY, f[n][m])

"""Task 02"""
def halved_score(x, y, match_score = 1, mismatch_score = -1, gap_penalty = 1):
    """calculates partly scores of the matrix"""
    n, m = len(x), len(y)
    mat = []
    # init matrix
    for i in range(n+1):
        mat.append([0]*(m+1))
    # Initialization
    for j in range(m+1):
        mat[0][j] = -gap_penalty*j
    for i in range(1, n+1):
        mat[i][0] = -gap_penalty*i
        # Recursion
        for j in range(1, m+1):
            if x[n-i] == y[m-j]:
                mat[i][j] = max(mat[i-1][j-1] + match_score,
                            mat[i-1][j] - gap_penalty,
                            mat[i][j-1] - gap_penalty)
            else:
                mat[i][j] = max(mat[i - 1][j - 1] + mismatch_score,
                                mat[i - 1][j] - gap_penalty,
                                mat[i][j - 1] - gap_penalty)
        # Now clear row from memory.
        mat[i-1] = []
    return mat[n]


def nw_linear_space(x, y):
    """Needleman-Wunsch algorithm using only linear space, also called the Hirschberg algorithm"""
    # This is the main Hirschberg routine.
    n, m = len(x), len(y)
    if n<2 or m<2:
        # In this case we just use the N-W algorithm.
        return nw_basic(x, y)
    else:
        # Make partitions, call subroutines.
        nr = int(np.round(n/2, 0))
        F, B = halved_score(x[:nr], y), halved_score(x[nr:], y)
        partition = [F[j] + B[m-j] for j in range(m+1)]
        cut = partition.index(max(partition))
        print(cut)
        # Clear all memory now, so that we don't store data during recursive calls.
        F, B, partition = [], [], []
        # Now make recursive calls.
        call_left = nw_linear_space(x[0:nr], y[0:cut])
        call_right = nw_linear_space(x[nr:], y[cut:])
        # Return result in format: [1st alignment, 2nd alignment, similarity]
        return [call_left[r] + call_right[r] for r in range(3)]

gap_penalty = 1

"""Task 03"""
def computeF (i, j, c1, c2):
    """calculates an entry of F recursively"""
    if i==0 and j==0:
        return 0
    elif i==0:
        return -j*gap_penalty
    elif j==0:
        return -i*gap_penalty
    else:
        return max(computeF(i-1, j-1, c1, c2) + is_match(c1,c2),
                   computeF(i-1, j, c1, c2) - gap_penalty,
                   computeF(i, j-1, c1, c2) - gap_penalty)

def nw_wo_table(x, y):
    """calculates all entries in a matrix for an optimal global alignment recursively
    :returns the score"""
    n, m = len(x), len(y)
    # Initialization
    for i in range(n):
        for j in range(m):
            s = computeF(i, j, x[i], y[j])
    return s