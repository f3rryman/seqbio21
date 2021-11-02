import numpy as np

def is_match(c1, c2, match, mismatch):
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
            score[0] = f[i, j] + is_match(x[i], y[j], match, mismatch)
            score[1] = f[i, j+1] - gap
            score[2] = f[i+1, j] - gap
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


def hirschberg_score(x, y):
    pass


def nw_hirschberg(x, y, match=1, mismatch=-1, gap=1):
    # https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
    seqX = ""
    seqY = ""
    if len(x) == 0:
        for i in range(len(y)):
            seqX += '-'
            seqY += y[i]
    elif len(y) == 0:
        for i in range(len(x)):
            seqX += x[i]
            seqY += '-'
    elif len(x) == 1 or len(y) == 1:
        (seqX, seqY, s) = nw_basic(x, y)
    else:
        xlen = len(x)
        xmid = len(x)/2
        ylen = len(y)

        scoreL = hirschberg_score(x[0:xmid], y)
        scoreR = hirschberg_score(x[xmid+1:xlen], y[::-1])
        ymid = np.argmax(scoreL + scoreR[::-1])

        (seqX, seqY) = nw_hirschberg(x[0:xmid], y[0:ymid]) + nw_hirschberg(x[xmid+1:xlen], y[ymid+1:ylen])
    return (seqX, seqY)


def computeF(i, j):
    if (i == 0) and (j == 0):
        return 0
    elif i == 0:
        return -j*1
    elif j == 0:
        return -i*1
    else:
        return np.max(computeF(i - 1, j - 1) + 1,
                      computeF(i - 1, j) - 1,
                      computeF(i, j - 1) - 1)


def nw_woTable(x, y, match=1, mismatch=-1, gap=1):
    n, m = len(x), len(y)
