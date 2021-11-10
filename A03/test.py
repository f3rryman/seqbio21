def calc_2constraints(s1, s2, lidx, l2idx, constraints):
    """
    calcs indexes of crossovers between two sequences and prints them
    :param s1: seq 1
    :param s2: seq 2
    :return: None, prints indexes to console
    """
    count = 0
    n, m = len(s1), len(s2)
    for i in range(n):  # upper left
        for j in range(1, m):  # lower right
            for l in range(i + 1, n):  # upper right
                for k in range(m):  # lower left
                    if k < j:
                        constraints.append("x{}{}_{}{} + x{}{}_{}{} < 1;"
                                      .format(lidx, i, l2idx, j, lidx, l, l2idx, k))
                        count += 1
    return (constraints, count)


def calc_2constraints(s1, s2, lidx, l2idx, constraints):
    """
    calcs constraints between two sequences
    :param s1: seq 1
    :param s2: seq 2
    :return: list of constraints
    """
    count = 0
    n, m = len(s1), len(s2)
    for i in range(n):  # upper left
        for j in range(m):  # lower right
            for l in range(i, n):  # upper right
                for k in range(j):  # lower left
                    if i != l or j != k:
                        constraints.append("x{}{}_{}{} + x{}{}_{}{} < 1;"
                                      .format(lidx, i, l2idx, j, lidx, l, l2idx, k))
                        count += 1

                        """
                        elif i == j and l != k:
                        cycles.append("X{}{}_{}{} + X{}{}_{}{} < 1;"
                                      .format(lidx, i, l2idx, j, lidx, l, l2idx, k))
                        count += 1

                    elif i != j and l == k:
                        cycles.append("X{}{}_{}{} + X{}{}_{}{} < 1;"
                                      .format(lidx, i, l2idx, j, lidx, l, l2idx, k))
                        count += 1
                        """
    return (constraints, count)


for a in range(n):  # upper left
    for b in range(m):  # lower right
        for c in range(a, n):  # upper right
            for d in range(m):  # lower left
                if (a != c or b != d) and d <= b:
                    cycles.append("X{}{}_{}{} + X{}{}_{}{} < 1;"
                                  .format(lidx, a, l2idx, b, lidx, c, l2idx, d))
                    count += 1

for a in range(n):  # upper left
    for b in range(m):  # lower right
        for c in range(a, n):  # upper right
            for d in range(b):  # lower left
                if a != c or b != d:
                    cycles.append("X{}{}_{}{} + X{}{}_{}{} < 1;"
                                  .format(lidx, a, l2idx, b, lidx, c, l2idx, d))
                    count += 1