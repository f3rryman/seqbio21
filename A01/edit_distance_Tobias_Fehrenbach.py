import sys
import fasta_Tobias_Fehrenbach as f
import numpy as np
import pandas as pd

""" variables """
inPath = sys.argv[1]
header = []
seq = []
result = []

""" read file """
data = f.read(inPath)
for i in range(len(data)):
    header.append(data[i][0])
    seq.append(data[i][1])

# calculate the edit distance d(s, t) with d(s, t) = # of positions where s, t differ from each other
# when both sequences are not of the same length, those sequences will be skipped
matrix = np.zeros([len(header), len(header)], int)
for idx, i in enumerate(seq):
    for idx2, j in enumerate(seq):
        d = 0
        h1 = header[idx]
        h2 = header[idx2]
        if len(i) != len(j):    # Throw a warning -> program doesnt do anything -> dont fill with gaps
            print("Sequences of {} and {} are not of the same length. Skipping edit distance calculation.")
            continue
        else:
            d = sum(1 for (l, r) in zip(i, j) if l != r)

        #print("idx: {}, idx2: {}, d: {}".format(idx, idx2, d))
        result.append([(h1 + " - " + h2, d)])
        matrix[idx, idx2] = d

out = pd.DataFrame(matrix, columns=header, index=header)
# print matrix without truncation
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
#print(result)
#print(matrix)
print("d(s, t) of all sequences from given file:")
print(out)