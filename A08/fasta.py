""" FastA file handler
    :functions: read(path)"""

def read(path):
    """reads in fasta file"""
    heads = []
    head = ""
    seqs = []
    seq = ""
    with open(path, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if seq != "": # append seq + head to their lists
                    seqs.append(seq)
                    heads.append(head)
                    seq = ""
                head = line # save current line to head if it startswith >
            else:
                seq += line # add current line to seq -> necessary to read in each sequence line
        seqs.append(seq) # add last seq to seqs array
        heads.append(head) # add last header to heads array
        data = [(heads[i], seqs[i]) for i in range(0, len(heads))] # combine to data

    f.close()
    return heads, seqs
