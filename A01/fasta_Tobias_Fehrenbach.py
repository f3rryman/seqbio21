# Assignment 01 Tobias Fehrenbach

# Task 02
# reads a fasta file
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

        return data

# writes a fasta file
def write(fasta_pairs, file_name = None):
    if (file_name != None):
        with open(file_name, "w") as f:
            f.write("".join("%s \n%s\n" % x for x in fasta_pairs))
    else:
        print("Could not write to file: No file path specified!\nSequences to write to file:\n{}".format(fasta_pairs))