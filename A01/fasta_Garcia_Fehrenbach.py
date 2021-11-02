""" FastA file handler
    :functions: read(path), write(data, file_path=None)"""

"""Task 02"""

def read(path):
    """reads in fasta file
        :returns: data as list of tuple (header, sequence)"""
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
def write(data, file_path = None):
    """writes data as tuple in form of (header, sequence)
        to a file if file_name is specified.
        else prints data to console"""
    if (file_path != None):
        with open(file_path, "w") as f:
            f.write("".join("%s \n%s\n" % x for x in data))
    else:
        print("Could not write to file: No file path specified!")
        print("Header : Sequence")
        print("".join("%s : %s\n" % x for x in data))