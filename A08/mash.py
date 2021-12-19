import argparse
import math
from datetime import datetime
import os.path
import mmh3 # hash function
import bisect
import fasta


class Mash:
    def __init__(self, k:int, s:int):
        self.k = k
        self.s = s
        self.input = []

    def create_kmers(self, k: int, input: []):
        """
        Creates kmers of size k from list of contigs of a genome.
        :param k: size of kmer
        :param input: list of sequences
        :return: list of all kmers of the genome
        """
        kmer_list = []  # list for all kmers of file
        for seq in input:  # for sequence in input file
            kmers = []
            for i in range(len(seq) - k + 1):  # create kmers
                kmer = seq[i:i + k]  # each kmer ranges from i to i+k
                kmers.append(kmer)  # kmers for each contig
            kmer_list.extend(kmers)  # kmers for whole genome
        return kmer_list

    def calc_sketch(self, k: int, s: int, input: []):
        """
        Calculates sketches of given input sequence input with k-mer size k and sketch size s.
        :param k: k-mer size
        :param s: sketch size
        :param input: input sequence
        :return: sorted list q of sketches -> descending order
        """
        q = sorted([]) # empty sorted list
        # create kmers
        kmer_list = self.create_kmers(k, input) # list for all kmers of file

        # create sketches
        for kmer in kmer_list: # go through kmers in list
            h = mmh3.hash64(kmer, 100) # hash current kmer, default seed is 0, seed 100 seems fine????
            if not q: # if q is empty -> append the hash
                q.append(h)
            elif h < q[-1] or len(q) < s: # if q not empty and if hash < max hash or q have not reached s
                if h not in q: # check if hash already exists in q
                    bisect.insort(q, h) # if not: insort the hash
                    if len(q) > s: # if q exceeds sketch size -> remove the highest value
                        q.remove(q[-1])
        #print(q[-1]) # is max(q) and q[-1] the same? yes, q[-1] increases the runtime
        #print(max(q))
        return q

    def calc_jaccard(self, sa: [], sb: [], s: int):
        """
        Calculates the jaccard index for two sketches.
        :param sa: sketch 1
        :param sb: sketch 2
        :return: jaccard index(sa, sb)
        """
        # creating sets of sa, sb destroys the order, which is not relevant for us since we work only with the lengths
        # further on
        s1 = set(sa)
        s2 = set(sb)

        intersect = len(s1.intersection(s2))
        uni = len(s1.union(s2))
        if uni > s: # cut union to s hashes????
            uni = len(list(s1.union(s2))[0:s])

        return intersect/uni

    def calc_distance(self, k: int, s: int, sa: [], sb: []):
        """
        Calculates the mash distance for two sketches with k-mer size k.
        :param k: k-mer size
        :param sa: sketch 1
        :param sb: sketch 2
        :return: mash distance(k, sa, sb)s
        """
        jis = self.calc_jaccard(sa, sb, s)
        if jis != 0:
            d_mash = -1/k * math.log(2*jis/(1+jis))
        else:
            d_mash = 1 # log of 0 is not defined -> distance is "maximal" -> no similarity between sequences -> 1????
        return d_mash

    def xchange_base(self, seqs: []):
        """
        Function for debugging. Exchanges the bases in a given genome = list of contigs. If we input the same sequence
        twice in mash and then use this function on one of them. The dissimilarity should get higher -> from 0 to xx
        :param seqs: list of contigs
        :return: list of contigs with exchanged bases
        """
        new_seqs = []
        for seq in seqs:
            new_seq = seq.replace("A", "T")
            #new_seq = new_seq.replace("G", "C")
            new_seqs.append(new_seq)
        return new_seqs


def main():
    """
    Main to run.
    :return:
    """
    parser = argparse.ArgumentParser(description="Simple implementation of the Mash algorithm. Computes sketches, "
                                                 "the jaccard index and Mash distances from a set of FastA files. "
                                                 "K-mer size and sketch size can be specified.")
    parser.add_argument("command",
                        type=str,
                        help="Specifies what should be calculated. Either sketch, jaccard or distances.")
    parser.add_argument("-k",
                        type=int,
                        help="Specify k-mer size. Default = 21.",
                        default=21)
    parser.add_argument("-s",
                        type=int,
                        help="Specify sketch size. Default = 1000.",
                        default=1000)
    parser.add_argument("input_files",
                        type=str,
                        nargs="+",
                        help="Input FastA files with sequences from which sketches, jaccard indices and distances are "
                             "calculated.")
    args = parser.parse_args()
    command = args.command
    k = args.k
    s = args.s
    input = args.input_files

    # handle fasta files
    filenames = []
    head_list = []
    seq_list = []
    for idx in range(len(input)):
        filenames.append(os.path.basename(input[idx]))
        head_list.append(fasta.read(input[idx])[0])
        seq_list.append(fasta.read(input[idx])[1])

    # check number of k-mers L - k+1
    # for debugging
    # for s1 in seq_list:
    #     L = 0
    #     for s2 in s1:
    #         # handling Ns within the sequences???
    #         #if s2.__contains__("N"):
    #         #    s2.replace("N", "")
    #         L += len(s2)-k+1
    #     print("Number of k-mers in genome: {}".format(L))

    # create mash
    mash = Mash(k, s)

    # more debugging
    # exchange a base to compare seq with similar seq with exchanged base
    #seq_list[1] = mash.xchange_base(seq_list[1])
    #print(seq_list[0][0][:10])
    #print(seq_list[1][0][:10])

    # calculate sketches of each file
    print("Calculating sketches. This may take a while ...")
    print(datetime.now())
    sketches = []
    for seq in seq_list: # seq is a list of contigs per file -> whole genome in contigs
        sketch = mash.calc_sketch(k=k, s=s, input=seq)  # calc sketch per file, sketch is list of s hashes
        sketches.append(sketch) # append sketch per file to list sketches
    print(datetime.now())

    # Output all sketches: name of file -> hash value per line -> descending order -> highest to lowest
    if command == "sketch":
        for idx in range(len(filenames)):
            print("Sketch of file: {}".format(filenames[idx]))
            for sketch in sketches:
                for this_hash in sketch[::-1]:
                    print(this_hash)

    if command == "jaccard" or "distances":
        if command=="jaccard":
            # for each pair of input files
            # compute jaccard index
            # store jaccard indices to matrix
            print("Matrix of Jaccard indices:")
            jaccard_matrix = []
            for row in range(len(sketches)):
                jaccard_matrix.append([])
                for col in range(len(sketches)):
                    jaccard_matrix[row].append(mash.calc_jaccard(sketches[row], sketches[col], s))
            # print jaccard index to console
            print("\t", len(jaccard_matrix), end="")
            for row in range(len(jaccard_matrix)):
                print()
                print(filenames[row], " ", end="")
                for entry in jaccard_matrix[row]:
                    print(entry, " ", end="")

        else:
            # for each pair of input files
            # compute mash distances
            # store mash distances to matrix
            print("Matrix of mash distances: ")
            dist_matrix = []
            for row in range(len(sketches)):
                dist_matrix.append([])
                for col in range(len(sketches)):
                    dist_matrix[row].append(mash.calc_distance(k, s, sketches[row], sketches[col]))
            # print mash distances to console
            print("\t", len(dist_matrix), end="")
            for row in range(len(dist_matrix)):
                print()
                print(filenames[row], " ", end="")
                for entry in dist_matrix[row]:
                    print(entry, " ", end="")


if __name__ == "__main__":
    main()
