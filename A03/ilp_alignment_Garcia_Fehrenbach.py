import argparse

""" Task 02 Simple mixed cycles """
s1 = "CGTA"
s2 = "CGAT"
s3 = "GATAT"
seqs = [s1, s2, s3]
len_seqs = len(seqs)


def calc_2constraints(seqs):
    """
    for each seq in seqs this function sets up crossovers between each two seq
    it appends the inequality to a list
    it calculates the number of inequalities
    example:
    00 01
    10 11
    -> x00_11 + x01_10 < 1
    lidx  a  lidx  c
    l2idx d  l2idx b
    -> xlidx a_l2idx b + xlidx c_l2idx d < 1
    :param seqs: list of sequences to compare
    :return: (inequalities, number of inequalities
    """
    count = 0
    cycles = []
    # number of seqs
    nr_seqs = len(seqs)
    for lidx in range(nr_seqs):
        for l2idx in range(1, nr_seqs):
            if lidx < l2idx: # compare seq 0 - 1, 0 - 2, 1 - 2
                #print(lidx, l2idx)
                n = len(seqs[lidx])
                m = len(seqs[l2idx])

                # cross overs between two seq
                for a in range(n):  # upper left
                    for b in range(m):  # lower right
                        for c in range(a, n):  # upper right
                            for d in range(m):  # lower left
                                if (a != c or b != d) and d <= b:
                                    cycles.append("X{}{}_{}{} + X{}{}_{}{} < 1;"
                                                  .format(lidx, a, l2idx, b, lidx, c, l2idx, d))
                                    count += 1

    cycles = sorted(cycles)
    return (cycles, count)


def calc_3constraints(seqs):
    """
    calculates all constraints between three sequences given a list of three sequences and returns them as inequalities
    :param seqs: list of sequences to compare
    :return: (inequalities, number of inequalities)
    """
    count = 0
    cycles = []
    # number of seqs
    nr_seqs = len(seqs)

    n = len(seqs[0])
    m = len(seqs[1])
    o = len(seqs[2])

    for a in range(n):
        for b in range(m):
            for c in range(o):
                for d in range(c-1):
                    for e in range(a+1, n):
                        cycles.append("X{}{}_{}{} + X{}{}_{}{} + X{}{}_{}{} < 2;"
                                      .format(0, a, 1, b, 1, b, 2, c, 2, d, 0, e))
                        count += 1

                for d in range(c+1, o):
                    for e in range(a-1):
                        cycles.append("X{}{}_{}{} + X{}{}_{}{} + X{}{}_{}{} < 2;"
                                      .format(0, a, 1, b, 1, b, 2, c, 2, d, 0, e))
                        count += 1

    cycles = sorted(cycles)
    return (cycles, count)


def calc_binary_constraints(seqs):
    """
    Calculates each edge between each char in each sequence and returns it as constraint
    e.g.
    X00_10<1;
    X00_11<1;
    X00_12<1;
    :param seqs: list of sequences
    :return: (bin_cons, count_of_bincons)
    """
    bin_cons = []
    count = 0
    nr_seqs = len(seqs)
    for lidx in range(nr_seqs):
        for l2idx in range(1, nr_seqs):
            if lidx < l2idx:
                n, m = len(seqs[lidx]), len(seqs[l2idx])
                for i in range(n):
                    for j in range(m):
                        #print("X{}{}_{}{}<1;".format(lidx, i, l2idx, j))
                        bin_cons.append("X{}{}_{}{}<1;".format(lidx, i, l2idx, j))
                        count += 1
    bin_cons = sorted(bin_cons)
    return (bin_cons, count)


def is_match(chr1, chr2):
    """
    compares two chars returns true if match
    :param chr1: char
    :param chr2: char
    :return: if chr1 == chr2 return true
    """
    if chr1 == chr2:
        return True
    else:
        return False


def setup_objective_function(seqs):
    """
    sets up the objective function in the form:
    max 1*X11_21+4*X12_21+ ... ;
    with match = 4, and mismatch = 1
    :param seqs: list of sequences
    :return: objective function
    """
    obj_entr = []
    obj_fct = "max: "
    nr_seqs = len(seqs)
    for lidx in range(nr_seqs):
        for l2idx in range(1, nr_seqs):
            if lidx < l2idx:
                n, m = len(seqs[lidx]), len(seqs[l2idx])
                for i in range(n):
                    for j in range(m):
                        if is_match(seqs[lidx][i], seqs[l2idx][j]):
                            obj_entr.append("4*X{}{}_{}{}".format(lidx, i, l2idx, j))
                        else:
                            obj_entr.append("1*X{}{}_{}{}".format(lidx, i, l2idx, j))
    obj_entr = sorted(obj_entr)
    for i in obj_entr:
        obj_fct += i  + "+"
    obj_fct = obj_fct[:-1] + ";"
    return obj_fct


def all_edges(seqs):
    """
    Calculates each edge between each char in each sequence
    :param seqs: list of sequences
    :return: (bin_cons, count_of_bincons)
    """
    es = []
    count = 0
    nr_seqs = len(seqs)
    for lidx in range(nr_seqs):
        for l2idx in range(1, nr_seqs):
            if lidx < l2idx:
                n, m = len(seqs[lidx]), len(seqs[l2idx])
                for i in range(n):
                    for j in range(m):
                        #print("X{}{}_{}{}<1;".format(lidx, i, l2idx, j))
                        es.append("X{}{}_{}{}".format(lidx, i, l2idx, j))
                        count += 1
    es = sorted(es)
    return (es, count)


def edges_as_int(es):
    es_as_int = "".join("int ")
    for i in es:
        es_as_int += i + ", "
    es_as_int = es_as_int[:-2] + ";"
    return es_as_int


# writes a fasta file
def write(obj_fct, cycles, bin_cons, es_as_int, file_path = None):
    """
    writes it parameters to an .lp file for solving the inequalities with lp_solve "http://lpsolve.sourceforge.net/"
    :param obj_fct: objective function: e.g. max 1*X11_21+4*X12_21+ ... ;
    :param cycles: simple mixed cycles constraints: e.g. x00_11 + x01_10 < 1;
                                                         ...
    :param bin_cons: binary constraints: e.g.   X00_10<1;
                                                X00_11<1;
                                                X00_12<1;
    :param es_as_int: all edges have to be defined as integers in the last row: e.g. int X00_10, X00_11, X00_12, ...;
    :param file_path: file path, has to end with .lp
    :return:
    """
    if (file_path != None) and file_path[-3:] == ".lp":
        with open(file_path, "w") as f:
            f.write("".join(obj_fct+ "\n\n"))
            f.write("".join("%s\n" % c for c in cycles))
            f.write("\n\n")
            f.write("".join("%s\n" % bc for bc in bin_cons))
            f.write("\n\n")
            f.write("".join(es_as_int))
            print(file_path + "\n successfully written.")
    else:
        print("EITHER NO FILEPATH SPECIFIED OR WRONG ENDING. FILE HAS TO END WITH '.lp'")



""" MAIN """
# creating parser for out file specification -> in.lp
parser = argparse.ArgumentParser(description="Creates an .lp file for three sequences."
                                             "With lp_solve /file/path/to/in.lp"
                                             "from directory of lp_solve")
parser.add_argument("Path",
                    metavar="path",
                    type=str,
                    help="Output path of lp file.")


args = parser.parse_args()
inPath = args.Path

# set up the objective function
obj_fct = setup_objective_function(seqs)

# set up inequalities for all simple mixed cycles
cycles2, no_of_cycles2 = calc_2constraints(seqs)
cycles3, no_of_cycles3 = calc_3constraints(seqs)

cycles = cycles2 + cycles3
no_of_cycles = no_of_cycles2 + no_of_cycles3

# set up binary constraints
bin_cons, bin_cons_count = calc_binary_constraints(seqs)

# set up all edges the graph contains
es, no_of_es = all_edges(seqs)
es_as_int = edges_as_int(es)

# write lp to file specified as first argument
write(obj_fct, cycles, bin_cons, es_as_int, inPath)
print("Number of constraints:\n2s constraints: {}\n3s constraints: {}\ncombined: {}"
      .format(no_of_cycles2, no_of_cycles3, no_of_cycles))
