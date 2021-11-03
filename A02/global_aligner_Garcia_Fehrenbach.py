import os
import sys
import fasta
import argparse
import nw
import tracemalloc
import time


"""Handling command line arguments"""
parser = argparse.ArgumentParser(description="Aligns two sequences from fasta file with Needleman-Wunsch variations.")
parser.add_argument("Path",
                    metavar="path",
                    type=str,
                    help="Input path of fasta file.")
parser.add_argument("--mode",
                    action="store",
                    type=int,
                    help="0: NW basic (default), 1: NW with linear space, 2: NW without table",
                    default=0)
parser.add_argument("--stat",
                    action="store_true",
                    help="If set prints out time and memory consumption the specified algorithm needed for calculation.")

args = parser.parse_args()
inPath = args.Path
mode = args.mode
stat = args.stat

# break if no file path is specified
if not os.path.isfile(inPath):
    print("Specified path does not exist.")
    sys.exit()

"""get data & parse needleman wunsch, this will only work for the first 2 seqs specified in fasta"""
data = fasta.read(inPath)
x = data[0][1]
y = data[1][1]

"""run calculation with statistics"""
if stat:
    if mode == 0:
        tracemalloc.start()
        start_time = time.time()
        sltn = nw.nw_basic(x, y)
        current, peak = tracemalloc.get_traced_memory()
        print("The optimal global alignment of {} vs. {} is:\n{}\n{}"
              .format(data[0][0], data[1][0], sltn[0], sltn[1]))
        print("With a score of {}.".format(sltn[2]))
        print("The algorithm needed {} s to align both sequences.".format(time.time()-start_time))
        print("Memory usage: current = {} bytes, peak = {} bytes".format(current, peak))
        tracemalloc.reset_peak()

    elif mode == 1:
        tracemalloc.start()
        start_time = time.time()
        sltn = nw.nw_linear_space(x,y)#nw.nw_hirschberg(x, y)
        current, peak = tracemalloc.get_traced_memory()
        print("The optimal global alignment using the Hirschberg algorithm is:\n{}\n{}".format(sltn[0], sltn[1]))
        print("With a score of {}.".format(sltn[2]))
        print("The algorithm needed {} s to align".format(time.time() - start_time))
        print("Memory usage: current = {} bytes, peak = {} bytes".format(current, peak))
        tracemalloc.reset_peak()

    elif mode == 2:
        tracemalloc.start()
        start_time = time.time()
        sltn = nw.nw_wo_table(x, y)
        current, peak = tracemalloc.get_traced_memory()
        #print("A possible optimal global alignment using NW without table is:\n{}\n{}".format(sltn[0], sltn[1]))
        print("With a score of {}.".format(sltn))
        print("The algorithm needed {} s to align".format(time.time() - start_time))
        print("Memory usage: current = {} bytes, peak = {} bytes".format(current, peak))
        tracemalloc.reset_peak()

    else:
        print("No viable mode. Mode must be either 0, 1 or 2")

else:
    """run calculation without statistics"""
    if mode == 0:
        sltn = nw.nw_basic(x, y)
        print("The optimal global alignment of {} vs. {} is:\n{}\n{}"
              .format(data[0][0], data[1][0], sltn[0], sltn[1]))
        print("With a score of {}.".format(sltn[2]))

    elif mode == 1:
        sltn = nw.nw_hirschberg(x, y)
        print("The optimal global alignment using the Hirschberg algorithm is:\n{}\n{}".format(sltn[0], sltn[1]))
        print("With a score of {}.".format(sltn[2]))

    elif mode == 2:
        sltn = nw.nw_wo_table(x, y)
        #print("A possible optimal global alignment using NW without table is:\n{}\n{}".format(sltn[0], sltn[1]))
        print("With a score of {}.".format(sltn))

    else:
        print("No viable mode. Mode must be either 0, 1 or 2")