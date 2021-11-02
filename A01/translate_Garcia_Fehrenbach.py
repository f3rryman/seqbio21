import sys
import fasta_Garcia_Fehrenbach as f

""" Task 04 """

# Codon Dictionary for DNA to AA translation
DNA_Codons = {
    # "M" - START, "_" - STOP
    "CGA":"A", "CGG":"A", "CGT":"A", "CGC":"A",
    "ACA":"C", "ACG":"C",
    "CTA":"D", "CTG":"D",
    "CTT":"E", "CTC":"E",
    "AAA":"F", "AAG":"F",
    "CCA":"G", "CCG":"G", "CCT":"G", "CCC":"G",
    "GTA":"H", "GTG":"H",
    "TAT":"I", "TAA":"I", "TAG":"I",
    "TTT":"K", "TTC":"K",
    "AAT":"L", "AAC":"L", "GAA":"L", "GAG":"L", "GAT":"L", "GAC":"L",
    "TAC":"M",
    "TTA":"N", "TTG":"N",
    "GGA":"P", "GGG":"P", "GGT":"P", "GGC":"P",
    "GTT":"Q", "GTC":"Q",
    "GCA":"R", "GCG":"R", "GCT":"R", "GCC":"R", "TCT":"R", "TCC":"R",
    "AGA":"S", "AGG":"S", "AGT":"S", "AGC":"S", "TCA":"S", "TCG":"S",
    "TGA":"T", "TGG":"T", "TGT":"T", "TGC":"T",
    "CAA":"V", "CAG":"V", "CAT":"V", "CAC":"V",
    "ACC":"W",
    "ATA":"Y", "ATG":"Y",
    "ATT":"_", "ATC":"_", "ACT":"_"
}


def translate(seq, init_pos=0):
    """Translate a DNA sequence into AA sequence
        :returns: AA sequence as string"""
    sAA = [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) - 2, 3)] # calcs single AAs from seq
    seqAA = "".join(sAA) # join sAA list to AA seq
    return seqAA

# variables
inPath = sys.argv[1]
if len(sys.argv) == 3:
    outPath = sys.argv[2]
header = []
seq = []
tseq = []

# read file
data = f.read(inPath)
for i in range(len(data)):
    header.append(data[i][0])
    seq.append(data[i][1])

# translate
for s in seq:
    tseq.append(translate(s))

out = [(header[i], tseq[i]) for i in range(0, len(header))]

# write to console if outPath is not specified as 2nd argument
if len(sys.argv) == 3:
    f.write(out, outPath)
    print("Translated sequences written to {}".format(outPath))
else:
    f.write(out)
