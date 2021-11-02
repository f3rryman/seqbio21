import sys
import fasta_Tobias_Fehrenbach as f

""" Task 04 """

# Codon Dictionary for DNA to AA translation
DNA_Codons = {
    # "M" - START, "_" - STOP
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TGT":"C", "TGC":"C",
    "GAT":"D", "GAC":"D",
    "GAA":"E", "GAG":"E",
    "TTT":"F", "TTC":"F",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "CAT":"H", "CAC":"H",
    "ATA":"I", "ATT":"I", "ATC":"I",
    "AAA":"K", "AAG":"K",
    "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "ATG":"M",
    "AAT":"N", "AAC":"N",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TGG":"W",
    "TAT":"Y", "TAC":"Y",
    "TAA":"_", "TAG":"_", "TGA":"_"
}


def translate(seq, init_pos=0):
    """Translate a DNA sequence into AA sequence"""
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
    with open(outPath, "w") as f:
        f.write("".join("%s \n%s\n" % x for x in out))
else:
    print("Header : Translated Sequence")
    print("".join("%s : %s\n" % x for x in out))
