import fasta_Garcia_Fehrenbach as f

""" Task 03 """
# read fasta file & print in to console
# ./assign-01-dna.fasta
path = input("Type file path of existing fasta & press ENTER\n")
fasta = f.read(path)
f.write(fasta)