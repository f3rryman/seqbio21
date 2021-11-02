import fasta_Tobias_Fehrenbach as f

""" Task 03 """
# read fasta file & print in to console
path = input("Type file path of existing fasta & press ENTER\n")
fasta = f.read(path)
f.write(fasta, "/home/tobif/Documents/MasterBioinfo/Semester03/SequenceBioinformatics/assignments/assignment01/fasta_test.fasta")