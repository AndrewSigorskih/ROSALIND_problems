from Bio.Seq import Seq
from Bio import SeqIO

instream = SeqIO.parse("rosalind_splc.txt","fasta")
target = str(next(instream).seq)

for intron in instream:
    target = target.replace(str(intron.seq),"")

print(Seq(target).translate(stop_symbol=""))


