import os
from Bio import Entrez, SeqIO
# how to get EMBOSS: http://emboss.open-bio.org/html/adm/ch01s01.html
with open("rosalind_need.txt", "r") as f:
    ids = f.readline().strip("\n").split()

Entrez.email = "your_name@your_mail_server.com"
handle = Entrez.efetch(db="nucleotide", id=[", ".join(ids)], rettype="fasta")
for record in SeqIO.parse(handle, "fasta"):
    with open(f"{record.name}.fasta", "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta")

emboss_pth = "tools/EMBOSS-6.6.0/emboss"
gapopen = 10.0
gapextend = 1.0
out_name = "needle_out.fasta"
cmd = f"{emboss_pth}/needle -asequence {ids[0]}.fasta -bsequence {ids[1]}.fasta\
    -datafile EDNAFULL -gapopen {gapopen} -gapextend {gapextend} -endopen {gapopen} -endextend {gapextend}\
        -outfile {out_name} -endweight Y"
os.system(cmd)
for id in ids:
    os.remove(f"{id}.fasta")
with open(out_name, "r") as f:
    for line in f:
        if line.startswith("# Score:"):
            res = line.strip().split(" ")[-1].split(".")[0]
            break
os.remove(out_name)
print(res)