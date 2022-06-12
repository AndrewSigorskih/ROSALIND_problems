import re
import requests
from io import StringIO
from Bio import SeqIO

pat = re.compile(r"(?=(N[^P][ST][^P]))")

def return_pos(id):
    handle = requests.get(f"http://www.uniprot.org/uniprot/{id}.fasta").text
    handle = StringIO(handle)
    record = SeqIO.read(handle, "fasta")
    return [m.start()+1 for  m in re.finditer(pat, str(record.seq))]


with open("rosalind_mprt.txt", 'r') as f, open("out.txt", "w") as o:
    ids = (line.strip() for line in f)
    for id in ids:
        res = return_pos(id)
        if res:
            print(id, file=o)
            print(*res, file=o)