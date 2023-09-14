from Bio import SeqIO, Seq
from itertools import product
import re

START = "ATG"
STOP = re.compile('(TAA|TAG|TGA)')

def get_orfs(s):
    uniq_prots = set()
    for i in range(0, len(s), 3):
        if (len(s) - i) < 3:
            break
        if (s[i:i+3] == START):
            for j in range(i, len(s), 3):
                if bool(re.match(STOP, s[j:j+3])):
                    uniq_prots.add(Seq.translate(s[i:j+3], \
                        to_stop=True, stop_symbol=""))
                    break
    return uniq_prots

if __name__ == "__main__":
    frw = SeqIO.read("rosalind_orf.txt", "fasta").seq
    rev = frw.reverse_complement()
    res = set()
    for i, s in product(range(3), [frw, rev]):
        res |= get_orfs(str(s)[i:])
    with open("out.txt", "w") as f:
        print(*list(res), sep="\n", file=f)
