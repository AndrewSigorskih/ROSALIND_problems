from Bio import Seq
import warnings
from itertools import product
warnings.filterwarnings("ignore")

with open("rosalind_orfr.txt", "r") as f:
    seq = Seq.Seq(f.readline().strip())
    rseq = seq.reverse_complement()
orfs = []
for i, s in product(range(3), [seq, rseq]):
    products = s[i:].translate().__str__().split('*')
    products = [p[p.find("M"):] for p in products if "M" in p]
    if products:
        orfs.append(max(products, key=len))
print(max(orfs, key=len))