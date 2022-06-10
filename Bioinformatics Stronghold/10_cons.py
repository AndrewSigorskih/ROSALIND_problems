import pandas as pd
from Bio import SeqIO

order = ("A", "C", "G", "T")

records = SeqIO.parse("rosalind_cons.txt", "fasta")
seqs = pd.DataFrame([list(record.seq) for record in records])

profile = pd.concat((seqs[i].value_counts() for i in seqs.columns), axis=1)
profile = profile.fillna(0).reindex(order).astype(int)

cons = "".join([profile[i].idxmax() for i in profile.columns])

with open("out.txt", "w") as f:
    print(cons, end="\n", file=f)
    for c in order:
        print(f"{c}:", *profile.loc[c].values, \
            sep=" ", end="\n", file=f)