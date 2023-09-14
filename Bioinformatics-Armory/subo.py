import os
from Bio import SeqIO
pth_to_lalign = "./tools/fasta-36.3.8h/bin/lalign36"

def hamming(s1, s2, stop=3):
    hammingnum = 0
    for x, y in zip(s1, s2):
        if x != y:
            hammingnum += 1
        if hammingnum > stop:
            break
    return hammingnum

def calc_inexact(a, b):
    c = 0
    for i in range(len(a) - len(b)):
        if hamming(a[i:i+len(s)], s)<=3:
            c += 1
    return c  

def run_lalign(name1, name2):
    print(name1, name2)
    res = 0
    cmd = f'{pth_to_lalign} -C 13 {name1}.fasta {name2}.fasta | grep "100" -A 3 | tail -n 1 > out.txt'
    #print(cmd)
    os.system(cmd)
    with open("out.txt", "r") as f:
        output = f.readlines()[-1].strip().split()[-1]
    return output


if __name__ == '__main__':
    records = list(SeqIO.parse("rosalind_subo.txt", "fasta"))
    for record in records:
        SeqIO.write(record, f"{record.name}.fasta", "fasta")
    
    s = run_lalign(records[0].name, records[1].name)
    print(s)
    print(calc_inexact(records[0], s), calc_inexact(records[1], s))
    