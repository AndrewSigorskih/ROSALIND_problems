codonNums = {
    "I" : 3,
    "L" : 6,
    "V" : 4,
    "F" : 2,
    "M" : 1,
    "C" : 2,
    "A" : 4,
    "G" : 4,
    "P" : 4,
    "T" : 4,
    "S" : 6,
    "Y" : 2,
    "W" : 1,
    "Q" : 2,
    "N" : 2,
    "H" : 2,
    "E" : 2,
    "D" : 2,
    "K" : 2,
    "R" : 6,
    "Stop": 3
}

with open("rosalind_mrna.txt", "r") as f:
    prot = f.readline().strip()

res = codonNums["Stop"]
for c in prot:
    res *= codonNums[c]
print(res % 1_000_000)