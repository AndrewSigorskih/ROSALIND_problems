from collections import Counter
from itertools import product

def main():
    with open("rosalind_conv.txt", "r") as f:
        s1, s2 = map(lambda x: \
            [float(j) for j in x.strip().split()],\
                (l for l in f))
    cnt = Counter((round(i - j, 5) for i, j in product(s1, s2)))
    res, num = cnt.most_common(1)[0]
    with open("out.txt", "w") as o:
        print(num, res, sep="\n", file=o)

if __name__ == "__main__":
    main()