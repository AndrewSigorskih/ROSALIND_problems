from functools import lru_cache
from Bio import SeqIO

import sys
sys.setrecursionlimit(10_000)

def levenstein_distance(a, b):
    
    @lru_cache(None)  # for memorization
    def min_dist(s1, s2):
        if s1 == len(a) or s2 == len(b):
            return len(a) - s1 + len(b) - s2

        # no change required
        if a[s1] == b[s2]:
            return min_dist(s1 + 1, s2 + 1)

        return 1 + min(
            min_dist(s1, s2 + 1),      # insert character
            min_dist(s1 + 1, s2),      # delete character
            min_dist(s1 + 1, s2 + 1),  # replace character
        )

    return min_dist(0, 0)

def main():
    seq1, seq2 = (x.seq for x in SeqIO.parse("rosalind_edit.txt", "fasta"))
    with open("out.txt", "w") as o:
        print(levenstein_distance(seq1, seq2), file=o)

if __name__ == "__main__":
    main()