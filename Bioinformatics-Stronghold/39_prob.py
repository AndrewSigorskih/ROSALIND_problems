from math import log10

def count_prob(S : str, GC : float):
    return sum(log10(GC/2) if x in ("G", "C") \
        else log10((1 - GC) / 2) for x in S)

def main():
    with open("rosalind_prob.txt", "r") as f:
        s = f.readline().strip()
        A = list(map(float, f.readline().strip().split()))
    with open("out.txt", "w") as o:
        print(*(f"{count_prob(s, gc):.3f}" for gc in A), file=o)

if __name__ == "__main__":
    main()