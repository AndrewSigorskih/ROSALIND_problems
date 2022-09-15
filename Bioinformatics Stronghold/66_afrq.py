from math import sqrt

def main():
    with open("rosalind_afrq.txt", "r") as f:
        A = list(map(float, f.readline().split()))
    with open("out.txt", "w") as o:
        print(*[f"{(2*sqrt(i)-i):.3f}" for i in A], sep=' ', file=o)

if __name__ == "__main__":
    main()