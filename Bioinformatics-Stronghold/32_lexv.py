from itertools import product, chain
def main():
    with open("rosalind_lexv.txt", "r") as f:
        alph = f.readline().strip().replace(" ", "")
        n = int(f.readline().strip())
    strings = ["".join(x) for x in \
        list(chain(*[product(alph, repeat=k) for k in range(1, n+1)]))]
    with open("out.txt", "w") as o:
        print(*sorted(strings, key=lambda word: [alph.index(c) for c in word]), \
            sep="\n", file=o)


if __name__ == "__main__":
    main()