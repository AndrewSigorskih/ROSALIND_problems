# surprisingly easy one
def main():
    with open("rosalind_ebin.txt", "r") as f:
        n = int(f.readline())
        P = list(map(float, f.readline().split()))
    with open("out.txt", "w") as o:
        print(*(round(n*p, 3) for p in P), file=o)

if __name__ == "__main__":
    main()