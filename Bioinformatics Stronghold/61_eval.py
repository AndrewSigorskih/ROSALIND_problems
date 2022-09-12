def main():
    with open("rosalind_eval.txt", "r") as f:
        n = int(f.readline())
        s = f.readline().strip()
        A = map(float, f.readline().strip().split())
    at = s.count('A') + s.count('T')
    gc = s.count('G') + s.count('C')
    res = (pow((1-x)/2, at) * pow(x/2, gc) * (n-len(s)+1) for x in A)
    with open("out.txt", "w") as o:
        print(*res, sep=" ", file=o)

if __name__ == "__main__":
    main()