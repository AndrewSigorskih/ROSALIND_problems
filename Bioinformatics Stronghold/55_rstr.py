def main():
    with open("rosalind_rstr.txt", "r") as f:
        N, x = f.readline().strip('\n').split()
        s = f.readline().strip('\n')
    N = int(N)
    x = float(x)
    at = s.count('A') + s.count('T')
    gc = s.count('G') + s.count('C')    
    prob = pow(x/2, gc) * pow((1-x)/2, at)

    with open("out.txt", "w") as o:
        print(f"{(1 - pow(1-prob, N)):.3f}", file=o)

if __name__ == "__main__":
    main()