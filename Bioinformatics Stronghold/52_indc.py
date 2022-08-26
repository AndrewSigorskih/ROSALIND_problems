from math import factorial, log10
P = .5

def C_n_k(n: int, k: int) -> int:
    return factorial(n) // (factorial(k) * factorial(n-k))

def indc(n: int) -> list:
    prob = 0
    res = []
    for i in range(2*n, 0, -1):
        prob += C_n_k(n*2, i) * ((P**i) * (1-P)**(2*n - i))
        res.append(f"{log10(prob):.3f}")
    return res[::-1]

def main():
    with open("rosalind_indc.txt", "r") as f:
        n = int(f.readline().strip('\n'))
    
    with open("out.txt", "w") as o:
        print(*indc(n), sep=" ", file=o)

if __name__ == "__main__":
    main()