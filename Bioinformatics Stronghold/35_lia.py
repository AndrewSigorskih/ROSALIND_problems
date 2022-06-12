from math import factorial

def C_n_k(n, k):
    return factorial(n) // (factorial(k) * factorial(n-k))

def calc_prob(n, k, p=.25):
    return (C_n_k(n, k)) * (p ** k) * (1 - p) ** (n - k)

def main():
    with open("rosalind_lia.txt", "r") as f:
        k, N = map(int, f.readline().strip().split())
    
    prob = sum(calc_prob(2**k, i) for i in range(N, 2**k + 1))
    print(f"{prob:.3f}")
        
if __name__ == "__main__":
    main()