from math import factorial
def C_n_k(n, k):
    return (factorial(n) // (factorial(k) * factorial(n-k)))

with open("rosalind_aspc.txt", "r") as f:
    N, m = map(int, f.readline().strip().split())
print(sum(C_n_k(N, k) for k in range(m, N+1)) % 1_000_000)
