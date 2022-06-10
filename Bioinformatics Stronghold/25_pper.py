from math import factorial

with open("rosalind_pper.txt", "r") as f:
    N, k = map(int, f.readline().strip().split())
print(int(factorial(N) / factorial(N-k)) % 1_000_000)