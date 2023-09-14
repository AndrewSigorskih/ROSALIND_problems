with open("rosalind_sset.txt", "r") as f:
    N = int(f.readline().strip())
print(2**N % 1_000_000)