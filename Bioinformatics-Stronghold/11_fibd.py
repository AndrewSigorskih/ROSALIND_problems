with open("rosalind_fibd.txt", "r") as f:
    n, m = map(int, f.readline().strip().split())

f = [0] * (n + 1)
f[0] = 1
for i in range (2, n + 1):
    f[i] = f[i-2] + f[i-1] - f[i - m - 1]

print(f[n] + f[n-1])