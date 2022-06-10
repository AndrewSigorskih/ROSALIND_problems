from math import factorial
from itertools import permutations

with open("rosalind_perm.txt", "r") as f:
    n = int(f.readline().strip())
with open("out.txt", "w") as f:
    print(factorial(n), file=f)
    for perm in permutations(range(1, n+1)):
        print(*perm, file=f)