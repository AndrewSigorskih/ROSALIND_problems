from math import factorial
import itertools

signs = (-1, 1)
with open("rosalind_sign.txt", "r") as f:
    n = int(f.readline().strip())

f = open("out.txt", "w")
print(factorial(n) * 2**n, file=f)
for item in itertools.permutations(range(1, n+1)):
    for perm in itertools.product(signs, repeat=n):
        print(*(x*y for x,y in zip(perm, item)), file=f)
f.close()