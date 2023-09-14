with open("rosalind_seto.txt", "r") as f:
    N = int(f.readline().strip())
    A, B = map(eval, (x.strip() for x in f))
U = set(range(1, N+1))
with open("out.txt", "w") as f:
    print(A.union(B), file=f)
    print(A.intersection(B), file=f)
    print(A - B, file=f)
    print(B - A, file=f)
    print(U - A, file=f)
    print(U - B, file=f)