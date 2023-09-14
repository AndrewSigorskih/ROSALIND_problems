with open("rosalind_iev.txt", "r") as f:
    v, w, x, y, z, _ = map(int, f.readline().strip().split())
print(2*(v+w+x) + 1.5*y + z)