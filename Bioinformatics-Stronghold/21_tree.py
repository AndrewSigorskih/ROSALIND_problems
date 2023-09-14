with open("rosalind_tree.txt", "r") as f:
    N = int(f.readline().strip())
    nedges = sum((1 for line in f))
print(N-1-nedges)