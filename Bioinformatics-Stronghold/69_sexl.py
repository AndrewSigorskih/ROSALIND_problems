def main():
    with open("rosalind_sexl.txt", "r") as f:
        A = list(map(float, f.readline().split()))
    with open("out.txt", "w") as o:
        print(*[f"{(2*(i - i**2)):.3f}" for i in A], sep=' ', file=o)

if __name__ == "__main__":
    main()