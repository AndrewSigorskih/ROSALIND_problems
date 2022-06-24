# python 3.8 or higher
from math import prod

def main():
    with open("rosalind_cunr.txt", "r") as f:
        n = int(f.readline().strip())
    print(prod(range((2*n - 5), 0, -2)) % 1_000_000)

if __name__ == "__main__":
    main()