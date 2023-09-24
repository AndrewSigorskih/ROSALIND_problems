def main():
    with open("rosalind_cstr.txt", "r") as f:
        lines = [x.strip() for x in f]
    with open("out.txt", "w") as o:
        for col in zip(*lines):
            mask = [1 if x == col[0] else 0 for x in col]
            if (sum(mask) > 1) and (sum(mask) < len(lines)-1):
                print(*mask, sep='', file=o)

if __name__ == "__main__":
    main()