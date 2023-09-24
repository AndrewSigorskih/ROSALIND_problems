def pdpl(data: list) -> list:
    res = [0, max(data)]
    for item in set(data):
        diff = [abs(item - x) for x in res]
        if set(diff).issubset(set(data)):
            res.append(item)
            for d in diff:
                data.remove(d)
    return sorted(res)

def main():
    with open("rosalind_pdpl.txt", "r") as f:
        data = list(map(int, (x for x in f.readline().split())))
    with open("out.txt", "w") as o:
        print(*pdpl(data), file=o)

if __name__ == "__main__":
    main()