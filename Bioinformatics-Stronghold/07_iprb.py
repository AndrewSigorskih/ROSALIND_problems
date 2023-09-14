def firstLaw(k, m, n):
    N = float(k + m + n)
    return 1 - ( m*n + .25*m*(m-1) + n*(n-1) ) / ( N*(N-1) )

if __name__ == "__main__":
    with open("rosalind_iprb.txt", "r") as f:
        k, m, n = map(int, f.readline().strip().split())
    print(f"{firstLaw(k, m, n):.5f}")