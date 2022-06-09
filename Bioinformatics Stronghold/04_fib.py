def fib(n, k):
    a, b = 1, 1
    for i in range(2, n):
        a, b = b, k*a + b
    return b

if __name__ == "__main__":
    with open("rosalind_fib.txt", "r") as f:
        n, k = map(int, f.readline().strip().split())
    print(fib(n,k))