import numpy as np
from math import comb
from itertools import product

def WrightFisher_drift(N: int, m: int, g: int) -> float:
    q = m/(2*N) # dominant allele
    p = 1 - q   # recessive_allele

    # probability of exactly i copies of recessive allele in gen 1
    prob = [comb(2*N, i) * (q**(i)) * (p**(2*N-i)) for i in range(1, 2*N+1)]

    # probability of exactly t copies of recessive allele in next gen
    for gen in range(1, g):
        gen_prob = []
        # probability of exactly t copies of recessive allele in current gen
        for t in range(1, 2*N+1):
            prob_t = [comb(2*N, t) * ((i/(2*N))**(t)) * ((1-(i/(2*N)))**(2*N-t)) for i in range(1, 2*N+1)]
            gen_prob.append(sum([prob_t[j] * prob[j] for j in range(2*N)]))
        prob = gen_prob

    return np.log10(1-sum(prob))

def foun(N: int, m: int, A: list) ->np.array:
    k = len(A)
    B = np.zeros((m, k))
    for i, j in product(range(m), range(k)):
        B[i,j] = WrightFisher_drift(N, A[j], i+1)
    return B

def main():
    with open("rosalind_foun.txt", "r") as f:
        N, m = map(int, f.readline().split())
        A = list(map(int, f.readline().split()))
    res = foun(N, m, A)
    with open("out.txt", "w") as o:
        for line in res:
            print(*line, sep=' ', file=o)
if __name__ == "__main__":
    main()