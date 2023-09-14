from math import comb

def WrightFisher_drift(N: int, m: int, g: int, k: int) -> float:
    p = m/(2*N) # dominant allele
    q = 1 - p   # recessive_allele

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

    return sum(prob[k-1:])

def main():
    with open("rosalind_wfmd.txt", "r") as f:
        N, m, g, k = map(int, f.readline().split())
    with open("out.txt", "w") as o:
        print(f"{WrightFisher_drift(N, m, g, k):.3f}", file=o)
if __name__ == "__main__":
    main()