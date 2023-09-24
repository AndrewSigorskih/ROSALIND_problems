def N_statistic(contigs: list, lengths: list, 
                tot_len: int, N: float = .5) -> int:
    for L in lengths:
        cur_sum = sum(len(c) for c in contigs if len(c) >= L)
        if (cur_sum / tot_len >= N):
            return L

def asmq(contigs: list) -> tuple:
    lengths = [len(c) for c in contigs]
    tot_len = sum(lengths)
    lengths = sorted(set(lengths), reverse=True)
    return N_statistic(contigs, lengths, tot_len, N=.5),\
        N_statistic(contigs, lengths, tot_len, N=.75)

def main():
    with open("rosalind_asmq.txt", "r") as f:
        contigs = [line.strip() for line in f]
    res = asmq(contigs)
    with open("out.txt", "w") as o:
        print(*res, file=o)
if __name__ == "__main__":
    main()