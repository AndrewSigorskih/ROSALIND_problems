from Bio import SeqIO

def hamming_dist(s1, s2):
    if (len(s1) != len(s2)):
        return -1
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def correct_reads(records: list)->list:
    good_reads = set()
    bad_reads = list()
    result = []
    # split groups
    for i in range(len(records)):
        s = records[i]
        s_rc = s.reverse_complement()
        if (s in good_reads) or s_rc in good_reads:
            continue
        for j in range(i+1, len(records)):
            if (hamming_dist(s, records[j]) == 0) or \
                (hamming_dist(s_rc, records[j]) == 0):
                good_reads.add(s)
                good_reads.add(s_rc)
                break
        if (s not in good_reads):
            bad_reads.append(s)
    # process bad reads
    for item in bad_reads:
        for read in good_reads:
            if (hamming_dist(item, read) == 1):
                result.append(f"{item}->{read}")
    return result

def main():
    records = [x.seq for x in SeqIO.parse("rosalind_corr.txt", "fasta")]    
    res = correct_reads(records)
    with open("out.txt", "w") as o:
        print(*res, sep='\n', file=o)

if __name__ == "__main__":
    main()