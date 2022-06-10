def hamming_dist(s1, s2):
    if (len(s1) != len(s2)):
        return -1
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

if __name__ == "__main__":
    with open('rosalind_hamm.txt') as infile:
        str1 = infile.readline().strip()
        str2 = infile.readline().strip()
    print(hamming_dist(str1, str2))