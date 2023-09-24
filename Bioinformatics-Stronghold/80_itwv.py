import numpy as np
from collections import Counter

# recursive function to check if motifs are interwoven
def is_interwoven(s1:str, s2:str, text_substr: str) -> bool:
    if (len(text_substr) == 0):
        return True
    if ((s1[0] == text_substr[0]) and (s2[0] == text_substr[0])):
        return is_interwoven(s1[1:], s2, text_substr[1:]) or is_interwoven(s1, s2[1:], text_substr[1:])
    elif (s1[0] == text_substr[0]):
        return is_interwoven(s1[1:], s2, text_substr[1:])
    elif (s2[0] == text_substr[0]):
        return is_interwoven(s1, s2[1:], text_substr[1:])
    return False

def itwv_matrix(text:str, dna: list) -> np.array:
    mat = np.zeros((len(dna), len(dna)), dtype=int)
    for i in range(len(dna)):
        for j in range(i, len(dna)):
            s1 = dna[i]
            s2 = dna[j]
            tot_len = len(s1) + len(s2)
            # only consider substrings that have same nt composition
            counts = Counter(s1 + s2)
            for index in range(len(text)-tot_len+1):
                substr = text[index:index+tot_len]
                substr_counts = Counter(substr)
                # adding extra character to stop recursion at and of motifs
                if (counts == substr_counts) and is_interwoven(s1+"$", s2+"$", substr):
                    mat[i][j] = 1
                    mat[j][i] = 1
                    break
    return mat

def main():
    with open("rosalind_itwv.txt", "r") as f:
        text = f.readline().strip()
        dna = [line.strip() for line in f]
    result = itwv_matrix(text, dna)
    with open("out.txt", "w") as o:
        for line in result:
            print(*line, file=o, sep=' ')

if __name__ == "__main__":
    main()