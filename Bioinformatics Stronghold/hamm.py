with open('rosalind_hamm.txt') as infile:
    str1 = infile.readline().strip('\n')
    str2 = infile.readline().strip('\n')
distance = 0
for base1, base2 in zip(str1, str2):
    if base1 != base2:
        distance += 1
print(distance)
