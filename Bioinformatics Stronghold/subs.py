with open('rosalind_subs.txt') as infile:
    target = infile.readline().strip('\n')
    query = infile.readline().strip('\n')
pos = 0
found = ""
while pos < len(target):
    pos = target.find(query, pos)
    if pos >= 0:
        found += (str(pos + 1) + " ")
    if pos == -1:
            break
    pos += 1
with open('out.txt', 'w') as out:
    out.write(found.strip(" "))
