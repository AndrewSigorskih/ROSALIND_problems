import re

with open("rosalind_subs.txt", "r") as f:
    s, pat = (x.strip() for x in f)
idx = (m.start()+1 for m in re.finditer(f"(?={pat})", s))
print(*idx)
