# "Meet-in-the-middle"-like approach
def get_all_reversals(s: list) -> list:
    # get all possible single reversals of array
    result = []
    for i in range(len(s)-1):
        for j in range(i+1, len(s)):
            r_list = s[i:j+1]
            result.append(s[:i] + list(reversed(r_list)) + s[j+1:])
    return result

def reversal_distance(s1: set, s2: set, distance: int) -> int:
    if s1 & s2: # if al least one reversal in common
        return distance
    new_s1 = set()

    for s in s1:
        reverse_arrays = get_all_reversals(list(s))
        for r in reverse_arrays:
            new_s1.add(tuple(r))

    new_s2 = set()
    for s in s2:
        reverse_arrays = get_all_reversals(list(s))
        for r in reverse_arrays:
            new_s2.add(tuple(r))

    distance += 2 # both arrays reversed

    if (s1 & new_s2) or (s2 & new_s1):
        return distance - 1
    
    if new_s1 & new_s2:
        return distance
    
    distance = reversal_distance(new_s1, new_s2, distance)
    return distance

def main():
    with open("rosalind_rear.txt", "r") as f:
        lines = [line.strip().split(" ") for line in f.readlines() if len(line) > 2]
    result = []
    for i in range(0, len(lines), 2):
        l1 = [int(l) for l in lines[i]]
        l2 = [int(l) for l in lines[i+1]]

        s1, s2 = set(), set()
        s1.add(tuple(l1))
        s2.add(tuple(l2))

        result.append(reversal_distance(s1, s2, 0))

    with open("out.txt", "w") as o:
        print(*result, file=o)

if __name__ == "__main__":
    main()