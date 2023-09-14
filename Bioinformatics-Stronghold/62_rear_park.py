# Park’s Exact Greedy Algorithm
# Source: https://medium.com/@matthewwestmk/87c62d690eef
def get_reverse(seq: list, i: int, j: int) -> list:
    return seq[:i] + seq[i:j][::-1] + seq[j:]

def findBreakpoints(seq:list, targ: list) -> list:
    return [
        i + 1 for i in range(len(seq) - 1) \
            if (abs(targ.index(seq[i]) - targ.index(seq[i + 1])) != 1)
    ]

def findMinimumBreakpointReversals(seqs: list, targ: list) -> list:
    reversals = []
    for seq in seqs:
        breakpoints = findBreakpoints(seq[0], targ)
        for i in range(len(breakpoints)-1):
            for j in range(i+1, len(breakpoints)):
                reversals.append((get_reverse(seq[0], breakpoints[i], breakpoints[j]), \
                                  seq[1] + [(breakpoints[i]-1, breakpoints[j]-1)])
                                )
    min_bp = len(targ)
    min_reversals = []
    for reversal in reversals:
        num_breakpoints = len(findBreakpoints(reversal[0], targ))
        if num_breakpoints < min_bp:
            min_bp = num_breakpoints
            min_reversals = [reversal]
        elif num_breakpoints == min_bp:
            min_reversals.append(reversal)
    return min_reversals

def sort_by_reversals(seq: list, targ: list) -> tuple:
    count = 0
    seq = ["-"] + seq + ["+"]
    targ = ["-"] + targ + ["+"]
    cur_seqs = [(seq, [])]
    while targ not in [current[0] for current in cur_seqs]:
        cur_seqs = findMinimumBreakpointReversals(cur_seqs, targ)
        count += 1
    return count, cur_seqs

def main():
    with open("rosalind_rear.txt", "r") as f:
        lines = [list(map(int, line.strip().split(" "))) \
                for line in f.readlines() \
                    if len(line) > 2]
    result = []
    for i in range(0, len(lines), 2):
        result.append(sort_by_reversals(lines[i], lines[i+1])[0])
    
    with open("out.txt", "w") as o:
        print(*result, file=o)

if __name__ == "__main__":
    main()