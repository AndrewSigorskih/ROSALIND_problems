masses = { 
        'A' : 71.03711,
        'C' : 103.00919,
        'D' : 115.02694,
        'E' : 129.04259,
        'F' : 147.06841,
        'G' : 57.02146,
        'H' : 137.05891,
        'I' : 113.08406,
        'K' : 128.09496,
        'L' : 113.08406,
        'M' : 131.04049,            
        'N' : 114.04293,
        'P' : 97.05276,
        'Q' : 128.05858,
        'R' : 156.10111,
        'S' : 87.03203,
        'T' : 101.04768,
        'V' : 99.06841,
        'W' : 186.07931,
        'Y' : 163.06333 
        }

def get_key(inp, i):
    res_key, res_val = min(masses.items(), \
            key=lambda x: abs((inp[i] - inp[i-1]) - x[1]))
    return res_key

def main():
    with open("rosalind_spec.txt", "r") as f:
        inp = tuple(map(float, (x.strip() for x in f)))
    with open("out.txt", "w") as o:
        print(*(get_key(inp, i) for i in range(1, len(inp))), sep="", file=o)
            
if __name__ == "__main__":
    main()