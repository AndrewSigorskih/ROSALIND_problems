# Python implementation to
# find longest increasing
# subsequence
# O(n Log n) time.
# Credit: Anant Agarwal, 
# https://www.geeksforgeeks.org/construction-of-longest-monotonically-increasing-subsequence-n-log-n/amp/

# Binary search
def GetCeilIndex(arr, T, l, r, key):
    while (r - l > 1):
        m = l + (r - l)//2
        if (arr[T[m]] >= key):
            r = m
        else:
            l = m
    return r

def LongestIncreasingSubsequence(arr, n):
    tailIndices =[0 for i in range(n + 1)]  
    prevIndices =[-1 for i in range(n + 1)]  
    len = 1
    for i in range(1, n):
        if (arr[i] < arr[tailIndices[0]]):
            tailIndices[0] = i
        elif (arr[i] > arr[tailIndices[len-1]]):
            prevIndices[i] = tailIndices[len-1]
            tailIndices[len] = i
            len += 1
        else:
            pos = GetCeilIndex(arr, tailIndices, -1,
                                   len-1, arr[i])
            prevIndices[i] = tailIndices[pos-1]
            tailIndices[pos] = i

    res = []
    i = tailIndices[len-1]
    while(i >= 0):
        res.append(arr[i])
        i = prevIndices[i]
    return res[::-1]

def main():
    with open("rosalind_lgis.txt", "r") as f:
        n = int(f.readline().strip())
        arr = list(map(int, f.readline().strip().split()))
    with open("out.txt", "w") as o:
        res = LongestIncreasingSubsequence(arr, n)
        print(*res, sep=" ", file=o)
        res = LongestIncreasingSubsequence([-x for x in arr], n)
        print(*(-x for x in res), sep=" ", file=o)
    
if __name__ == "__main__":
    main()