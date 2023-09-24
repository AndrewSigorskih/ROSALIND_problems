from itertools import combinations
from typing import List, Tuple, Set

def posToBits(leaves: Tuple[int]) -> int:
    res = 0
    for i, l in enumerate(leaves):
        if not l:
            continue
        mask = 1 << i
        res |= mask
    return res
        
def bitsToPos(bits: int, n: int) -> List[int]:
    res = [0 for _ in range(n)]
    for i in range(n):
        mask = 1 << i
        if (bits & mask):
            res[i] = 1
    return res

def get_all_quartets(branch: List[int]) -> Set:
        result = set()
        ones = [i for i, leaf in enumerate(branch) if leaf]
        zeros = [i for i, leaf in enumerate(branch) if not leaf]
        if (len(ones) < 2) or (len(zeros) < 2):
            return result
        for i, j in combinations(ones, 2):
            for k,l in combinations(zeros, 2):
                result.add(
                    frozenset(
                        [frozenset([i,j]), frozenset([k,l])]
                    )
                )
        return result

class Branch:
    __slots__ = ('_branch', '_size', '_leavesNum')
    def __init__(self, leaves: Tuple[int]):
        self._branch = posToBits(leaves)
        self._size = len(leaves)
        self._leavesNum = sum(leaves)

    def __str__(self) -> str:
        return "".join((str(x) for x in self.to_list()))

    @property
    def branch(self) -> int:
        return self._branch
    
    @property
    def size(self) -> int:
        return self._size
    
    @property
    def leavesNum(self) -> int:
        return self._leavesNum
    
    def to_list(self) -> List[int]:
        return bitsToPos(self._branch, self._size)
