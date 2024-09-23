import numpy as np
import copy
import argparse


def get_pos_point(a: str, i: int, b: str, j: int, match: int, mismatch: int):
    return match if a[i] == b[j] else mismatch


def get_alignments(M, a: str, i: int, b: str, j: int, gap: int, match: int, mismatch: int, mem):
    if i == 0 and j == 0:
        return [["", ""]]
    
    alignments = []
    if i > 0 and j > 0 and M[i-1, j-1] + get_pos_point(a, i - 1, b, j - 1, match, mismatch) == M[i, j]:
        if not (i - 1, j - 1) in mem:
            mem[(i - 1, j - 1)] = get_alignments(M, a, i - 1, b, j - 1, gap, match, mismatch, mem)
        
        aux = copy.deepcopy(mem[(i - 1, j - 1)])
        for k in range(len(aux)):
            aux[k][0] = aux[k][0] + a[i - 1]
            aux[k][1] = aux[k][1] + b[j - 1]
        
        alignments.extend(aux)
    if i > 0 and M[i - 1, j] + gap == M[i, j]:
        if not (i - 1, j) in mem:
            mem[(i - 1, j)] = get_alignments(M, a, i - 1, b, j, gap, match, mismatch, mem)
        
        aux = copy.deepcopy(mem[(i - 1, j)])
        for k in range(len(aux)):
            aux[k][0] = aux[k][0] + a[i - 1]
            aux[k][1] = aux[k][1] + "-"
        
        alignments.extend(aux)
    if j > 0 and M[i, j - 1] + gap == M[i, j]:
        if not (i, j - 1) in mem:
            mem[(i, j - 1)] = get_alignments(M, a, i, b, j - 1, gap, match, mismatch, mem)
        
        aux = copy.deepcopy(mem[(i, j - 1)])
        for k in range(len(aux)):
            aux[k][0] = aux[k][0] + "-"
            aux[k][1] = aux[k][1] + b[j - 1]
        
        alignments.extend(aux)
    
    return alignments


def dp_align(a: str, m: int, b: str, n: int, gap: int, match: int, mismatch: int):
    scores = np.zeros((m + 1, n + 1), dtype=np.int16)
    scores[1:m + 1, 0] = gap * np.arange(1, m + 1)
    scores[0, 1:n + 1] = gap * np.arange(1, n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            scores[i, j] = max(scores[i, j - 1] + gap,
                               scores[i - 1, j - 1] + get_pos_point(a, i - 1, b, j - 1, match, mismatch),
                               scores[i - 1, j] + gap)
    
    return [scores[m, n], get_alignments(scores, a, m, b, n, gap, match, mismatch, {})]


GAP = -5
MATCH = 3
MISMATCH = -2


parser = argparse.ArgumentParser()
parser.add_argument("a", help="First sequence to be aligned", type=str)
parser.add_argument("b", help="Second sequence to be aligned", type=str)
parser.add_argument("--gap", default=GAP, help="Gap score", type=int)
parser.add_argument("--match", default=MATCH, help="Match score", type=int)
parser.add_argument("--mismatch", default=MISMATCH, help="Mismatch score", type=int)
args = parser.parse_args()

result = dp_align(args.a, len(args.a), args.b, len(args.b), args.gap, args.match, args.mismatch)
print(f"Score - {result[0]}\nAlignments - {result[1]}")