from peasyprofiller.profiller import profiller as pprof
from random import choice
import numpy as np
import argparse


def get_pos_point(a: str, i: int, b: str, j: int, match: int, mismatch: int):
    return match if a[i] == b[j] else mismatch


def brute_force_align(a: str, ai: int, af: int, b: str, bi: int, bf: int, gap: int, match: int, mismatch: int):
    if ai > af:
        return gap * (bf - bi + 1)
    if bi > bf:
        return gap * (af - ai + 1)
    
    max_points = gap + brute_force_align(a, ai + 1, af, b, bi, bf, gap, match, mismatch)
    for bk in range(bi, bf + 1):
        aux = brute_force_align(a, ai + 1, af, b, bk + 1, bf, gap, match, mismatch)
        aux1 = gap * (bk - bi) + get_pos_point(a, ai, b, bk, match, mismatch) + aux
        aux2 = gap * (bk - bi + 2) + aux
        max_points = max(max_points, aux1, aux2)
    return max_points


def simple_brute_force_align(a: str, m: int, b: str, n: int, gap: int, match: int, mismatch: int):
    if m == 0:
        return gap * n
    if n == 0:
        return gap * m
    return max(simple_brute_force_align(a, m - 1, b, n - 1, gap, match, mismatch) + get_pos_point(a, m - 1, b, n - 1, match, mismatch),
               simple_brute_force_align(a, m - 1, b, n, gap, match, mismatch) + gap,
               simple_brute_force_align(a, m, b, n - 1, gap, match, mismatch) + gap)


def dp_align(a: str, m: int, b: str, n: int, gap: int, match: int, mismatch: int):
    scores = np.zeros((m + 1, n + 1), dtype=np.int16)
    scores[1:m + 1, 0] = gap * np.arange(1, m + 1)
    scores[0, 1:n + 1] = gap * np.arange(1, n + 1)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            scores[i, j] = max(scores[i, j - 1] + gap,
                               scores[i - 1, j - 1] + get_pos_point(a, i - 1, b, j - 1, match, mismatch),
                               scores[i - 1, j] + gap)
    
    return scores[m, n]


pprof.start("Sequence generation")
NUCLEOTIDES = ["A", "T", "C", "G"]
GAP = -5
MATCH = 3
MISMATCH = -2

parser = argparse.ArgumentParser()
parser.add_argument("n", help="Size of random sequences to be generated", type=int)
args = parser.parse_args()

a = "".join(choice(NUCLEOTIDES) for i in range(args.n))
b = "".join(choice(NUCLEOTIDES) for i in range(args.n))
pprof.stop("Sequence generation")

if args.n <= 12:
    pprof.start("Brute")
    print(brute_force_align(a, 0, args.n - 1, b, 0, args.n - 1, GAP, MATCH, MISMATCH))
    pprof.stop("Brute")
else:
    print("Skipping brute force (n > 12)")

if args.n <= 10:
    pprof.start("Simple brute")
    print(simple_brute_force_align(a, args.n, b, args.n, GAP, MATCH, MISMATCH))
    pprof.stop("Simple brute")
else:
    print("Skipping simple brute force (n > 10)")

pprof.start("DP")
print(dp_align(a, args.n, b, args.n, GAP, MATCH, MISMATCH))
pprof.stop("DP")

save_path = f"results/global_{args.n}"
pprof.save_csv(save_path)
pprof.plot(save_path)