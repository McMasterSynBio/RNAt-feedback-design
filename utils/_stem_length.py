import os, sys, argparse
import numpy as np

PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}


def _traceback_pairs(backtrack, i: int, j: int) -> list[tuple[int, int]]:
    if i >= j:
        return []

    kind, value = backtrack[i][j]
    if kind == "skip":
        return _traceback_pairs(backtrack, i, j - 1)
    if kind == "pair":
        t = value
        return (
            _traceback_pairs(backtrack, i, t - 1)
            + [(t, j)]
            + _traceback_pairs(backtrack, t + 1, j - 1)
        )
    return []

def _stem_length(seq: str) -> tuple[list[tuple[int,int]], int]:
    """Returns the length of the stem in a given sequence."""

    assert len(seq) > 4, "Sequence must be at least 5 nucleotides long to have a stem."

    n = len(seq)
    opt = np.zeros((n+1, n+1), dtype=int)
    backtrack = [[("done", 0) for _ in range(n+1)] for _ in range(n+1)]

    for k in range (4, n):
        for i in range(1, n - k + 1):
            j = i + k
            opt_b = opt[i, j - 1]
            backtrack[i][j] = ("skip", 0)
            for t in range(i, j - 3):
                if (seq[t-1], seq[j-1]) in PAIRS:
                    if opt[i, t-1] + 1 + opt[t + 1, j - 1] > opt_b:
                        backtrack[i][j] = ("pair", t)
                        opt_b = opt[i, t-1] + 1 + opt[t + 1, j - 1]         
            opt[i, j] = opt_b

    return _traceback_pairs(backtrack, 1, n), opt[1, n]

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the stem length of a given sequence."
    )
    parser.add_argument("sequence", type=str, help="The sequence to calculate the stem length for.")
    args = parser.parse_args()

    pairs, k = _stem_length(args.sequence)
    print(f"Stem length: {k} with pairs {pairs}")
    