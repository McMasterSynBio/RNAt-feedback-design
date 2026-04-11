import os, sys, argparse
from ._stem_length import _stem_length

def _loop_size(seq: str) -> int:
    """Returns the size of the loop in a given sequence."""
    bps, n = _stem_length(seq)
    return bps[-1][1] - bps[-1][0] - 1 if len(bps) > 0 else 10**9 # Some v large value arbitrarily

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the loop size of a given sequence."
    )
    parser.add_argument("sequence", type=str, help="The sequence to calculate the loop size for.")
    args = parser.parse_args()

    loop_size = _loop_size(args.sequence)
    print(f"Loop size: {loop_size}")
    