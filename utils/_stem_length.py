import os, sys, argparse

def _stem_length(seq: str) -> int:
    """Returns the length of the stem in a given sequence."""
    return 0

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the stem length of a given sequence."
    )
    parser.add_argument("sequence", type=str, help="The sequence to calculate the stem length for.")
    args = parser.parse_args()

    stem_length = _stem_length(args.sequence)
    print(f"Stem length: {stem_length}")
    