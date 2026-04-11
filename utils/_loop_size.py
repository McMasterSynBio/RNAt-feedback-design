import os, sys, argparse

def _loop_size(seq: str) -> int:
    """Returns the size of the loop in a given sequence."""
    return 0

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the loop size of a given sequence."
    )
    parser.add_argument("sequence", type=str, help="The sequence to calculate the loop size for.")
    args = parser.parse_args()

    loop_size = _loop_size(args.sequence)
    print(f"Loop size: {loop_size}")
    