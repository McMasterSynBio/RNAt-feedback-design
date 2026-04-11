import os, sys, argparse

def _au_gc_comp(seq: str) -> int:
    """Returns the AU/GC composition of a given sequence."""
    if len(seq) == 0:
        return 0
    
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Calculate the AU/GC composition of a given sequence."
    )
    parser.add_argument("sequence", type=str, help="The sequence to calculate the AU/GC composition for.")
    args = parser.parse_args()

    au_gc_comp = _au_gc_comp(args.sequence)
    print(f"AU/GC composition: {au_gc_comp}")
    