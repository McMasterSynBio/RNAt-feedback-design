"""RNA thermosensor sequence design pipeline.

Designs RNA sequences that melt at a target temperature while satisfying
structural constraints (GC composition, loop size, stem length).

Usage examples:

    # Quick 3-run design targeting Tm = 37°C (default)
    python main.py --num-runs 3

    # Custom target Tm, sequence length, 20 runs for diversity
    python main.py --target-tm 42 --flank-seq-length 40 --num-runs 20
"""

import argparse, os, sys
from random import randint

from pipelines import run_design_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Design RNA thermosensor sequences with a target melting temperature."
    )

    # Core design parameters
    parser.add_argument("--flank-one-length",   type=int,   default=20,    help="Flank sequence length (nt)")
    parser.add_argument("--flank-two-length",   type=int,   default=5,    help="Flank sequence length (nt)")
    parser.add_argument("--kozak-seq", type=str, default="AAAAAA", help="Kozak sequence to include in the design")
    parser.add_argument("--tail-seq", type=str, default="GCUUCAGGU", help="Tail sequence to include in the design")

    parser.add_argument("--num-runs",     type=int,   default=10,    help="Number of independent NUAD runs")
    parser.add_argument("--target-tm",    type=float, default=37.0,  help="Target melting temperature (°C)")

    # Search control
    parser.add_argument("--max-iterations", type=int, default=50,  help="Max NUAD iterations per run")
    parser.add_argument("--out-dir",      type=str,   default="results/opt", help="Output directory")

    args = parser.parse_args()

    # Safety assertions
    assert 8 < args.flank_one_length < 31, "Flank sequence length must be between 9 and 30 nt."
    assert 0 < args.flank_two_length < 11, "Flank sequence length must be between 1 and 10 nt."

    # dirs within /results
    n_sim = sum(1 for entry in os.scandir(args.out_dir) if entry.is_dir())
    out_root = os.path.join(args.out_dir, f"sim_{n_sim+1:02d}")
    os.makedirs(out_root, exist_ok=True)

    # generate random seed
    random_seed = randint(0, 2**16 - 1)


    run_design_pipeline(
        kozak=args.kozak_seq,
        tail=args.tail_seq,
        flank_one_length=args.flank_one_length,
        flank_two_length=args.flank_two_length,
        num_runs=args.num_runs,
        target_tm=args.target_tm,
        max_iterations=args.max_iterations,
        out_directory=out_root,
        base_random_seed=random_seed,
    )


if __name__ == "__main__":
    main()

