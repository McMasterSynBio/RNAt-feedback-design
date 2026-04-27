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
    parser.add_argument("--flank-seq-length",   type=int,   default=30,    help="Flank sequence length (nt)")
    parser.add_argument("--num-runs",     type=int,   default=10,    help="Number of independent NUAD runs")
    parser.add_argument("--target-tm",    type=float, default=37.0,  help="Target melting temperature (°C)")

    # Search control
    parser.add_argument("--max-iterations", type=int, default=50,  help="Max NUAD iterations per run")
    parser.add_argument("--out-dir",      type=str,   default="results/design", help="Output directory")

    args = parser.parse_args()

    # Safety assertions
    assert 8 < args.flank_seq_length < 31, "Flank sequence length must be between 9 and 30 nt."

    # dirs within /results
    n_sim = sum(1 for entry in os.scandir(args.out_dir) if entry.is_dir())
    out_root = os.path.join(args.out_dir, f"sim_{n_sim+1:02d}")
    os.makedirs(out_root, exist_ok=True)

    # generate random seed
    random_seed = randint(0, 2**16 - 1)


    run_design_pipeline(
        flank_length=args.flank_seq_length,
        num_runs=args.num_runs,
        target_tm=args.target_tm,
        max_iterations=args.max_iterations,
        out_directory=out_root,
        base_random_seed=random_seed,
    )


if __name__ == "__main__":
    main()

