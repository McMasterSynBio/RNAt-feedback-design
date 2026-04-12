"""RNA thermosensor sequence design pipeline.

Designs RNA sequences that melt at a target temperature while satisfying
structural constraints (GC composition, loop size, stem length).

Usage examples:

    # Quick 3-run design targeting Tm = 37°C (default)
    python main.py --num-runs 3

    # Custom target Tm, sequence length, 20 runs for diversity
    python main.py --target-tm 42 --seq-length 40 --num-runs 20
"""

import argparse

from designer import run_design_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Design RNA thermosensor sequences with a target melting temperature."
    )

    # Core design parameters
    parser.add_argument("--seq-length",   type=int,   default=30,    help="Sequence length (nt)")
    parser.add_argument("--num-runs",     type=int,   default=10,    help="Number of independent NUAD runs")
    parser.add_argument("--target-tm",    type=float, default=37.0,  help="Target melting temperature (°C)")
    parser.add_argument("--tm-tolerance", type=float, default=1.0,   help="Tm tolerance before penalty (°C)")
    parser.add_argument("--tm-weight",    type=float, default=5.0,   help="Weight of Tm constraint vs structural")

    # Structural constraint bounds
    parser.add_argument("--gc-lo",        type=float, default=0.4,   help="Min GC fraction")
    parser.add_argument("--gc-hi",        type=float, default=0.6,   help="Max GC fraction")
    parser.add_argument("--loop-lo",      type=int,   default=3,     help="Min loop size (nt)")
    parser.add_argument("--loop-hi",      type=int,   default=10,    help="Max loop size (nt)")
    parser.add_argument("--stem-lo",      type=int,   default=5,     help="Min stem length (bp)")
    parser.add_argument("--stem-hi",      type=int,   default=15,    help="Max stem length (bp)")

    # Search control
    parser.add_argument("--max-iterations", type=int, default=None,  help="Max NUAD iterations per run")
    parser.add_argument("--out-dir",      type=str,   default="results", help="Output directory")
    parser.add_argument("--random-seed",  type=int,   default=42,    help="Base random seed")

    args = parser.parse_args()

    run_design_pipeline(
        seq_length=args.seq_length,
        num_runs=args.num_runs,
        target_tm=args.target_tm,
        tm_tolerance=args.tm_tolerance,
        tm_weight=args.tm_weight,
        gc_lo=args.gc_lo,
        gc_hi=args.gc_hi,
        loop_lo=args.loop_lo,
        loop_hi=args.loop_hi,
        stem_lo=args.stem_lo,
        stem_hi=args.stem_hi,
        max_iterations=args.max_iterations,
        out_directory=args.out_dir,
        base_random_seed=args.random_seed,
    )


if __name__ == "__main__":
    main()

