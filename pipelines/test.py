import os, sys, argparse

import numpy as np
import pandas as pd

from constraints import CompositionConstraint, LoopConstraint, StemConstraint, KozakExposureConstraint
from .helpers.exposure_window import exposure_window_prob, ExposureWindow

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Get KOZAK exposure probability window average and bottlenecks for sequences."
    )

    # Core design parameters
    parser.add_argument("input_dir", type=str, default="./data/sequences.csv", help="Input directory containing sequences CSV")
    parser.add_argument("--T_lo", type=float, default=25.0, help="Lower bound of temperature range for exposure evaluation")
    parser.add_argument("--T_hi", type=float, default=75.0, help="Upper bound of temperature range for exposure evaluation")
    parser.add_argument("--T_step", type=float, default=0.5, help="Step size for temperature range for exposure evaluation")
    args = parser.parse_args()
    
    # csv filename
    csv_filename = os.path.basename(args.input_dir).split(".")[0]

    df = pd.read_csv(args.input_dir)

    # Evaluate each sequence against the constraints
    os.makedirs(f"./results/test/{csv_filename}", exist_ok=True)
    with open(f"./results/test/{csv_filename}/results.txt", "w") as fh:
        for row in df.itertuples(index=False):
            name = row.name if hasattr(row, "name") else row.rank
            seq = str(row.sequence).upper()
            kozak_start = int(row.kozak_start)
            kozak_end = int(row.kozak_end)
            try:
                right_stem_start = int(row.right_stem_start)
                right_stem_end = int(row.right_stem_end)
                start_codon_start = int(row.start_codon_start)
                start_codon_end = int(row.start_codon_end)
            except:
                right_stem_start = None
                right_stem_end = None
                start_codon_start = None
                start_codon_end = None

            print(f"Sequence: {name}", file=fh)
            print(f"RNA     : {seq}", file=fh)
            print(f"Kozak    : {seq[kozak_start:kozak_end]}", file=fh)
            if right_stem_start and start_codon_start:
                print(f"Right Stem: {seq[right_stem_start:right_stem_end]}", file=fh)
                print(f"Start Codon: {seq[start_codon_start:start_codon_end]}", file=fh)
            print("-" * 50, file=fh)
            window = ExposureWindow(t_lo=args.T_lo, t_hi=args.T_hi, t_step=args.T_step)
            kozak_probs_avg, kozak_probs_bn = exposure_window_prob(seq, exposure_params=window, pos=(kozak_start, kozak_end))
            if right_stem_start and start_codon_start:
                right_stem_probs_avg, right_stem_probs_bn = exposure_window_prob(seq, exposure_params=window, pos=(right_stem_start, right_stem_end))
                start_codon_probs_avg, start_codon_probs_bn = exposure_window_prob(seq, exposure_params=window, pos=(start_codon_start, start_codon_end))
            print(f"{'Temp:':<40}" + ", ".join([f"{t:10.1f}" for t, _ in kozak_probs_avg.items()]), file=fh)
            print(f"{'Kozak Exposure Avgv:':<40}" + ", ".join([f"{float(n):10.6f}" for _, n in kozak_probs_avg.items()]), file=fh)
            print(f"{'Kozak Exposure Bottleneck:':<40}" + ", ".join([f"{float(n):10.6f}" for _, n in kozak_probs_bn.items()]), file=fh)
            if right_stem_start and start_codon_start:
                print(f"{'Right Stem Exposure Avgv:':<40}" + ", ".join([f"{float(n):10.6f}" for _, n in right_stem_probs_avg.items()]), file=fh)
                print(f"{'Right Stem Exposure Bottleneck:':<40}" + ", ".join([f"{float(n):10.6f}" for _, n in right_stem_probs_bn.items()]), file=fh)
                print(f"{'Start Codon Exposure Avgv:':<40}" + ", ".join([f"{float(n):10.6f}" for _, n in start_codon_probs_avg.items()]), file=fh)
                print(f"{'Start Codon Exposure Bottleneck:':<40}" + ", ".join([f"{float(n):10.6f}" for _, n in start_codon_probs_bn.items()]), file=fh)
            print("-" * 50, file=fh)
            print(file=fh)