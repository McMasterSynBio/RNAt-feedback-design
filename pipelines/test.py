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

    KOZAK = "AUGG"
    df = pd.read_csv(args.input_dir)

    # Evaluate each sequence against the constraints
    os.makedirs(f"./results/test/{csv_filename}", exist_ok=True)
    with open(f"./results/test/{csv_filename}/results.txt", "w") as fh:
        for row in df.itertuples(index=False):
            name = row.name
            seq = str(row.sequence).upper()
            kozak = ""
            try:
                kozak = str(row.kozak).upper()
            except AttributeError:
                kozak = KOZAK
                if KOZAK not in seq:
                    seq = seq + KOZAK

            print(f"Sequence: {name}", file=fh)
            print(f"RNA     : {seq}", file=fh)
            print(f"Kozak    : {kozak}", file=fh)
            print("-" * 50, file=fh)
            window = ExposureWindow(t_lo=args.T_lo, t_hi=args.T_hi, t_step=args.T_step)
            probs_avg, probs_bn = exposure_window_prob(seq, exposure_params=window, kozak_pattern=kozak)
            print(f"{'Temp:':<22}" + ", ".join([f"{t:10.1f}" for t, _ in probs_avg.items()]), file=fh)
            print(f"{'Exposure Avgv:':<22}" + ", ".join([f"{float(n):10.6f}" for _, n in probs_avg.items()]), file=fh)
            print(f"{'Exposure Bottleneck:':<22}" + ", ".join([f"{float(n):10.6f}" for _, n in probs_bn.items()]), file=fh)
            print("-" * 50, file=fh)
            print(file=fh)