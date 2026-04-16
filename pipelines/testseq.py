import os, sys, numpy as np, pandas as pd

from constraints import CompositionConstraint, LoopConstraint, StemConstraint, KozakExposureConstraint
from .helpers.mean_exposure import _exposure_window_means

def test_sequences_from_csv():
    
    KOZAK = 'AUGG'
    # Load test sequences from CSV
    df = pd.read_csv("./data/sequences.csv")
    sequences = [str(seq).upper() for seq in df["sequence"].tolist()]
    sequences = [(seq + KOZAK if KOZAK not in seq else seq) for seq in sequences]

    # Initialize constraints
    comp_constraint = CompositionConstraint()
    loop_constraint = LoopConstraint()
    stem_constraint = StemConstraint()
    exposure_constraints = [KozakExposureConstraint(target=Tm) for Tm in np.arange(30, 50, 0.5)]

    # Evaluate each sequence against the constraints
    os.makedirs("./results/test", exist_ok=True)
    with open("./results/test/results.txt", "w") as fh:
        for seq in sequences:
            print(f"Testing sequence: {seq}", file=fh)
            print(f"GC/AU: {comp_constraint.evaluate((seq,), None).value} | \
                    Loop size: {loop_constraint.evaluate((seq,), None).value} | \
                    Stem length: {stem_constraint.evaluate((seq,), None).value}"
                    , file=fh
            )
            for exposure_constraint in exposure_constraints:
                print(exposure_constraint.evaluate((seq,), None).value, file=fh)
                mean_below, mean_above = _exposure_window_means(seq, exposure_constraint)
                print(f"Mean exposure below target: {mean_below}", file=fh)
                print(f"Mean exposure above target: {mean_above}", file=fh)
            print("-" * 50, file=fh)

if __name__ == "__main__":
    test_sequences_from_csv()