import os

import numpy as np
import pandas as pd

from constraints import CompositionConstraint, LoopConstraint, StemConstraint, KozakExposureConstraint
from .helpers.mean_exposure import _exposure_window_weight

def test_sequences_from_csv():
    KOZAK = "AUGG"
    # Load test sequences from CSV
    df = pd.read_csv("./data/sequences.csv")
    records = []
    for row in df.itertuples(index=False):
        seq = str(row.sequence).upper()
        seq = seq + KOZAK if KOZAK not in seq else seq
        records.append((str(row.name), seq))

    # Initialize constraints
    comp_constraint = CompositionConstraint()
    loop_constraint = LoopConstraint()
    stem_constraint = StemConstraint()
    exposure_constraints = [KozakExposureConstraint(target=Tm) for Tm in np.arange(30, 50, 0.5)]

    # Evaluate each sequence against the constraints
    os.makedirs("./results/test", exist_ok=True)
    with open("./results/test/results.txt", "w") as fh:
        for name, seq in records:
            comp_value = comp_constraint.evaluate((seq,), None).value
            loop_value = loop_constraint.evaluate((seq,), None).value
            stem_value = stem_constraint.evaluate((seq,), None).value

            print(f"Sequence: {name}", file=fh)
            print(f"RNA     : {seq}", file=fh)
            print(
                f"Base metrics -> GC/AU: {comp_value:.4f} | Loop size: {loop_value} | Stem length: {stem_value}",
                file=fh,
            )
            print(
                "target_C | exposure_error | mean_exposure_below | mean_exposure_above",
                file=fh,
            )
            print("-" * 68, file=fh)
            for exposure_constraint in exposure_constraints:
                result = exposure_constraint.evaluate((seq,), None)
                mean_below, mean_above = _exposure_window_weight(seq, exposure_constraint)
                below_str = "NA" if mean_below is None else f"{mean_below:.4f}"
                above_str = "NA" if mean_above is None else f"{mean_above:.4f}"
                print(
                    f"{exposure_constraint.target:8.1f} | "
                    f"{result.excess:14.4f} | "
                    f"{below_str:19} | "
                    f"{above_str:19}",
                    file=fh,
                )
            print("-" * 50, file=fh)
            print(file=fh)

if __name__ == "__main__":
    test_sequences_from_csv()