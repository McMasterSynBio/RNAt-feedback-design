import nupack, numpy as np
import argparse, sys, os

from utils._kozak import locate_kozak
from data.kozak_pattern import KOZAK

def _kozak_exposure_probability_temp(seq: str, T: float, kozak_pattern: str = KOZAK) -> float:
    """
    Returns the probability that the Kozak sequence is exposed at given temperature input.
    """
    model = nupack.Model(celsius=T, material='rna')
    pm = nupack.pairs(strands=[seq], model=model)
    mat = pm.to_array()
    # Unpack binding probabilities
    kozak_pos = locate_kozak(seq, kozak_pattern)
    if kozak_pos is None:
        raise ValueError("Kozak sequence not found in the input sequence.")
    ks, ke = kozak_pos
    # Average exposure probability along with bottleneck (min)
    probs = np.diag(mat)[ks:ke]
    avg, bn = np.mean(probs), np.min(probs)
    return avg, bn



if __name__ == "__main__":
    """
    Testing module for calculating the Kozak exposure probability at a given temperature.
    """

    parser = argparse.ArgumentParser(description="Calculate the Kozak exposure probability at a given temperature.")
    parser.add_argument("sequence", type=str, help="The RNA sequence to analyze.")
    parser.add_argument("--temp", type=float, default=32.0, help="The temperature (°C) at which to calculate the Kozak exposure probability.")
    args = parser.parse_args()

    T = args.temp
    avg, bn = _kozak_exposure_probability_temp(args.sequence, T)
    print(f"Kozak exposure probability at {T}°C: avg={avg:.4f}, bn={bn:.4f}")