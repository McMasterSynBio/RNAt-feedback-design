"""Melting temperature (Tm) estimation via NUPACK thermodynamics.

Computes the temperature at which the MFE hairpin structure's free energy
crosses zero — i.e. where the folded secondary structure ceases to be
thermodynamically favourable (the melting transition).
"""

import nupack
from scipy.optimize import brentq


def _structure_energy_at(seq: str, structure: nupack.Structure, celsius: float) -> float:
    """Return free energy (kcal/mol) of a specific *structure* at *celsius* °C."""
    strand = nupack.Strand(seq, name="probe")
    model = nupack.Model(celsius=celsius, material="rna")
    return float(nupack.structure_energy([strand], structure, model=model))


def estimate_tm(
    seq: str,
    t_lo: float = 10.0,
    t_hi: float = 90.0,
    reference_celsius: float = 25.0,
) -> float:
    """Estimate the melting temperature of *seq* by bisection on ΔG_structure(T).

    1. Finds the MFE structure at *reference_celsius* (the "folded" state).
    2. Sweeps `structure_energy` of that fold across temperatures.
    3. Returns the temperature where ΔG_structure crosses zero.

    If ΔG_structure is positive at *t_lo* (never stable), returns *t_lo*.
    If ΔG_structure is negative at *t_hi* (very stable), returns *t_hi*.
    """
    strand = nupack.Strand(seq, name="probe")
    ref_model = nupack.Model(celsius=reference_celsius, material="rna")
    mfe_results = nupack.mfe([strand], model=ref_model)
    mfe_struct = mfe_results[0].structure

    dG_lo = _structure_energy_at(seq, mfe_struct, t_lo)
    dG_hi = _structure_energy_at(seq, mfe_struct, t_hi)

    # No sign change → Tm is outside the search window
    if dG_lo >= 0:
        return t_lo
    if dG_hi <= 0:
        return t_hi

    # Final structure energy at mfe optimal, varying temperatures
    # Lowest energy value (most stable) is likely most relevant to the structure.
    tm = brentq(
        lambda t: _structure_energy_at(seq, mfe_struct, t),
        t_lo, t_hi, xtol=0.1,
    )
    return float(tm)
