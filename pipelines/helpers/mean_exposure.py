import numpy as np
from constraints.exposure import KozakExposureConstraint
from utils.therm._exposure import _kozak_exposure_probability_temp

def _exposure_window_means(
    seq: str,
    exposure_constraint: KozakExposureConstraint,
) -> tuple[float | None, float | None]:
    """
    Calculate the mean exposure probabilities below and above the target temperature.

    Args:
        seq: RNA sequence.
        exposure_constraint: KozakExposureConstraint object.

    Returns:
        A tuple containing the mean exposure probabilities below and above the target temperature.
    """
    temps = np.arange(
        exposure_constraint.t_lo,
        exposure_constraint.t_hi + exposure_constraint.t_step,
        exposure_constraint.t_step,
    )

    weighted_probs = []
    for temp in temps:
        try:
            avg, bn = _kozak_exposure_probability_temp(seq, temp)
        except ValueError:
            return None, None
        weighted_probs.append(exposure_constraint.alpha * avg + (1 - exposure_constraint.alpha) * bn)

    weighted_probs = np.array(weighted_probs)
    below = temps < exposure_constraint.target
    above = ~below

    mean_below = float(np.mean(weighted_probs[below])) if np.any(below) else None
    mean_above = float(np.mean(weighted_probs[above])) if np.any(above) else None
    return mean_below, mean_above