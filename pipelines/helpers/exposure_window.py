from dataclasses import dataclass

import numpy as np
from constraints.exposure import KozakExposureConstraint
from utils.therm._exposure import _kozak_exposure_probability_temp
from data.kozak_pattern import KOZAK

@dataclass
class ExposureWindow:
    t_lo: float
    t_hi: float
    t_step: float


def exposure_window_prob(
    seq: str,
    exposure_constraint: KozakExposureConstraint = None,
    exposure_params: ExposureWindow = None,
    kozak_pattern: str = KOZAK
) -> tuple[float | None, float | None]:
    """
    Calculate the mean exposure probabilities below and above the target temperature.

    Args:
        seq: RNA sequence.
        exposure_constraint: KozakExposureConstraint object.

    Returns:
        A tuple containing the mean exposure probabilities below and above the target temperature.
    """

    assert exposure_constraint is not None or exposure_params is not None, "Either exposure_constraint or exposure_params must be provided."
    
    temps = None
    if exposure_constraint is not None:
        temps = np.arange(
            exposure_constraint.t_lo,
            exposure_constraint.t_hi + exposure_constraint.t_step,
            exposure_constraint.t_step,
        )
    else:
        temps = np.arange(
            exposure_params.t_lo,
            exposure_params.t_hi + exposure_params.t_step,
            exposure_params.t_step,
        )

    probs_avg = {}
    probs_bn = {}
    for temp in temps:
        try:
            avg, bn = _kozak_exposure_probability_temp(seq, temp, kozak_pattern=kozak_pattern)
            probs_avg[temp] = avg
            probs_bn[temp] = bn
        except ValueError:
            return None, None

    return probs_avg, probs_bn