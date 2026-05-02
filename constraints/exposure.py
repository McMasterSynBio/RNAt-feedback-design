import nuad.constraints as nc
import numpy as np

from utils.therm._exposure import _kozak_exposure_probability_temp


class KozakExposureConstraint(nc.StrandConstraint):
    """Penalises sequences whose predicted Kozak exposure deviates from a target probability.

    Uses NUPACK's partition function (ΔG_ensemble) across a temperature sweep
    to locate the Kozak exposure probability.

    Inputs:
        target: Desired Kozak exposure probability (0-100C).
        weight: Relative importance of this constraint in the overall design.
        t_step: Temperature step size for the sweep (°C).
        alpha: Weighting factor between average exposure and bottleneck probability (0-1).
    """

    def __init__(
        self,
        rs_pos: tuple[int,int],
        target: float = 37.0,
        weight: float = 1.0,
        t_step: float = 0.5,
        alpha: float = 0.5
    ):
        
        assert 0 < alpha < 1, "Alpha must be between 0 and 1."

        super().__init__(
            short_description="Kozak exposure constraint",
            description=(
                f"Penalises sequences whose Kozak exposure deviates from "
                f"{target}."
            ),
            weight=weight,
            evaluate=self._evaluate,
        )
        self.target = target
        self.t_lo = target - 10
        self.t_hi = target + 10
        self.t_step = t_step
        self.alpha = alpha

        # right stem position
        self.pos = rs_pos

        # Constraint hardcoded params
        # (for future step constraint update)
        self.theta_low = 0.2
        self.theta_high = 0.8
        self.delta = 2.0
        self.W_max = 4.0
        self.lambda_low = 1.0
        self.lambda_high = 1.0
        self.lambda_center = 2.0
        self.lambda_width = 2.0

    # For later debugging purposes
    @property
    def right_stem(self, seq: str) -> tuple[int,int]:
        return seq[self.pos[0]:self.pos[1]]

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:
        
        seq = seqs[0].replace("T", "U")
        avgs, bns = [], []
        temps = np.arange(self.t_lo, self.t_hi + self.t_step, self.t_step)
        for T in temps:
            avg, bn = _kozak_exposure_probability_temp(seq, T, self.pos)
            avgs.append(avg); bns.append(bn)
        # Consider bottlenecks and averages in error compute
        avgs, bns = np.array(avgs), np.array(bns)
        wp = self.alpha * avgs + (1 - self.alpha) * bns
        # Penalize premature melting (T < Tm) and reward exposure post T > Tm
        below = temps < self.target; above = ~below
        eps = self.t_step # avoid zero div
        penalty_below = np.sum(wp[below] * (self.target - temps[below]))
        penalty_above = np.sum((1.0 - wp[above]) / (temps[above] - self.target + eps))
        # Final value compute
        excess = max(0.0, penalty_below + penalty_above)
        value = np.mean(wp)
        # Return results
        return nc.Result(
            excess=excess,
            value=value,
            unit=None,
            summary=f"Fold error: {excess:.4f} with average weighted exposure probability {value:.4f}",
        )
