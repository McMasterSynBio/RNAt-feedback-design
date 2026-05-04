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
        self.t_lo = target - 20
        self.t_hi = target + 20
        self.t_step = t_step
        self.alpha = alpha

        # right stem position
        self.pos = rs_pos

        # Constraint hardcoded params
        # (for future step constraint update)
        self.theta_low = 0.1
        self.theta_high = 0.9
        self.delta = 1.5
        self.W_max = 2.5
        self.lambda_low = 1.5
        self.lambda_high = 1.5
        self.lambda_center = 2.0
        self.lambda_width = 3.0

    # For later debugging purposes
    @property
    def right_stem(self, seq: str) -> tuple[int,int]:
        return seq[self.pos[0]:self.pos[1]]

    def _weighted_exposure_curve(self, seq: str) -> tuple[np.ndarray, np.ndarray]:
        temps = np.arange(self.t_lo, self.t_hi + self.t_step, self.t_step)
        avgs, bns = [], []
        for temp in temps:
            avg, bn = _kozak_exposure_probability_temp(seq, temp, self.pos)
            avgs.append(avg)
            bns.append(bn)

        avgs = np.array(avgs)
        bns = np.array(bns)
        weighted_probs = self.alpha * avgs + (1 - self.alpha) * bns
        return temps, weighted_probs

    def _interpolate_crossing(
        self,
        temps: np.ndarray,
        values: np.ndarray,
        threshold: float,
    ) -> float | None:
        exact_hits = np.where(np.isclose(values, threshold))[0]
        if exact_hits.size:
            return float(temps[exact_hits[0]])

        crossings: list[float] = []
        for index in range(len(values) - 1):
            left = values[index]
            right = values[index + 1]
            if left < threshold <= right:
                fraction = (threshold - left) / (right - left)
                crossing = temps[index] + fraction * (temps[index + 1] - temps[index])
                crossings.append(float(crossing))

        if not crossings:
            return None

        return min(crossings, key=lambda crossing: abs(crossing - self.target))

    def _evaluate_step(self, seq: str) -> tuple[float, float, str]:
        temps, weighted_probs = self._weighted_exposure_curve(seq)

        below = temps <= (self.target - self.delta)
        above = temps >= (self.target + self.delta)

        penalty_low = float(
            np.sum(np.maximum(0.0, weighted_probs[below] - self.theta_low) ** 2)
        ) if np.any(below) else 0.0
        penalty_high = float(
            np.sum(np.maximum(0.0, self.theta_high - weighted_probs[above]) ** 2)
        ) if np.any(above) else 0.0

        t10 = self._interpolate_crossing(temps, weighted_probs, 0.1)
        t50 = self._interpolate_crossing(temps, weighted_probs, 0.5)
        t90 = self._interpolate_crossing(temps, weighted_probs, 0.9)

        missing_center = (self.t_hi - self.t_lo) ** 2
        penalty_center = missing_center if t50 is None else (t50 - self.target) ** 2

        if t10 is None or t90 is None:
            penalty_width = self.W_max ** 2
            width = None
        else:
            width = max(0.0, t90 - t10)
            penalty_width = max(0.0, width - self.W_max) ** 2

        excess = (
            self.lambda_low * penalty_low
            + self.lambda_high * penalty_high
            + self.lambda_center * penalty_center
            + self.lambda_width * penalty_width
        )
        value = float(np.mean(weighted_probs))

        summary = (
            f"Step error: {excess:.4f}; low={penalty_low:.4f}; high={penalty_high:.4f}; "
            f"center={penalty_center:.4f}; width={penalty_width:.4f}; "
            f"T10={('NA' if t10 is None else f'{t10:.2f}')}; "
            f"T50={('NA' if t50 is None else f'{t50:.2f}')}; "
            f"T90={('NA' if t90 is None else f'{t90:.2f}')}; "
            f"mean_exposure={value:.4f}"
        )
        return excess, value, summary

    def _evaluate_legacy(self, seq: str) -> tuple[float, float, str]:
        """Deprecated broad-window exposure score kept for comparison/debugging."""
        temps, weighted_probs = self._weighted_exposure_curve(seq)
        below = temps < self.target
        above = ~below
        eps = self.t_step
        penalty_below = np.sum(weighted_probs[below] * (self.target - temps[below]))
        penalty_above = np.sum((1.0 - weighted_probs[above]) / (temps[above] - self.target + eps))
        excess = max(0.0, float(penalty_below + penalty_above))
        value = float(np.mean(weighted_probs))
        summary = (
            f"Deprecated fold error: {excess:.4f} with average weighted exposure probability {value:.4f}"
        )
        return excess, value, summary

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:
        seq = seqs[0].replace("T", "U")
        excess, value, summary = self._evaluate_step(seq)
        return nc.Result(
            excess=excess,
            value=value,
            unit=None,
            summary=summary,
        )
