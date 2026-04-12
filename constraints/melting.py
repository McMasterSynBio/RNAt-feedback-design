import nuad.constraints as nc

from utils._melting import _estimate_tm


class MeltingConstraint(nc.StrandConstraint):
    """Penalises sequences whose predicted Tm deviates from a target temperature.

    Uses NUPACK's partition function (ΔG_ensemble) across a temperature sweep
    to locate the melting temperature, then scores:

        excess = max(0, |Tm - target| - tolerance)
    """

    def __init__(
        self,
        target: float = 37.0,
        tolerance: float = 1.0,
        weight: float = 1.0,
        t_lo: float = 10.0,
        t_hi: float = 90.0,
    ):
        super().__init__(
            short_description="Melting temperature constraint",
            description=(
                f"Penalises sequences whose Tm deviates from "
                f"{target}°C beyond ±{tolerance}°C."
            ),
            weight=weight,
            evaluate=self._evaluate,
        )
        self.target = target
        self.tolerance = tolerance
        self.t_lo = t_lo
        self.t_hi = t_hi

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:
        seq = seqs[0].replace("T", "U")
        tm = _estimate_tm(seq, t_lo=self.t_lo, t_hi=self.t_hi)
        deviation = abs(tm - self.target)
        excess = max(0.0, deviation - self.tolerance)

        return nc.Result(
            excess=excess,
            value=tm,
            unit="°C",
            summary=f"Tm={tm:.1f}°C (target {self.target}°C, excess={excess:.2f})",
        )
