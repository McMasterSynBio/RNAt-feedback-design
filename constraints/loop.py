import nuad.constraints as nc

from ..utils._loop_size import _loop_size

class LoopConstraint(nc.StrandConstraint):

    def __init__(self, weight=1, lo=0, hi=100, unit="nt"):
    
        super().__init__(
            short_description="Loop size constraint",
            description="A constraint that penalizes sequences with loop sizes outside the specified range.",
            weight=weight,
            evaluate=self._evaluate,
        )

        self.w = weight
        self.lo = lo
        self.hi = hi
        self.unit = unit

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:

        seq = seqs[0]
        v = _loop_size(seq)
        e = max(0, self.lo - v) + max(0, v - self.hi)

        return nc.Result(
            excess=e,
            value=v,
            unit=self.unit,
            summary=f"Loop size: {v} (excess: {e})"
        )