import nuad.constraints as nc

from utils._loop_size import _loop_size

class LoopConstraint(nc.StrandConstraint):
    """
    Derivative of nuad.constraints.StrandConstraint.
    Objective:
        Keeps the loop size of the sequence within a desired range.
    Note:
        The high and low bounds can be determined by the biological relevance/feasibility of sequences within them.
        (i.e. 1st to 3rd quartiles statistically)
    """

    def __init__(self, weight=1, lo=4, hi=12, unit="nt"):
    
        super().__init__(
            short_description="Loop size constraint",
            description="A constraint that penalizes sequences with loop sizes outside the specified range.",
            weight=weight,
            evaluate=self._evaluate,
        )
        self.lo = lo
        self.hi = hi
        self.unit = unit

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:

        seq = seqs[0].replace("T", "U")
        v = _loop_size(seq)
        e = max(0, self.lo - v) + max(0, v - self.hi)

        return nc.Result(
            excess=e * self.weight,
            value=v,
            unit=self.unit,
            summary=f"Loop size: {v} (excess: {e})"
        )