import nuad.constraints as nc

from utils._stem_length import _stem_length

class StemConstraint(nc.StrandConstraint):
    """
    Derivative of nuad.constraints.StrandConstraint.
    Objective:
        Keeps the stem length of the sequence within a desired range.
    Note:
        The high and low bounds can be determined by the biological relevance/feasibility of sequences within them.
        (i.e. 1st to 3rd quartiles statistically)
    """


    def __init__(self, weight=1, lo=3, hi=25, unit="nt"):
    
        super().__init__(
            short_description="Stem size constraint",
            description="A constraint that penalizes sequences with stem sizes outside the specified range.",
            weight=weight,
            evaluate=self._evaluate,
        )
        self.lo = lo
        self.hi = hi
        self.unit = unit

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:

        seq = seqs[0].replace("T", "U")
        bps, v = _stem_length(seq)
        e = max(0, self.lo - v) + max(0, v - self.hi)

        return nc.Result(
            excess=e * self.weight,
            value=v,
            unit=self.unit,
            summary=f"Stem size: {v} (excess: {e}) with pairs {bps}"
        )