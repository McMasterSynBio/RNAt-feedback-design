import nuad.constraints as nc

from utils._composition import _au_gc_comp

class CompositionConstraint(nc.StrandConstraint):

    def __init__(self, weight=1, lo=0, hi=100, unit="nt"):
    
        super().__init__(
            short_description="Composition constraint",
            description="A constraint that penalizes sequences with AU/GC composition outside the specified range.",
            weight=weight,
            evaluate=self._evaluate,
        )

        self.w = weight
        self.lo = lo
        self.hi = hi
        self.unit = unit

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:

        seq = seqs[0]
        v = _au_gc_comp(seq)
        e = max(0, self.lo - v) + max(0, v - self.hi)

        return nc.Result(
            excess=e,
            value=v,
            unit=self.unit,
            summary=f"AU/GC composition: {v} (excess: {e})"
        )