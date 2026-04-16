import nuad.constraints as nc

from utils._composition import _au_gc_comp

class CompositionConstraint(nc.StrandConstraint):
    """
    A derivative of nuad.constraints.StrandConstraint. 
    Objective:
        Compute the GC composition content and keeps it within the range passed on to it.
    Note:
        The ideal range of choice can be determined base on the use case and biological feasilibity/relevance.
    """

    def __init__(self, weight=1, lo=0.2, hi=0.8, unit="nt"):
    
        super().__init__(
            short_description="Composition constraint",
            description="A constraint that penalizes sequences with AU/GC composition outside the specified range.",
            weight=weight,
            evaluate=self._evaluate,
        )
        self.lo = lo
        self.hi = hi
        self.unit = unit

    def _evaluate(self, seqs: tuple[str, ...], strand: nc.Strand | None) -> nc.Result:

        seq = seqs[0].replace("T", "U")
        v = _au_gc_comp(seq)
        e = 0 if self.lo <= v <= self.hi else 1e5

        return nc.Result(
            excess=e * self.weight,
            value=v,
            unit=self.unit,
            summary=f"AU/GC composition: {v} (excess: {e})"
        )