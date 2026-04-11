import nuad.constraints as nc
import nuad.search as ns

from constraints import CompositionConstraint, LoopConstraint, StemConstraint

def main():
    print("Hello from rnat-feedback-design!")

def test():
    # Run NUPACK on template sequence
    pool = nc.DomainPool(name="test_pool", length=30)
    rnat = nc.Domain(name="rnat", pool=pool)
    strand = nc.Strand(name="test_strand", domains=[rnat])
    design = nc.Design(strands=[strand])

    constraints = [
        CompositionConstraint(weight=1, lo=0.4, hi=0.6),
        LoopConstraint(weight=1, lo=3, hi=10),
        StemConstraint(weight=1, lo=5, hi=15),
    ]

    params = ns.SearchParameters(
        constraints=[
            CompositionConstraint(weight=1, lo=0.4, hi=0.6),
            LoopConstraint(weight=1, lo=3, hi=10),
            StemConstraint(weight=1, lo=5, hi=15),
        ],
        # weigh_violations_equally=True,
        out_directory="search_results",
        restart=False,
        # force_overwrite=True,
        report_only_violations=False,
        random_seed=1,
    )

    ns.search_for_sequences(
        design=design, 
        params=params,
    )


if __name__ == "__main__":
    test()
