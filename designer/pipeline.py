"""Multi-run design pipeline.

Orchestrates N independent NUAD searches to produce a diverse pool of RNA
thermosensor candidates.  Each run starts from a different random seed and
mutates via Hamming-distance steps, scored against the constraint set.
Results are collected, deduplicated, and ranked by proximity to target Tm.
"""

# MARK: Imports

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import nuad.constraints as nc
import nuad.search as ns
from utils.therm._melting import estimate_tm

from constraints import (
    CompositionConstraint,
    LoopConstraint,
    MeltingConstraint,
    StemConstraint,
    KozakExposureConstraint,
)

# MARK: Helpers

@dataclass
class DesignResult:
    """A single designed sequence with its constraint scores."""
    sequence: str
    tm: float | None
    total_excess: float
    run_id: int

def _write_settings(out_root: Path, **settings):
    settings_path = out_root / "settings.json"
    with open(settings_path, "w") as fh:
        json.dump(settings, fh, indent=2)



def _build_constraints(
    target_tm: float,
) -> list[nc.Constraint]:
    return [
        KozakExposureConstraint(target=target_tm, weight=6.0),
        CompositionConstraint(),
        LoopConstraint(),
        StemConstraint(),
    ]


# MARK: Single Run
def _run_single(
    run_id: int,
    seq_length: int,
    constraints: list[nc.Constraint],
    out_root: Path,
    random_seed: int,
    max_iterations: int | None,
) -> DesignResult | None:
    """Execute one NUAD search and return the best sequence found."""
    pool = nc.DomainPool(name=f"pool_{run_id}", length=seq_length)
    domain = nc.Domain(name=f"rnat_{run_id}", pool=pool)
    strand = nc.Strand(name=f"strand_{run_id}", domains=[domain])
    design = nc.Design(strands=[strand])

    run_dir = str(out_root / f"run_{run_id}")

    params = ns.SearchParameters(
        constraints=constraints,
        out_directory=run_dir,
        restart=False,
        force_overwrite=True,
        report_only_violations=False,
        random_seed=random_seed,
        max_iterations=max_iterations,
    )

    ns.search_for_sequences(design=design, params=params)

    # Read back the best sequence
    design_file = Path(run_dir) / "design_best.json"
    if not design_file.exists():
        return None

    with open(design_file) as fh:
        data = json.load(fh)

    # Extract sequence from design JSON and convert to RNA alphabet
    domains = data.get("domains", [])
    if not domains:
        return None
    best_seq = domains[0].get("sequence", "").replace("T", "U")

    tm = estimate_tm(best_seq)

    return DesignResult(
        sequence=best_seq,
        tm=tm,
        total_excess=0.0,  # will be filled by ranking
        run_id=run_id,
    )


# MARK: Main Pipeline
def run_design_pipeline(
    seq_length: int = 30,
    num_runs: int = 10,
    target_tm: float = 37.0,
    max_iterations: int | None = None,
    out_directory: str = "results",
    base_random_seed: int = 42,
) -> list[DesignResult]:
    """Run *num_runs* independent NUAD searches and return ranked results.

    Parameters
    ----------
    seq_length : int
        Length of the RNA sequence to design.
    num_runs : int
        Number of independent search runs (diversity).
    target_tm : float
        Desired melting temperature in °C.
    loop_lo, loop_hi : int
        Allowed loop size range (nt).
    stem_lo, stem_hi : int
        Allowed stem length range (bp).
    max_iterations : int | None
        Cap on NUAD iterations per run (None = until convergence).
    out_directory : str
        Root directory for all run outputs.
    base_random_seed : int
        Base seed — each run uses base_random_seed + run_id.
    """
    _write_settings(
        out_root=Path(out_directory),
        seq_length=seq_length,
        num_runs=num_runs,
        target_tm=target_tm,
        max_iterations=max_iterations,
        base_random_seed=base_random_seed,
    )

    constraints = _build_constraints(target_tm)

    out_root = Path(out_directory)
    out_root.mkdir(parents=True, exist_ok=True)

    results: list[DesignResult] = []

    for run_id in range(num_runs):
        print(f"\n{'='*60}")
        print(f"  Run {run_id + 1} / {num_runs}")
        print(f"{'='*60}")

        result = _run_single(
            run_id=run_id,
            seq_length=seq_length,
            constraints=constraints,
            out_root=out_root,
            random_seed=base_random_seed + run_id,
            max_iterations=max_iterations,
        )
        if result is not None:
            results.append(result)

    # Deduplicate by sequence
    seen = set()
    unique: list[DesignResult] = []
    for r in results:
        if r.sequence not in seen:
            seen.add(r.sequence)
            unique.append(r)

    # Sort by |Tm − target| (closest to target first)
    unique.sort(key=lambda r: abs((r.tm or 999) - target_tm))

    # Write summary
    summary_path = out_root / "summary.csv"
    with open(summary_path, "w") as fh:
        fh.write("rank,sequence,Tm_C,delta_Tm,run_id\n")
        for rank, r in enumerate(unique, 1):
            delta = abs((r.tm or 999) - target_tm)
            fh.write(f"{rank},{r.sequence},{r.tm:.1f},{delta:.1f},{r.run_id}\n")

    print(f"\n{'='*60}")
    print(f"  Pipeline complete — {len(unique)} unique sequences")
    print(f"  Summary written to {summary_path}")
    print(f"{'='*60}\n")

    for rank, r in enumerate(unique[:10], 1):
        delta = abs((r.tm or 999) - target_tm)
        print(f"  #{rank}  Tm={r.tm:.1f}°C  (Δ={delta:.1f})  {r.sequence}")

    return unique
