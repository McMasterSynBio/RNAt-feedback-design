"""Multi-run design pipeline.

Orchestrates N independent NUAD searches to produce a diverse pool of RNA
thermosensor candidates. Each run starts from a different random seed and
mutates via Hamming-distance steps, scored against the constraint set.
Results are collected, deduplicated, and ranked by exposure-based
optimization score.
"""

# MARK: Imports

from __future__ import annotations

import json, re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import nuad.constraints as nc
import nuad.search as ns

from .helpers.mean_exposure import exposure_window_weight
from utils.therm._exposure import _kozak_exposure_probability_temp

from constraints import (
    CompositionConstraint,
    LoopConstraint,
    StemConstraint,
    KozakExposureConstraint,
)

# MARK: Helpers

@dataclass
class DesignResult:
    """A single designed sequence with its constraint scores."""
    sequence: str
    total_excess: float
    exposure_error: float | None
    exposure_at_target: float | None
    mean_exposure_below: float | None
    mean_exposure_above: float | None
    run_id: int

def _write_settings(out_root: Path, **settings):
    settings_path = out_root / "settings.json"
    with open(settings_path, "w") as fh:
        json.dump(settings, fh, indent=2)


def _parse_report_metrics(report_path: Path) -> tuple[float | None, float | None]:
    total_excess = None
    exposure_error = None

    total_pattern = re.compile(r"total score of constraint violations:\s*([0-9]+(?:\.[0-9]+)?)")
    exposure_pattern = re.compile(r"Fold error:\s*([0-9]+(?:\.[0-9]+)?)")

    with open(report_path) as fh:
        for line in fh:
            if total_excess is None:
                match = total_pattern.search(line)
                if match:
                    total_excess = float(match.group(1))
            if exposure_error is None and "Fold error:" in line:
                match = exposure_pattern.search(line)
                if match:
                    exposure_error = float(match.group(1))

    return total_excess, exposure_error


def _target_exposure(seq: str, target_tm: float, rs_pos: tuple[int,int]) -> float | None:
    try:
        avg, bn = _kozak_exposure_probability_temp(seq, target_tm, rs_pos)
    except ValueError:
        return None
    return 0.5 * (avg + bn)



def _build_constraints(
    target_tm: float,
    rs_pos: tuple[int,int]
) -> list[nc.Constraint]:
    """
    Build a list of constraints for the RNA design.

    Args:
        target_tm: Target melting temperature.

    Returns:
        A list of NUAD constraints.
    """
    return [
        KozakExposureConstraint(rs_pos=rs_pos, target=target_tm, weight=10.0, alpha=0.999),
        # CompositionConstraint(),
        # LoopConstraint(),
        # StemConstraint(),
    ]


# MARK: Single Run
def _run_single(
    run_id: int,
    kozak: str,
    tail: str,
    flank_one_length: int,
    flank_two_length: int,
    constraints: list[nc.Constraint],
    out_root: Path,
    random_seed: int,
    max_iterations: int | None,
) -> DesignResult | None:
    """
    Execute one NUAD search and return the best sequence found.
    Args:
        run_id: Unique identifier for the run.
        kozak: Kozak sequence to include in the design.
        tail: Tail sequence to include in the design.
        flank_one_length: Length of the first RNA sequence flank to design.
        flank_two_length: Length of the second RNA sequence flank to design.
        constraints: List of NUAD constraints to apply.
        out_root: Root directory for output files.
        random_seed: Random seed for reproducibility.
        max_iterations: Maximum number of NUAD iterations.
    Returns:
        A DesignResult object with best sequence metrics.
    """
    pool_flank_one = nc.DomainPool(
        name=f"flank_one_{run_id}", 
        length=flank_one_length
    )

    pool_kozak = nc.DomainPool(
        name=f"kozak_{run_id}", 
        # length=len(kozak),
        possible_sequences=[kozak]
    )

    pool_flank_two = nc.DomainPool(
        name=f"flank_two_{run_id}", 
        length=flank_two_length
    )

    pool_start = nc.DomainPool(
        name=f"start_{run_id}",
        possible_sequences=["AUG"]
    )

    pool_tail = nc.DomainPool(
        name=f"tail_{run_id}",
        possible_sequences=[tail]
    )

    domain_flank_one = nc.Domain(name=f"flank_one_{run_id}", pool=pool_flank_one)
    domain_kozak = nc.Domain(name=f"kozak_{run_id}", pool=pool_kozak)
    domain_flank_two = nc.Domain(name=f"flank_two_{run_id}", pool=pool_flank_two)
    domain_start = nc.Domain(name=f"start_{run_id}", pool=pool_start)
    domain_tail = nc.Domain(name=f"tail_{run_id}", pool=pool_tail)

    strand = nc.Strand(
        name=f"strand_{run_id}", 
        domains=[domain_flank_one, domain_kozak, domain_flank_two, domain_start, domain_tail]
    )
    design = nc.Design(strands=[strand])

    run_dir = str(out_root / f"run_{run_id}")

    params = ns.SearchParameters(
        constraints=constraints,
        out_directory=run_dir,
        restart=False,
        force_overwrite=True,
        save_report_for_all_updates=True,
        random_seed=random_seed,
        max_iterations=max_iterations
    )

    ns.search_for_sequences(
        design=design, 
        params=params
    )

    # Read back the best sequence
    design_file = Path(run_dir) / "design_best.json"
    if not design_file.exists():
        return None

    with open(design_file) as fh:
        data = json.load(fh)

    report_file = Path(run_dir) / "report_best.txt"

    # Extract sequence from design JSON and convert to RNA alphabet
    domains = data.get("domains", [])
    if len(domains) < 5:
        return None
    best_seq = (
        domains[0].get("sequence", "").replace("T", "U") +  # flank_one
        domains[1].get("sequence", "").replace("T", "U") +  # kozak
        domains[2].get("sequence", "").replace("T", "U") +  # flank_two
        domains[3].get("sequence", "").replace("T", "U") +  # start
        domains[4].get("sequence", "").replace("T", "U")    # tail
    )

    total_excess, exposure_error = _parse_report_metrics(report_file) if report_file.exists() else (None, None)
    exposure_at_target = _target_exposure(best_seq, constraints[0].target, constraints[0].pos) if constraints else None
    mean_exposure_below = None
    mean_exposure_above = None
    if constraints and isinstance(constraints[0], KozakExposureConstraint):
        mean_exposure_below, mean_exposure_above = exposure_window_weight(best_seq, constraints[0])

    return DesignResult(
        sequence=best_seq,
        total_excess=float("inf") if total_excess is None else total_excess,
        exposure_error=exposure_error,
        exposure_at_target=exposure_at_target,
        mean_exposure_below=mean_exposure_below,
        mean_exposure_above=mean_exposure_above,
        run_id=run_id,
    )


# MARK: Main Pipeline
def run_design_pipeline(
    kozak: str,
    tail: str,
    flank_one_length: int = 25,
    flank_two_length: int = 5,
    num_runs: int = 10,
    target_tm: float = 37.0,
    max_iterations: int | None = None,
    out_directory: str = "results",
    base_random_seed: int = 42,
) -> list[DesignResult]:
    """Run *num_runs* independent NUAD searches and return ranked results.

    Parameters
    ----------
    kozak : str
        Kozak sequence to include in the design.
    tail : str
        Tail sequence to include in the design.
    flank_one_length : int
        Length of the first RNA sequence flank to design.
    flank_two_length : int
        Length of the second RNA sequence flank to design.
    num_runs : int
        Number of independent search runs (diversity).
    target_tm : float
        Desired melting temperature in °C.
    loop_lo, loop_hi : int
        Allowed loop size range (nt).
    stem_lo, stem_hi : int
        Allowed stem length range (bp).
    max_iterations : int | None
        Cap on NUAD iterations per run.
    out_directory : str
        Root directory for all run outputs.
    base_random_seed : int
        Base seed — each run uses base_random_seed + run_id.
    """
    _write_settings(
        out_root=Path(out_directory),
        kozak=kozak,
        tail=tail,
        flank_one_length=flank_one_length,
        flank_two_length=flank_two_length,
        total_length=flank_one_length + len(kozak) + flank_two_length + 3 + len(tail),
        num_runs=num_runs,
        target_tm=target_tm,
        max_iterations=max_iterations,
        base_random_seed=base_random_seed,
    )

    rs_start = flank_one_length
    rs_end = rs_start + len(kozak) + flank_two_length + 3
    constraints = _build_constraints(target_tm, rs_pos=(rs_start, rs_end))

    out_root = Path(out_directory)
    out_root.mkdir(parents=True, exist_ok=True)

    results: list[DesignResult] = []

    for run_id in range(num_runs):
        print(f"\n{'='*60}")
        print(f"  Run {run_id + 1} / {num_runs}")
        print(f"{'='*60}")

        result = _run_single(
            run_id=run_id,
            kozak=kozak,
            tail=tail,
            flank_one_length=flank_one_length,
            flank_two_length=flank_two_length,
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

    # Rank by the actual exposure-based optimization score.
    unique.sort(
        key=lambda r: (
            r.total_excess,
            float("inf") if r.exposure_error is None else r.exposure_error,
            float("-inf") if r.mean_exposure_above is None else -r.mean_exposure_above,
            float("inf") if r.mean_exposure_below is None else r.mean_exposure_below,
        )
    )

    # Write summary
    summary_path = out_root / "summary.csv"
    with open(summary_path, "w") as fh:
        fh.write(
            "rank,sequence,kozak,total_excess,exposure_error,exposure_at_target,"
            "mean_exposure_below,mean_exposure_above,run_id\n"
        )
        for rank, r in enumerate(unique, 1):
            exposure_error = "" if r.exposure_error is None else f"{r.exposure_error:.4f}"
            exposure_at_target = "" if r.exposure_at_target is None else f"{r.exposure_at_target:.4f}"
            mean_exposure_below = "" if r.mean_exposure_below is None else f"{r.mean_exposure_below:.4f}"
            mean_exposure_above = "" if r.mean_exposure_above is None else f"{r.mean_exposure_above:.4f}"
            fh.write(
                f"{rank},{r.sequence},{r.sequence[rs_start:rs_end]},{r.total_excess:.2f},{exposure_error},{exposure_at_target},"
                f"{mean_exposure_below},{mean_exposure_above},{r.run_id}\n"
            )

    print(f"\n{'='*60}")
    print(f"  Pipeline complete — {len(unique)} unique sequences")
    print(f"  Summary written to {summary_path}")
    print(f"{'='*60}\n")

    for rank, r in enumerate(unique[:10], 1):
        exposure_error = "NA" if r.exposure_error is None else f"{r.exposure_error:.4f}"
        below_str = "NA" if r.mean_exposure_below is None else f"{r.mean_exposure_below:.4f}"
        above_str = "NA" if r.mean_exposure_above is None else f"{r.mean_exposure_above:.4f}"
        print(
            f"  #{rank}  score={r.total_excess:.2f}  exposure_error={exposure_error}  "
            f"below={below_str}  above={above_str} right_stem={r.sequence[rs_start:rs_end]}  {r.sequence}"
        )

    return unique
