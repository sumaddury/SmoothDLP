"""Benchmark B -- "worst case" as two targeted structural stress tests.

Bit length alone (Benchmark A) already explains most of the variance in
solve() time, per infra::ProblemParams's sub-exponential smoothness-bound
heuristic. Rather than guess at a vague "worst case", this benchmark isolates
two *specific, code-grounded* structural axes documented in CLAUDE.md/
infra.h that could make two same-bit-length primes behave differently:

  (1) p-1's factorization shape. crtSolve Hensel-lifts the linear solve up
      through every prime-power factor of p-1 (infra.h / lin_alg.h). We hold
      bit length fixed (60 bits, cheap enough for several trials) and
      compare two structural extremes against Benchmark A's unstructured
      60-bit trials as a control:
        - safe primes (p = 2q+1): the minimal-factor-count extreme -- just
          {2, q}, one of which is a single huge-modulus Hensel lift.
        - smooth-(p-1) primes: p-1 factors entirely over small primes --
          many small-modulus Hensel lifts instead of one big one.

  (2) The Wiedemann/SparseElimination backend switch. lin_alg.cpp dispatches
      on factor-base size k via a hard-coded WIEDEMANN_THRESHOLD = 4000.
      Separately computing k = pi(B) from infra::computeSmoothnessBound shows
      k crosses 4000 right around 80 bits (k(76) ~ 2990, k(84) ~ 5410) -- so
      a prime landing just above the threshold could be disproportionately
      slower than the smooth bit-length trend from Benchmark A predicts, for
      reasons that are an implementation artifact (backend switch) rather
      than genuine problem difficulty. We run one trial each at 76 and 84
      bits to check for a discontinuity. This pair is deliberately small
      (not a full sweep) since each trial here is expensive and it's meant
      to be read alongside Benchmark A's 70/80-bit points, not standalone.
"""

from __future__ import annotations

import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmarks import common as bc

PM1_SHAPE_BITS = 60
PM1_SHAPE_TRIALS = 5
PM1_SHAPE_TIMEOUT_S = 150.0

THRESHOLD_BITS = (76, 84)
THRESHOLD_TIMEOUT_S = 400.0


def _run_trial(label: str, bits: int, p: int, rng, timeout_s: float, trial_idx: int, t_start: float, verbose: bool) -> dict:
    factors = bc.pm1_factors(p)
    k = bc.factor_base_size(p)
    g, b = bc.sample_primitive_pair(p, factors, rng)
    result = bc.timed_solve(p, g, b, timeout_s)

    if verbose:
        status = (
            "TIMEOUT" if result.timed_out
            else f"{result.elapsed_s:.3f}s" if result.elapsed_s is not None
            else "ERROR"
        )
        print(
            f"[worst:{label}] bits={bits:>3} trial={trial_idx} k={k:>5} "
            f"num_pm1_factors={len(factors)} p={p} -> {status} "
            f"(t+{time.time() - t_start:6.1f}s elapsed)",
            flush=True,
        )

    return {
        "bench": "worst",
        "label": label,
        "bits": bits,
        "trial_idx": trial_idx,
        "p": p,
        "g": g,
        "b": b,
        "k": k,
        "num_pm1_factors": len(factors),
        "elapsed_s": result.elapsed_s,
        "timed_out": result.timed_out,
        "correct": result.correct,
    }


def run_pm1_shape(verbose: bool = True) -> list[dict]:
    rows: list[dict] = []
    t_start = time.time()

    safe_rng = random.Random(f"{bc.BENCH_SEED}-worst-safe")
    pair_rng = random.Random(f"{bc.BENCH_SEED}-worst-safe-pairs")
    for i in range(PM1_SHAPE_TRIALS):
        p = bc.safe_prime(PM1_SHAPE_BITS, safe_rng)
        rows.append(_run_trial("safe_prime", PM1_SHAPE_BITS, p, pair_rng, PM1_SHAPE_TIMEOUT_S, i, t_start, verbose))

    smooth_rng = random.Random(f"{bc.BENCH_SEED}-worst-smoothpm1")
    pair_rng2 = random.Random(f"{bc.BENCH_SEED}-worst-smoothpm1-pairs")
    for i in range(PM1_SHAPE_TRIALS):
        p = bc.smooth_minus_one_prime(PM1_SHAPE_BITS, smooth_rng)
        rows.append(_run_trial("smooth_pm1", PM1_SHAPE_BITS, p, pair_rng2, PM1_SHAPE_TIMEOUT_S, i, t_start, verbose))

    return rows


def run_threshold_straddle(verbose: bool = True) -> list[dict]:
    rows: list[dict] = []
    t_start = time.time()

    for bits in THRESHOLD_BITS:
        prime_rng = random.Random(f"{bc.BENCH_SEED}-worst-threshold-{bits}")
        pair_rng = random.Random(f"{bc.BENCH_SEED}-worst-threshold-pairs-{bits}")
        p = bc.random_prime(bits, prime_rng)
        label = f"threshold_{bits}b"
        rows.append(_run_trial(label, bits, p, pair_rng, THRESHOLD_TIMEOUT_S, 0, t_start, verbose))

    return rows


def run(verbose: bool = True) -> list[dict]:
    return run_pm1_shape(verbose=verbose) + run_threshold_straddle(verbose=verbose)


if __name__ == "__main__":
    rows = run()
    bc.write_csv(rows, bc.RESULTS_DIR / "worst_case.csv")
    print(f"\nwrote {len(rows)} rows to {bc.RESULTS_DIR / 'worst_case.csv'}")
