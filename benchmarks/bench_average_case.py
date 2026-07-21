"""Benchmark A -- average case: DLP.solve() wall time across a bit-length
sweep of *unstructured* random primes.

Math / sampling, see benchmarks/README.md for the full writeup:

  - infra::ProblemParams picks the smoothness bound via the classic
    L-notation heuristic ln B = (1/sqrt2) * sqrt(ln p * ln ln p), which makes
    the factor-base size k = pi(B) -- and hence the expected relation-
    collection + linear-solve cost -- a sub-exponential function of p's bit
    length. Bit length is therefore the dominant, near-deterministic driver
    of solve() cost; this benchmark sweeps it directly.
  - For each bit length we sample P independent random primes (rounding a
    uniform random bits-bit odd integer up to the next prime via
    smooth.is_prime -- no special structure), and for each prime, R
    independent primitive-root pairs (g, b). This P x R grid separates
    *between-prime* variance (does this p's own structure, e.g. its p-1
    shape, make it slow?) from *within-prime* variance (is it just luck in
    which relations get hit?) -- see benchmarks/README.md's variance
    decomposition discussion.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmarks import common as bc

# bits -> (P primes, R pairs per prime, per-trial timeout seconds).
# Timeouts are generous safety valves (they should essentially never fire
# under normal conditions) sized well above the typical time measured during
# design-time calibration, not a target time.
SWEEP: dict[int, tuple[int, int, float]] = {
    30: (4, 2, 60.0),
    40: (4, 2, 60.0),
    50: (3, 2, 90.0),
    60: (3, 2, 150.0),
    70: (2, 2, 250.0),
    80: (2, 1, 500.0),
}


def run(sweep: dict[int, tuple[int, int, float]] = SWEEP, verbose: bool = True) -> list[dict]:
    rows: list[dict] = []
    t_start = time.time()

    for bits, (num_primes, num_pairs, timeout_s) in sweep.items():
        prime_rng = __import__("random").Random(f"{bc.BENCH_SEED}-average-primes-{bits}")

        for prime_idx in range(num_primes):
            p = bc.random_prime(bits, prime_rng)
            factors = bc.pm1_factors(p)
            k = bc.factor_base_size(p)
            pair_rng = __import__("random").Random(f"{bc.BENCH_SEED}-average-pairs-{bits}-{prime_idx}")

            for pair_idx in range(num_pairs):
                g, b = bc.sample_primitive_pair(p, factors, pair_rng)
                result = bc.timed_solve(p, g, b, timeout_s)

                if verbose:
                    status = (
                        "TIMEOUT" if result.timed_out
                        else f"{result.elapsed_s:.3f}s" if result.elapsed_s is not None
                        else "ERROR"
                    )
                    print(
                        f"[average] bits={bits:>3} prime={prime_idx} pair={pair_idx} "
                        f"k={k:>5} p={p} -> {status} "
                        f"(t+{time.time() - t_start:6.1f}s elapsed)",
                        flush=True,
                    )

                rows.append({
                    "bench": "average",
                    "bits": bits,
                    "prime_idx": prime_idx,
                    "pair_idx": pair_idx,
                    "p": p,
                    "g": g,
                    "b": b,
                    "k": k,
                    "num_pm1_factors": len(factors),
                    "elapsed_s": result.elapsed_s,
                    "timed_out": result.timed_out,
                    "correct": result.correct,
                })

    return rows


if __name__ == "__main__":
    rows = run()
    bc.write_csv(rows, bc.RESULTS_DIR / "average_case.csv")
    print(f"\nwrote {len(rows)} rows to {bc.RESULTS_DIR / 'average_case.csv'}")
