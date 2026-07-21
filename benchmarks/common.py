"""Shared helpers for the DLP.solve() benchmark suite.

Pure Python, driving only the public `smooth` pybind11 surface (`is_prime`,
`factorize`, `sieve_to`, `DLP`) -- no C++ changes. See benchmarks/README.md
for the methodology these helpers implement.
"""

from __future__ import annotations

import csv
import math
import multiprocessing as mp
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import smooth

RESULTS_DIR = Path(__file__).parent / "results"

# Fixed top-level seed so the whole suite is reproducible run to run; every
# sampling call below derives its own random.Random from a seed *string*
# (not this constant directly) so difference benchmarks/trials don't share
# a stream.
BENCH_SEED = 20260705


# ---------------------------------------------------------------------------
# Prime / primitive-root sampling
# ---------------------------------------------------------------------------

def random_prime(bits: int, rng) -> int:
    """A prime with exactly `bits` bits, sampled by rounding a uniform random
    odd bits-bit integer up to the next prime -- the "no special structure"
    baseline every other constructor in this module is contrasted against."""
    n = rng.getrandbits(bits - 1) | (1 << (bits - 1)) | 1
    while not smooth.is_prime(n):
        n += 2
        if n.bit_length() > bits:
            n = rng.getrandbits(bits - 1) | (1 << (bits - 1)) | 1
    return n


def pm1_factors(p: int) -> list[int]:
    """Distinct prime factors of p-1 (smooth.factorize already returns full
    (prime, exponent) factorization; primitive-root testing only needs the
    primes themselves)."""
    return [q for q, _ in smooth.factorize(p - 1)]


def is_primitive_root(cand: int, p: int, pm1_prime_factors: list[int]) -> bool:
    order = p - 1
    return all(pow(cand, order // q, p) != 1 for q in pm1_prime_factors)


def sample_primitive_root(p: int, pm1_prime_factors: list[int], rng, exclude: Optional[int] = None) -> int:
    while True:
        cand = rng.randrange(2, p - 1)
        if cand == exclude:
            continue
        if is_primitive_root(cand, p, pm1_prime_factors):
            return cand


def sample_primitive_pair(p: int, pm1_prime_factors: list[int], rng) -> tuple[int, int]:
    """DLP.solve's documented precondition is that BOTH g and b are genuine
    primitive roots mod p -- this is the only valid input distribution."""
    g = sample_primitive_root(p, pm1_prime_factors, rng)
    b = sample_primitive_root(p, pm1_prime_factors, rng, exclude=g)
    return g, b


def factor_base_size(p: int) -> int:
    """Recomputes the factor-base size k the C++ side would build for this p
    (infra::computeSmoothnessBound: ln B = (1/sqrt2)*sqrt(ln p * ln ln p), then
    k = pi(B)). Purely for benchmark annotation/plotting -- the solver computes
    this itself internally, this is not fed back into it."""
    ln_p = math.log(p)
    ln_B = (1.0 / math.sqrt(2.0)) * math.sqrt(ln_p * math.log(ln_p))
    B = math.ceil(math.exp(ln_B))
    return len(smooth.sieve_to(B))


# ---------------------------------------------------------------------------
# Structural p-1 constructions (Benchmark B: does p-1's shape matter beyond
# what bit length alone predicts?)
# ---------------------------------------------------------------------------

def safe_prime(bits: int, rng) -> int:
    """p = 2q+1 with q prime: the minimal-factor-count extreme for p-1 (just
    {2, q}) -- crtSolve Hensel-lifts over only two factors, one of which is a
    single modulus almost as large as p itself."""
    qbits = bits - 1
    while True:
        q = random_prime(qbits, rng)
        p = 2 * q + 1
        if p.bit_length() == bits and smooth.is_prime(p):
            return p


def smooth_minus_one_prime(
    bits: int,
    rng,
    small_primes: tuple[int, ...] = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31),
) -> int:
    """A prime p whose p-1 factors entirely over `small_primes` (many
    small-modulus Hensel lifts) -- the opposite extreme from safe_prime.
    Construction: build p-1 as a product of random small-prime powers up to
    close to 2**bits, then search small multipliers m until base*m + 1 is
    prime (a standard Pocklington-style construction)."""
    target = 1 << bits
    while True:
        base = 1
        for q in small_primes:
            while base * q < target and rng.random() < 0.7:
                base *= q
        if base.bit_length() < bits - 8:
            continue
        for m in range(1, 20000):
            cand = base * m
            if cand.bit_length() < bits - 1:
                continue
            if cand.bit_length() > bits:
                break
            p = cand + 1
            if p.bit_length() == bits and smooth.is_prime(p):
                return p


# ---------------------------------------------------------------------------
# Timeout-wrapped solve() runner
# ---------------------------------------------------------------------------
#
# DLP.solve is bound without py::call_guard<py::gil_scoped_release>, so it
# holds the GIL for its entire (potentially long) C++ execution -- a Python
# thread-based timeout cannot interrupt it (signal/thread delivery only runs
# between bytecode instructions, never during a blocking C call). A separate
# OS process can always be killed regardless, so that's what bounds a single
# pathological trial's wall time here. The elapsed time reported is measured
# *inside* the child, around only the solve() call, so subprocess spawn
# overhead never pollutes a measurement -- it only eats into the timeout
# budget on the parent side.

def _solve_worker(p: int, g: int, b: int, out_queue) -> None:
    import smooth as _smooth  # re-imported in the fresh child process

    solver = _smooth.DLP(int(p))
    t0 = time.perf_counter()
    x = solver.solve(int(g), int(b))
    t1 = time.perf_counter()
    out_queue.put((t1 - t0, int(x)))


@dataclass
class SolveResult:
    elapsed_s: Optional[float]
    timed_out: bool
    correct: Optional[bool]


def timed_solve(p: int, g: int, b: int, timeout_s: float) -> SolveResult:
    ctx = mp.get_context("spawn")
    out_queue: mp.Queue = ctx.Queue()
    proc = ctx.Process(target=_solve_worker, args=(p, g, b, out_queue))
    proc.start()
    proc.join(timeout_s)

    if proc.is_alive():
        proc.terminate()
        proc.join()
        return SolveResult(elapsed_s=None, timed_out=True, correct=None)

    if out_queue.empty():
        # Abnormal exit (e.g. crash) without ever reaching the queue.put --
        # distinct from a timeout, but still "no measurement".
        return SolveResult(elapsed_s=None, timed_out=False, correct=None)

    elapsed, x = out_queue.get()
    correct = pow(g, x, p) == b
    return SolveResult(elapsed_s=elapsed, timed_out=False, correct=correct)


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------

def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        if not rows:
            return
        # Benchmark A and B rows don't share an identical schema (e.g. "label"
        # only exists on B's rows) -- union the keys (in first-seen order)
        # rather than assuming every row matches rows[0]'s fields.
        fieldnames = list(dict.fromkeys(k for row in rows for k in row.keys()))
        writer = csv.DictWriter(f, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)
