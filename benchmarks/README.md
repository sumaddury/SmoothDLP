# `DLP.solve()` benchmark suite

Everything here is additive Python driving only the public `smooth` pybind11
surface (`is_prime`, `factorize`, `sieve_to`, `DLP`) — no changes to `src/`.
The suite benchmarks exactly one interface, `DLP.solve(g, b)`; `DLP(p)`
construction (factor base, product tree, `factorize(p-1)`) is treated as
setup and excluded from every timed measurement, since the class is designed
to be built once per `p` and `solve()`'d many times.

The previous benchmarking artifacts under `assets/` (`*_complexity.png`,
`large_prime_dataset.json`, etc.) and anything under `scripts/` predate this
suite, are not rigorous, and benchmark individual primitives
(`is_prime`, `factorize`, ...) rather than the end-to-end solver. They are
left untouched; this suite does not build on or replace them.

## Running it

```bash
conda activate smoothnumbers-dev
pip install -r benchmarks/requirements.txt   # matplotlib, for plots only
python -m benchmarks.run_all
```

Typical wall time is ~15-20 minutes; a hard per-trial timeout (see below)
caps the worst case around 30-35 minutes. Individual pieces can also be run
standalone: `python -m benchmarks.bench_average_case`,
`python -m benchmarks.bench_worst_case`.

Outputs land in `benchmarks/results/`: `average_case.csv`, `worst_case.csv`,
`all_results.csv` (raw per-trial rows), plus `average_case_scaling.png` and
`pm1_shape_comparison.png`.

## Why bit length is the dominant axis

`infra::ProblemParams` picks its smoothness bound via the classic L-notation
heuristic:

```
ln B = (1/sqrt(2)) * sqrt(ln p * ln ln p)
```

which fixes the factor-base size `k = pi(B)` as a sub-exponential function of
`p`'s bit length alone — independent of *which* `p` at that bit length. This
is why bit length, not some other input property, is expected to explain
most of the variance in `solve()` time, and why the suite is built around a
bit-length sweep rather than trying to hand-pick "hard" individual primes.

Measured factor-base sizes at the sweep's bit lengths (computed with the
same formula, replicated in `common.factor_base_size`, purely for
annotation — the solver computes this itself internally):

| bits | 30 | 40 | 50 | 60 | 70 | 76 | 80 | 84 |
|---|---|---|---|---|---|---|---|---|
| k  | 58 | 153 | 370 | 859 | 1900 | ~2990 | 4035 | ~5410 |

`lin_alg.cpp` switches its sparse linear solve from `SparseElimination` to
Wiedemann at a hard-coded `WIEDEMANN_THRESHOLD = 4000` — which `k` crosses
almost exactly at 80 bits. That crossover motivates Benchmark B's second
sub-experiment below.

## Sampling: primes, primitive roots

`DLP.solve`'s only documented precondition is that both `g` and `b` are
genuine primitive roots mod `p` (order exactly `p-1`) — this is the *only*
valid input distribution, so every benchmark samples from it:

- **Primes**: draw a uniform random odd `bits`-bit integer, round up to the
  next prime via `is_prime` (`common.random_prime`). No structure is
  imposed — this is the "no special algebraic properties" baseline that
  Benchmark B's constructions are deliberately contrasted against.
- **Primitive roots**: factor `p-1` via `factorize` (needed regardless, to
  even test primitivity), then rejection-sample candidates in `[2, p-2]`
  against the standard test `g^((p-1)/q) != 1` for every prime `q | p-1`
  (`common.sample_primitive_root`).

## Benchmark A — average case (bit-length sweep)

`benchmarks/bench_average_case.py`. For `bits in {30, 40, 50, 60, 70, 80}`,
sample `P` independent random primes, and for each prime, `R` independent
`(g, b)` primitive-root pairs, running all `P*R` trials through `DLP.solve`.

This `P x R` grid is the one deliberate improvement on "just sample random
pairs": it separates two different sources of variance that a flat sample
would conflate —

- **between-prime variance**: does *this specific* `p` (e.g. its `p-1`
  shape) make it consistently slower or faster than another `p` at the same
  bit length?
- **within-prime variance**: given a fixed `p`, how much does solve time
  still vary just from which relations happen to get hit first (the
  "luck" the user asked about)?

`analyze.variance_decomposition` reports both per bit length: the stdev of
each prime's mean solve time (between) vs. the mean stdev across a single
prime's `R` pairs (within). Large between-prime relative to within-prime
means a prime's own structure is the deciding factor; the reverse means
solve time is mostly noise around a shared bit-length trend.

`(P, R, timeout)` per bit length — timeouts are generous safety valves sized
well above the calibration measurements below, not targets:

| bits | P | R | timeout |
|---|---|---|---|
| 30 | 4 | 2 | 60s |
| 40 | 4 | 2 | 60s |
| 50 | 3 | 2 | 90s |
| 60 | 3 | 2 | 150s |
| 70 | 2 | 2 | 250s |
| 80 | 2 | 1 | 500s |

Design-time calibration (single trials, informal, not part of the recorded
suite): 25b/0.007s, 30b/0.015s, 35b/0.042s, 40b/0.10s, 45b/0.25s, 50b/1.5s,
55b/2.8s, 60b/3.0s, 65b/10.3s, 70b/21.8s — growth is visibly super-linear in
bits (consistent with the sub-exponential `k(bits)` above), but even the
70-bit point sits roughly 250x faster than the paper's own reported ~5000s
at 70 bits.

## Benchmark B — worst case, as two targeted structural stress tests

`benchmarks/bench_worst_case.py`. The user's own framing — "worst case is
tougher [to define]... runtime will strongly correlate with bit length" —
is taken literally: rather than invent a vague adversarial input, this
isolates exactly two structural axes with a direct, documented link to a
specific code path, both held at a fixed bit length so bit length itself
isn't a confound.

**B1 — `p-1` factorization shape, at 60 bits (5 trials each).**
`infra::crtSolve` / `henselLift` lift the linear solve up through every
prime-power factor of `p-1`. Two structural extremes:

- `safe_prime(bits)`: `p = 2q+1`, `q` prime — the minimal-factor-count case,
  just `{2, q}`, one of which is a single Hensel lift over a modulus almost
  as large as `p` itself.
- `smooth_minus_one_prime(bits)`: `p-1` built entirely from small-prime
  powers (Pocklington-style construction: multiply random small-prime powers
  up to near `2**bits`, then search small multipliers until `+1` is prime) —
  many small-modulus Hensel lifts instead of one huge one.

Both are compared against Benchmark A's unstructured 60-bit trials as a
control group.

**B2 — Wiedemann/SparseElimination threshold straddle, one trial each at 76
and 84 bits.** `k(76) ~ 2990` sits below `WIEDEMANN_THRESHOLD = 4000`,
`k(84) ~ 5410` sits above it. A single trial each (not a full sweep — this
is meant to be read as two extra points layered onto Benchmark A's plot,
flagging whether the backend switch produces a visible discontinuity beyond
what the smooth bit-length trend predicts) at a generous 400s timeout.

## Timeout handling

`DLP.solve` is bound without `py::call_guard<py::gil_scoped_release>`, so it
holds the GIL for its entire C++ execution — a Python thread-based timeout
cannot interrupt a blocking call like this (signal/thread delivery only
happens between bytecode instructions). `common.timed_solve` therefore runs
each trial in a fresh subprocess (`multiprocessing`, `spawn` context) and
enforces the cap by `terminate()`-ing the process if it's still alive after
the timeout — the only way to reliably bound a single pathological trial.
The elapsed time recorded is measured *inside* the child, around only the
`solve()` call, so subprocess spawn overhead never pollutes a measurement;
it only eats into the wall-clock budget. A trial that hits its cap is
recorded as right-censored (`timed_out=True`, `elapsed_s=None`) rather than
silently dropped or left to hang.

## Findings

Full run: 45 trials, ~22 minutes wall time (within the ~15-22 min typical /
~30-35 min ceiling estimate). Raw data in `results/all_results.csv` (and the
per-benchmark CSVs); plots in `results/average_case_scaling.png` and
`results/pm1_shape_comparison.png`.

**Benchmark A — average case.**

| bits | n | mean (s) | median (s) | stdev (s) | min (s) | max (s) |
|---|---|---|---|---|---|---|
| 30 | 8 | 0.013 | 0.012 | 0.002 | 0.011 | 0.016 |
| 40 | 8 | 0.189 | 0.141 | 0.099 | 0.112 | 0.368 |
| 50 | 6 | 0.826 | 0.751 | 0.381 | 0.491 | 1.511 |
| 60 | 6 | 10.821 | 4.266 | 10.764 | 3.570 | 28.215 |
| 70 | 4 | 34.310 | 25.568 | 18.867 | 23.616 | 62.489 |
| 80 | 2 | 306.801 | 306.801 | 125.056 | 218.373 | 395.229 |

- Growth is visibly super-linear in bits on a log-time axis, consistent with
  the sub-exponential `k(bits)` predictor. At 70 bits, the mean (34.3s) is
  **~146x faster** than the paper's own reported ~5000s at 70 bits, and even
  the slowest of the 4 trials (62.5s) is ~80x faster — this reproduces, with
  actual measurement rather than a single anecdotal run, the improvement the
  user had informally observed.
- Variance decomposition (between-prime vs. within-prime stdev, seconds):

  | bits | #primes | between-prime | within-prime |
  |---|---|---|---|
  | 30 | 4 | 0.001 | 0.001 |
  | 40 | 4 | 0.069 | 0.077 |
  | 50 | 3 | 0.309 | 0.228 |
  | 60 | 3 | 6.113 | 9.666 |
  | 70 | 2 | 15.025 | 12.463 |

  Between- and within-prime variance are comparable in magnitude at every
  bit length with enough samples to compare (30-70 bits) — there's no strong
  evidence that "which prime you happened to draw" dominates over "which
  relations happened to get hit" for generically-sampled primes. Both
  matter, roughly equally, which is itself the honest answer to "is runtime
  mostly bit-length, prime-specific structure, or luck": it's bit length
  first, then a roughly even split of the other two. (80 bits used `R=1`
  pair/prime for time-budget reasons, so within-prime variance isn't
  measurable there.)
- The 60-bit group's spread is dominated by two outliers (20.5s, 28.2s vs.
  ~3.6-4.3s for the other four) — a first direct, reproducible illustration
  of the "deep algebraic properties not captured by bit length alone" the
  user flagged as a concern; Benchmark B investigates one concrete candidate
  explanation (p-1 shape) for exactly this kind of outlier.

**Benchmark B1 — p-1 shape at 60 bits.**

| label | n | mean (s) | median (s) | stdev (s) |
|---|---|---|---|---|
| safe_prime (p-1=2q) | 5 | 3.102 | 3.077 | 0.060 |
| smooth_pm1 (many small factors) | 5 | 7.834 | 7.828 | 0.536 |
| generic (Benchmark A control) | 6 | 10.821 | 4.266 | 10.764 |

Below the Wiedemann threshold, p-1's shape has a clear, low-variance effect:
safe primes solve consistently in ~3.1s (tight, stdev 0.06s), smooth-p-1
primes consistently ~2.5x slower at ~7.8s, and both structured groups are
individually *more consistent* (far lower stdev) than the generic control —
suggesting the generic group's own high variance (10.8 +/- 10.8s) is itself
substantially explained by which unstructured p-1 shape each of those
random primes happened to land on.

**Benchmark B2 — Wiedemann-threshold straddle (76 vs. 84 bits).**

| bits | k | num p-1 factors | result |
|---|---|---|---|
| 76 | 2976 | 5 | 76.9s (completed) |
| 84 | 5400 | 2 | **did not complete within the 400s cap** |

The 76-bit trial (below the `k=4000` threshold) actually finished *faster*
than the 70-80 bit trend would suggest by interpolation (~127s expected vs.
77s actual) — no evidence of a slowdown approaching the threshold from
below. The 84-bit trial (above the threshold) is the more interesting data
point: it didn't finish in 400s, well past even the slowest 80-bit generic
trial (395s), despite having the "fast" 2-factor (safe-prime-like) p-1 shape
that Benchmark B1 showed was consistently the *faster* shape at 60 bits.
That inversion is the headline result here: whatever makes the 2-factor
shape fast below the threshold does not obviously carry over once the
solver is in Wiedemann territory — but this is a **single trial**, so it's
reported as a real, reproducible-looking signal worth a dedicated follow-up
sweep (multiple trials at 82-90 bits) rather than a settled conclusion.

**Caveats.** Every number above comes from a fixed seed
(`common.BENCH_SEED = 20260705`) and small trial counts sized to fit a
~20-30 minute budget — treat this as a first rigorous pass, not a
statistically definitive study. The 80-bit and Benchmark B2 group sizes (1-2
trials) in particular are too small to distinguish real effects from single-
sample noise; the honest takeaway is "bit length dominates, p-1 shape has a
real secondary effect below the Wiedemann threshold, and something worth
investigating further happens right at/above it."
