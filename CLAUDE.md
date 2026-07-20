# SmoothNumbers

C++/pybind11 implementation of the **Double Index Calculus Algorithm** (Huang,
Zhang, Zhao, Peng, Liao, Wang, 2024 — `texts/2409.08784v1.pdf`), a discrete-log
solver for finite prime fields. Given prime `p`, generator-ish `g`, and target
`b`, the paper's trick is to collect smooth-number relations for `g` and `b`
*simultaneously* (instead of fully solving `g`'s side first): as soon as one
factor-base prime's discrete log is known on both sides, `x ≡ α·β⁻¹ (mod p-1)`
falls out directly, needing fewer relations than classic index calculus.

Research/educational scale only: the paper benchmarks 30–75 bit primes; this
project targets the same range, not real cryptographic sizes.

**Before touching math-heavy code, skim `texts/2409.08784v1.pdf` (the paper)
first.** Most other files in `texts/` (`SmoothNumbersPlan.pdf`, `canvas.pdf`,
and everything not the paper itself) are **stale planning docs from the
original Python-module design** (May 2025, package names like `gauss_dream.py`,
`discrete_log.py`) — the project was later rewritten from scratch in C++ with
different structure. Don't trust filenames/APIs from those PDFs; they describe
a plan that was abandoned, not the current code.

## Directory map

```
src/            C++ core, compiled into the smooth._core pybind11 extension
  types.h         u128 typedef, mpz_class<->u128 helpers, clz128/ctz128/popcount128
  montgomery.h/.cpp   128-bit modular-arithmetic backend (mul, reduce, Montgomery REDC) —
                      implemented, fully tested, and wired into setup.py; reachable from Python
                      transitively (gauss_dream/infra call it), not bound directly itself
  gauss_dream.cpp/h   sieveTo, isPrime (Miller-Rabin via Montgomery), factorize_naive, squfof, factorize
  smooth_algos.cpp/h  isSmooth, logDickman (Dickman rho), mp_ln, log_mul, psiApprox (uses dickman_table.bin)
  lin_alg.cpp/h       raw LinBox layer: solveModPrime/rankModPrime (GF(q) sparse solve/rank,
                      SparseElimination or Wiedemann) + henselLift (lifts a base GF(q) solve to
                      q^e via e-1 correction rounds). Not a setup.py source on its own -- pulled
                      into infra.cpp via #include "lin_alg.cpp" (unity build, see below); not
                      bound in core.cpp, so unreachable from Python.
  infra.cpp/h         buildProductTree/smoothCandidates/treeFactorize (Bernstein batch-smoothness,
                      bound directly in core.cpp) + ProblemParams/addRelations/crtSolve/
                      filterDetermined/modInv (the double-index-calculus relation-collection +
                      CRT/Hensel solve + Omega_g/Omega_b + final-combination layer, built on
                      lin_alg.cpp). Not bound as standalone Python entry points -- only reachable
                      from Python transitively through dlp::DLP (below).
  dlp.cpp/h           dlp::DLP -- the top-level driver that assembles infra.cpp's building blocks
                      into the paper's Algorithm 1 end to end (DLP(p).solve(g, b) -> x). Bound
                      directly in core.cpp as a Python class. Does NOT unity-include infra.cpp
                      (only includes infra.h), so it links against infra.cpp as an ordinary
                      separate translation unit -- see "Design notes: dlp::DLP and crtSolve's
                      concurrency constraint" below for why that distinction actually mattered
                      during development.
  core.cpp            pybind11 bindings -> smooth._core
  __init__.py          Python-facing `smooth` package surface
tests/
  test_gauss_dream.py, test_infra.py, test_smooth_algos.py, test_dlp.py   pytest, run against the
                        built extension
  conftest.py           reference/oracle implementations used by the above
  internal/             standalone C++ test suites, one per header with complete functions:
                        test_montgomery, test_gauss_dream, test_smooth_algos (GMP-only), plus
                        test_infra, test_lin_alg, and test_dlp (link LinBox/Givaro/NTL/FFLAS-FFPACK
                        too). `make -C tests/internal test` builds + runs all six. Independent of
                        setup.py/pybind11 -- these link straight against GMP/LinBox, no Python.
scripts/        gitignored scratch area: legacy debug/investigation artifacts from earlier LinBox
                rank-deficiency work. Not part of the build; historical only.
texts/          reference PDFs. Only 2409.08784v1.pdf (the paper) is authoritative; rest is
                background reading or the stale original plan (see warning above).
setup.py        builds smooth._core. Sources: core.cpp, gauss_dream.cpp, smooth_algos.cpp,
                infra.cpp, dlp.cpp, montgomery.cpp. (lin_alg.cpp is not listed -- it's
                unity-included by infra.cpp, see below, so it has no separate translation unit.)
environment.yml conda env `smoothnumbers-dev`: gmp, ntl, givaro, fflas-ffpack, linbox, boost, pybind11
```

### Unity-build requirement around `lin_alg.h`

Any two translation units that both transitively include `lin_alg.h`'s LinBox matrix/vector
headers and get linked into the same binary hit ~60 duplicate-symbol errors at link time —
LinBox's `commentator.h`/`args-parser.h` define several symbols out-of-line without `inline`.
The fix used throughout this codebase is `#include` the `.cpp` directly instead of the `.h`, so
everything stays in one translation unit: `infra.cpp` does `#include "lin_alg.cpp"` (not
`lin_alg.h`), and `tests/internal/test_lin_alg.cpp`/`test_infra.cpp` each `#include` their
matching `.cpp` for the same reason. Don't "clean this up" into separate compiled objects
without re-checking the link still succeeds.

### Compiler-selection gotcha

The plain `c++`/`g++` on this machine resolves to Homebrew's `g++-16` (libstdc++ ABI), which is
ABI-incompatible with the conda-forge-built LinBox/Givaro/NTL (libc++ ABI). `setup.py` already
handles this by pointing `CXX`/`CC` at the conda env's `*-apple-darwin*-clang++`. If you ever
compile anything by hand (outside `setup.py`/`tests/internal/Makefile`, which both do this
correctly already), use that same conda clang++, not a bare `c++`/`g++`.

## Build & test

```bash
conda activate smoothnumbers-dev            # provides GMP/NTL/LinBox/Givaro/FFLAS-FFPACK
python setup.py build_ext --inplace          # builds src/_core.cpython-*.so
python -m pytest tests/ -q                   # Python-level tests (skipped if smooth fails to import)

make -C tests/internal test                  # all 6 standalone C++ suites (needs the conda env
                                              # active for LinBox/Givaro/NTL -- see Makefile)
```

## Module status: what's ready / blocked / missing

The full pipeline described in the paper (§III of `2409.08784v1.pdf`) is now
implemented, tested, and reachable from Python end to end via `dlp::DLP` —
there is no longer a "not started" category. What remains is only the
narrower distinction between what's bound as its own Python entry point
versus what's only reachable transitively through `DLP.solve(...)`.

- **Done, tested, reachable from Python (bound in `core.cpp`)**:
  - `montgomery.h/.cpp` — all 9 functions, wired into `setup.py`'s
    `sources=[...]`. Not bound directly (no `mont::` Python entry points),
    but reachable transitively through everything below that calls it.
  - `gauss_dream.cpp` (`isPrime`, `factorize`, `sieveTo`) and
    `smooth_algos.cpp` (`isSmooth`, `logDickman`, `psiApprox`, `log_mul`) —
    bound and tested (`tests/test_gauss_dream.py` / `tests/test_smooth_algos.py`).
    `factorize_naive`, `squfof` (`gauss_dream.cpp`) and `mp_ln`
    (`smooth_algos.cpp`) remain fully implemented and tested **at the C++
    level only** (`tests/internal/test_gauss_dream.cpp` /
    `test_smooth_algos.cpp`) — deliberately removed from `core.cpp`/
    `__init__.py`, so `hasattr(smooth, "factorize_naive")` etc. is now
    `False`. Use `gauss::factorize_naive`/`gauss::squfof`/`salgo::mp_ln`
    directly from C++ if you need them; don't re-add pybind bindings for
    them without checking why they were pulled first.
  - `infra.cpp::buildProductTree` / `smoothCandidates` / `treeFactorize` —
    Bernstein batch-smoothness testing, tested (`tests/test_infra.py`).
  - `dlp::DLP` — `DLP(p).solve(g, b)` implements the paper's Algorithm 1 end
    to end: doubling-schedule relation collection (`infra::addRelations`),
    CRT/Hensel solve (`infra::crtSolve`), the Ω_g/Ω_b direct-verification
    construction (`infra::filterDetermined`), and the final
    `x ≡ α·β⁻¹ (mod p-1)` combination (`infra::modInv`). Bound in `core.cpp`
    via `py::class_<dlp::DLP>`; tested both at the C++ level
    (`tests/internal/test_dlp.cpp`) and via pytest (`tests/test_dlp.py`).
    Read "Design notes: dlp::DLP and crtSolve's concurrency constraint"
    below before touching its threading structure — it looks like ordinary
    `std::async` fork-join but has one load-bearing constraint that isn't
    obvious from the code alone.

- **Done and tested (`make -C tests/internal test`, see `test_lin_alg.cpp` /
  `test_infra.cpp`), but not bound as standalone Python entry points —
  reachable from Python only transitively through `dlp::DLP`**:
  - `lin_alg.cpp` (`solveModPrime`, `rankModPrime`, `henselLift`) — the raw
    GF(q) solve/rank layer plus Hensel lifting from a base solve up to q^e.
    See doc comments in `lin_alg.h` for the `SolveStatus`/`MAX_MODULUS`
    contract and why the field type is `Givaro::Modular<RecInt::ruint<7>>`
    specifically (a different-looking but "equivalent" LinBox field silently
    returns wrong solves above q ~ 2^32 — do not swap it, see the comment in
    `lin_alg.h` and the memory note on this).
  - `infra.cpp` (`ProblemParams`, `addRelations`, `crtSolve`,
    `filterDetermined`, `modInv`) — the relation-collection +
    CRT/Hensel-over-factors-of-(p-1) solve + overlap + combination layer,
    built on `lin_alg.cpp`. See doc comments in `infra.h` for the full
    contract, in particular: `crtSolve` returns an **empty vector** if
    `henselLift` fails on any single factor of p-1 (not an exception, not a
    partial result — check for empty), and `addRelations`' `base` argument
    has an **unchecked precondition** that it be a genuine primitive root
    mod p (see "Primitive-root precondition" below) — this is a real,
    common failure mode if violated, not a rare edge case. `dlp::DLP` is
    what actually calls all of these today; binding them individually as
    their own Python entry points (mirroring `buildProductTree`/
    `smoothCandidates`/`treeFactorize`) would only be needed for partial
    access (e.g. just relation collection, without a full solve) —
    `modInv` also needs its definition in `infra.cpp` to stay a normal
    (non-`inline`) function if it's ever touched again: it was briefly
    marked `inline` there, which is only safe for a function whose
    definition's translation unit also calls it (true of the file-local
    `modInvPrimePower`/`binExp` helpers, not true of `modInv` itself, a
    public API function nothing in `infra.cpp` calls) — with `inline` on
    it, the compiler is free to discard the symbol entirely whenever it's
    unreferenced in that one TU, which silently breaks linking for any
    consumer (like `dlp.cpp`) that links against `infra.cpp` as an ordinary
    separate object rather than unity-including it the way
    `test_infra.cpp` does.

## Design notes: dlp::DLP and crtSolve's concurrency constraint

`dlp::DLP::solve` runs its g-side and b-side pipelines
(`addRelations` + `crtSolve` + `filterDetermined`) concurrently via
`std::async`, but the two `crtSolve` calls are serialized against each other
through a file-local `static std::mutex` in `dlp.cpp` even though everything
else about the two sides runs in parallel. This was **not** a defensive
default — it was added after direct testing during development found that
two threads calling `crtSolve` concurrently crash the underlying
LinBox/Givaro/NTL/FFLAS-FFPACK stack nondeterministically (SIGABRT or
SIGTRAP, at unpredictable points — not reliably on the first call, sometimes
after dozens of successful calls), regardless of input. `addRelations` was
separately confirmed safe under the same kind of concurrent-call stress test
(pure GMP/Montgomery arithmetic on disjoint per-side state) — the
unsafety is specific to the LinBox-touching code path, not to concurrency in
this codebase generally.

The root cause was not fully diagnosed (most likely candidate: LinBox's
`commentator` singleton, or some other global/non-thread-local state
somewhere in that dependency chain — the same code that already causes the
unrelated duplicate-symbol issue documented under "Unity-build requirement"
above), but the empirical finding is solid and reproducible. The mutex is
intentionally a single **process-wide** static, not a per-`DLP`-instance
member — since the evidence points to genuinely global state in the
underlying library, not per-`ProblemParams` state, scoping the mutex down to
per-instance would silently reintroduce the crash for two `DLP` instances
(even over different primes) solving concurrently. Don't remove or
re-scope this mutex without re-running a concurrent-`crtSolve`-call stress
test first.

## Primitive-root precondition on `addRelations`' `base`

`addRelations(base, ...)` requires `base` to be a genuine primitive root mod
`p` (order exactly p-1), whether it's playing the "g" or "b" role — the two
roles go through the identical code path, confirmed by direct testing of
both. This is **not validated by the code** (matching this codebase's
general precondition philosophy: it is not `addRelations`' job to validate
its contract, the same way `std::lower_bound` does not verify its range is
sorted). If violated, `base^t mod p` only ever ranges over the proper
subgroup `base` generates, which caps what `crtSolve`'s linear system can
learn about factor-base primes outside that subgroup's structure. This
produces a **permanent rank ceiling on some factor of p-1** that persists no
matter how many relations are collected (confirmed directly: 100+ relations
past the minimum, still deficient) — distinct from ordinary statistical rank
deficiency, which is transient and reliably clears by roughly 3x the
factor-base size k (also confirmed directly, across multiple primes and
factor shapes; see `tests/internal/test_infra.cpp`'s k-relative
failure-rate/rank sweep). If `crtSolve` keeps returning empty well past a
few multiples of k, check whether `base` is actually a primitive root before
assuming more relations will help.

## Design notes: the `pow2mod` family

`montgomery.h` ended up with **two** functions in this family, not the
originally-anticipated three (`pow2mod_odd`/`pow2mod_any`/`pow2mod`
dispatcher) — the "any `n`" case and the dispatcher merged into one CRT-based
`pow2mod`:

- **`pow2mod_odd(base, bits, n, n_prime, r2)`** — computes
  `base^(2^bits) mod n` for odd `n` via Montgomery-form repeated squaring
  (`bits` squarings, not square-and-multiply — the exponent is always exactly
  a power of two). `n_prime = mont::inverse(n)`, `r2 = mont::r_squared_mod_n(n)`.
- **`pow2mod(base, bits, n, m_prime, r2_m)`** — computes the same quantity for
  **any** `n` (odd or even). Splits `n = 2^e · m` (`e = ctz(n)`, `m` odd),
  solves the `2^e` part with a cheap bitmask (no division at all), solves the
  `m` part by calling `pow2mod_odd` directly, then recombines with CRT
  (Garner's formula). **`m_prime`/`r2_m` are relative to `m` (n's odd part),
  not to `n` itself** — see `pow2mod_params_for()` in
  `tests/internal/test_montgomery.cpp` for exactly how a caller derives them.
  `n == 1` is special-cased to return `0` directly.

**Why the extra parameters (`n_prime`/`r2`/`m_prime`/`r2_m`) exist at all:**
`r_squared_mod_n` internally calls the comparatively expensive `mulmod_any`
(full 256-bit division via `reduce256`). Rather than have `pow2mod_odd`/
`pow2mod` recompute it internally on every call, the caller precomputes it
once and passes it in — the same hoisting pattern `gauss_dream.cpp::isPrime`
already uses (it computes `n_prime`/`r2` once outside its per-witness
Miller-Rabin loop, not once per witness).

**Why `pow2mod` doesn't need a separate dispatcher anymore:** the CRT split
handles both parities uniformly — when `n` is odd, `e = 0`, the masked part
is trivial, and the whole computation degenerates to exactly
`pow2mod_odd(base, bits, n, ...)` (this is asserted directly by
`test_pow2mod_matches_pow2mod_odd_for_odd_n` in the test suite). So there was
never a need for callers to branch on parity themselves — `pow2mod` alone
handles every `n`, and `pow2mod_odd` is exposed separately only because
`pow2mod` (and any other future caller that already knows its modulus is odd)
benefits from skipping the `ctz`/CRT machinery entirely.

### Hot path — where the optimization effort actually went

The one real call site today is `infra.cpp::smoothCandidates`:
`mont::pow2mod(r, bits, x)`, invoked once per candidate in a batch — batches
are the entire point of Bernstein's trick (potentially thousands to millions
of candidates per relation-collection run). Candidates `x` are essentially
uniform over `[1, p-1]` for a large prime `p`, so the even/odd split is close
to 50/50 — **both cases run regularly on every batch**, not just the odd one.
COLD would mean "not regularly hit in the critical path" — nothing in this
family qualifies:

- **HOT — `pow2mod_odd`**. The core workhorse: true Montgomery REDC squaring.
  Hit directly on every odd candidate, *and* indirectly as the inner engine
  for the even case too (via `pow2mod`'s CRT split) — so in practice it's
  exercised on close to 100% of candidates one way or another. This is why
  it's written as a tight loop with no redundant conversions in/out of
  Montgomery form.
- **WARM — `pow2mod`**. The even-`n` case is handled here, not with a dumb
  "square `bits` times via `mulmod_any`" loop — precisely because that case
  is regularly hit (~half of candidates), not rare, and deserved real design
  effort. Since `e = ctz(n)` is small most of the time
  (`P(e=k) ≈ 2^{-(k+1)}`), the odd cofactor `m` is usually nearly as large as
  `n`, so routing it through `pow2mod_odd` means the even case costs about as
  much as the odd case — most of it *is* the odd case underneath, via
  delegation rather than a from-scratch reimplementation.
