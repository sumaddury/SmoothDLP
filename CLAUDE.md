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
                      bound in core.cpp) + ProblemParams/addRelations/crtSolve (the double-index-
                      calculus relation-collection + CRT/Hensel solve layer, built on lin_alg.cpp,
                      not yet bound in core.cpp -- see "Module status" below)
  core.cpp            pybind11 bindings -> smooth._core
  __init__.py          Python-facing `smooth` package surface
tests/
  test_gauss_dream.py, test_infra.py, test_smooth_algos.py   pytest, run against the built extension
  conftest.py           reference/oracle implementations used by the above
  internal/             standalone C++ test suites, one per header with complete functions:
                        test_montgomery, test_gauss_dream, test_smooth_algos (GMP-only), plus
                        test_infra and test_lin_alg (link LinBox/Givaro/NTL/FFLAS-FFPACK too).
                        `make -C tests/internal test` builds + runs all five. Independent of
                        setup.py/pybind11 -- these link straight against GMP/LinBox, no Python.
scripts/        gitignored scratch area: legacy debug/investigation artifacts from earlier LinBox
                rank-deficiency work. Not part of the build; historical only.
texts/          reference PDFs. Only 2409.08784v1.pdf (the paper) is authoritative; rest is
                background reading or the stale original plan (see warning above).
setup.py        builds smooth._core. Sources: core.cpp, gauss_dream.cpp, smooth_algos.cpp,
                infra.cpp, montgomery.cpp. (lin_alg.cpp is not listed -- it's unity-included by
                infra.cpp, see below, so it has no separate translation unit.)
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

make -C tests/internal test                  # all 5 standalone C++ suites (needs the conda env
                                              # active for LinBox/Givaro/NTL -- see Makefile)
```

## Module status: what's ready / blocked / missing

- **Done, tested, reachable from Python (bound in `core.cpp`)**:
  - `montgomery.h/.cpp` — all 9 functions, wired into `setup.py`'s
    `sources=[...]`. Not bound directly (no `mont::` Python entry points),
    but reachable transitively through everything below that calls it.
  - `gauss_dream.cpp` (`isPrime`, `factorize_naive`, `squfof`, `factorize`,
    `sieveTo`) and `smooth_algos.cpp` (`isSmooth`, `logDickman`, `psiApprox`,
    `mp_ln`, `log_mul`) — fully working end-to-end, `tests/test_gauss_dream.py`
    / `tests/test_smooth_algos.py` pass against the built extension.
  - `infra.cpp::buildProductTree` / `smoothCandidates` / `treeFactorize` —
    Bernstein batch-smoothness testing, tested (`tests/test_infra.py`).

- **Done and tested (`make -C tests/internal test`, see `test_lin_alg.cpp` /
  `test_infra.cpp`), but not yet bound in `core.cpp` — unreachable from
  Python**:
  - `lin_alg.cpp` (`solveModPrime`, `rankModPrime`, `henselLift`) — the raw
    GF(q) solve/rank layer plus Hensel lifting from a base solve up to q^e.
    See doc comments in `lin_alg.h` for the `SolveStatus`/`MAX_MODULUS`
    contract and why the field type is `Givaro::Modular<RecInt::ruint<7>>`
    specifically (a different-looking but "equivalent" LinBox field silently
    returns wrong solves above q ~ 2^32 — do not swap it, see the comment in
    `lin_alg.h` and the memory note on this).
  - `infra.cpp` (`ProblemParams`, `addRelations`, `crtSolve`) — the
    relation-collection + CRT/Hensel-over-factors-of-(p-1) solve layer, built
    on `lin_alg.cpp`. See doc comments in `infra.h` for the full contract,
    in particular: `crtSolve` returns an **empty vector** if `henselLift`
    fails on any single factor of p-1 (not an exception, not a partial
    result — check for empty), and `addRelations`' `base` argument has an
    **unchecked precondition** that it be a genuine primitive root mod p (see
    "Primitive-root precondition" below) — this is a real, common failure
    mode if violated, not a rare edge case.
  - Binding this layer into `core.cpp` (mirroring the `buildProductTree`/
    `smoothCandidates`/`treeFactorize` pattern already there) is a
    reasonably mechanical next step whenever Python-side access is needed.

- **Not started — no code exists yet**:
  - The dual-channel double-index-calculus driver itself: generator search,
    smoothness-bound selection, g-side/b-side relation collection loops
    running concurrently, the overlap-finder comparing partial log tables
    across both sides, and the final `x ≡ α·β⁻¹ (mod p-1)` combination. This
    is the paper's actual novel contribution (§III of `2409.08784v1.pdf`)
    and is the biggest remaining chunk of work — `addRelations`/`crtSolve`
    are the building blocks it will be assembled from, one side at a time,
    but nothing wires the two sides together yet.

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
