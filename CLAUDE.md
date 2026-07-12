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
  montgomery.h/.cpp   128-bit modular-arithmetic backend (mul, reduce, Montgomery REDC) — implemented + fully tested, but not yet reachable from Python, see below
  gauss_dream.cpp/h   sieveTo, isPrime (Miller-Rabin via Montgomery), factorize_naive, squfof, factorize
  smooth_algos.cpp/h  isSmooth, logDickman (Dickman rho), mp_ln, log_mul, psiApprox (uses dickman_table.bin)
  infra.cpp/h         buildProductTree/smoothCandidates/treeFactorize (Bernstein batch-smoothness),
                      linSolve/crtSolve (LinBox sparse solve + CRT/Hensel over factors of p-1), rank_relation_gf2
  core.cpp            pybind11 bindings -> smooth._core
  __init__.py          Python-facing `smooth` package surface
tests/
  test_gauss_dream.py, test_infra.py, test_smooth_algos.py   pytest, run against the built extension
  conftest.py           reference/oracle implementations used by the above
  internal/             standalone C++ test suite for montgomery.{h,cpp} only (no LinBox/pybind11 needed)
    test_montgomery.cpp, Makefile  -- `make -C tests/internal test`
scripts/        gitignored scratch area: legacy debug artifacts (temp.txt, rank_test.txt,
                test_linbox_u64.cpp, data_collect.cpp, sweep_results.csv) from an earlier
                LinBox rank-deficiency investigation. Not part of the build; historical only.
texts/          reference PDFs. Only 2409.08784v1.pdf (the paper) is authoritative; rest is
                background reading or the stale original plan (see warning above).
notebooks/      prototyping notebooks (primality -> factoring -> smoothness), predate the C++ port
setup.py        builds smooth._core. Sources: core.cpp, gauss_dream.cpp, smooth_algos.cpp, infra.cpp
environment.yml conda env `smoothnumbers-dev`: gmp, ntl, givaro, fflas-ffpack, linbox, boost, pybind11
```

## Build & test

```bash
conda activate smoothnumbers-dev            # provides GMP/NTL/LinBox/Givaro/FFLAS-FFPACK
python setup.py build_ext --inplace          # builds src/_core.cpython-*.so
python -m pytest tests/ -q                   # Python-level tests (skipped if smooth fails to import)

make -C tests/internal test                  # standalone montgomery.{h,cpp} suite (GMP-only, no conda extension needed)
```

## Current status (2026-07-11, updated same day) — read before editing montgomery.*

**Update: `montgomery.h/.cpp` is now fully implemented and fully tested.**
`mulmod`, `pow2mod_odd`, and `pow2mod` all landed (the `pow2mod` family
collapsed from a planned three-piece odd/any/dispatcher split down to just
two functions — see "Design notes" below for the current shape). All 9
`mont::` functions have edge-case + randomized-vs-GMP-oracle tests in
`tests/internal/test_montgomery.cpp` (22 sections, 222k+ checks, all
passing via `make -C tests/internal test`). Points 1, 3, and 4 from the
original version of this status section are resolved; **point 2 (below) is
now the sole remaining blocker** for the whole project.

1. **`setup.py` still never lists `src/montgomery.cpp` as a build source —
   this is now the only reason the project doesn't build/import/test
   end-to-end.** The compiled `smooth._core` extension links "successfully"
   (macOS `-undefined dynamic_lookup` defers symbol resolution) but **fails
   at import time** with `symbol not found in flat namespace
   '__ZN4mont15r_squared_mod_nEo'` — every `mont::` symbol gauss_dream.cpp/
   infra.cpp reference is missing from the shared object, because
   `montgomery.cpp` itself is never compiled into it. `tests/conftest.py`
   uses `pytest.importorskip("smooth")`, but since the failure happens while
   *loading conftest itself*, pytest reports a collection **error**, not a
   skip — the whole suite still errors out at collection. **The fix is a
   one-line addition** to `setup.py`'s `sources=[...]` list; do this next.
2. Stray editor autosave files still sit untracked in `src/`
   (`##montgomery.cpp##`, `#montgomery.h#` — the emacs lock file and one
   autosave from earlier have since cleared on their own). Harmless, not
   part of the build, but consider adding an emacs-autosave `.gitignore`
   pattern so they stop showing up in `git status`.
3. **`infra.h` is missing declarations for `linSolve`/`crtSolve`/
   `rank_relation_gf2`** — all three are defined in `infra.cpp` but not
   declared in the header, not bound in `core.cpp`, and not covered by any
   test. They implement the CRT+Hensel-lift sparse solve over factors of
   `p-1` (LinBox `Wiedemann`/`SparseElimination`) and are load-bearing for the
   whole project's endgame, but currently dangling/unverified code with no
   way to exercise them short of writing a new C++ harness or wiring them
   into `core.cpp`.
4. **The actual "double" algorithm orchestration doesn't exist in `src/` at
   all**: no relation-collection loop, no generator-finder, no smoothness-
   bound heuristic, no g-side/b-side dual channel, no overlap-finder, no
   `x ≡ α·β⁻¹ (mod p-1)` step. (An earlier prototype of this — `dlp.cpp`,
   `bench_helpers.cpp` — existed before the "restructure and cleanup" commit
   and is gone now; `texts/PROJECT_BRIEFING.md` describes that old, now-
   nonexistent layout — treat it as historical context on the *math*, not
   as a map of current files.) This is the paper's actual contribution and is
   the biggest remaining chunk of work, but it's blocked on point 1 above —
   nothing downstream of montgomery is reachable from Python until that's
   fixed.

## Module status: what's ready / blocked / missing

(Note: this section tracks *implementation progress*, not call-frequency —
see "Design notes" below for the HOT/WARM/COLD performance classification.)

- **Done, but not yet reachable from Python — this is the active blocker**:
  - `montgomery.h/.cpp` — all 9 functions implemented and fully tested in
    isolation (`make -C tests/internal test`). Not yet added to `setup.py`'s
    `sources=[...]`, so `smooth._core` can't actually use it yet (see status
    point 1 above). This is now a build-config fix, not an algorithm task.

- **Implemented + has test coverage, but currently blocked/unverified until
  montgomery is wired into `setup.py`**:
  - `gauss_dream.cpp` (`isPrime`, `factorize_naive`, `squfof`, `factorize`,
    `sieveTo`) — fully written, has `tests/test_gauss_dream.py`, but `isPrime`
    calls `mont::mulmod` for any `n > 1,000,000`, so it's non-functional until
    montgomery is reachable (see status point 1).
  - `smooth_algos.cpp` (`isSmooth`, `logDickman`, `psiApprox`, `mp_ln`,
    `log_mul`) — fully written, has `tests/test_smooth_algos.py`, transitively
    blocked because `isSmooth` calls `isPrime`.
  - `infra.cpp::buildProductTree` / `treeFactorize` — done, tested
    (`tests/test_infra.py`), **not** blocked by montgomery (pure GMP product-
    tree/gcd work).
  - `infra.cpp::smoothCandidates` — implemented, tested, but directly calls
    `mont::pow2mod`, so it's blocked the same way `isPrime` is.
  - `infra.cpp::linSolve` / `crtSolve` / `rank_relation_gf2` — implemented
    (LinBox Wiedemann/SparseElimination + CRT/Hensel), but undeclared in
    `infra.h`, unbound in `core.cpp`, and **zero automated test coverage**.
    Correctness was never confirmed even in the prior iteration of this code
    (see `texts/PROJECT_BRIEFING.md` §9 for the historical rank-deficiency
    investigation this traces back to) — treat as unverified, not done.

- **Not started — no code exists yet**:
  - The dual-channel double-index-calculus driver itself: generator search,
    smoothness-bound selection, g-side/b-side relation collection loops, the
    overlap-finder comparing partial log tables, and the final
    `x ≡ α·β⁻¹ (mod p-1)` combination. This is the paper's actual novel
    contribution (§III of `2409.08784v1.pdf`) and doesn't exist anywhere in
    `src/` right now.

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
