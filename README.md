# SmoothNumbers

A from-scratch C++ implementation (exposed to Python via pybind11) of the
**Double Index Calculus Algorithm**, a method for solving the discrete
logarithm problem in finite prime fields, from:

> Huang, Zhang, Zhao, Peng, Liao, Wang. *"Double Index Calculus Algorithm:
> Faster Solving Discrete Logarithm Problem in Finite Prime Field."*
> arXiv:2409.08784 (2024).

This is a research/educational implementation. It targets the same scale the
paper itself benchmarks against (roughly 30–75 bit primes) — not real-world
cryptographic key sizes.

## The discrete logarithm problem

Fix a large prime `p`. The nonzero residues `{1, 2, ..., p-1}` under
multiplication mod `p` form a group. Pick a base `g` in that group. Given a
target `b`, the **discrete logarithm problem (DLP)** asks for the exponent
`x` such that:

```
g^x ≡ b (mod p)
```

Computing `g^x mod p` from `x` is cheap (fast exponentiation by repeated
squaring). Going backwards — recovering `x` from `g^x mod p` — is believed to
be computationally hard for a well-chosen large `p`. That asymmetry is the
security foundation underneath Diffie–Hellman key exchange, ElGamal
encryption, DSA signatures, and related schemes. This project is about
attacking DLP faster than the naive/generic approach, as a piece of
cryptanalysis research.

## Index calculus, briefly

The classical approach for this problem over prime fields is **index
calculus**. The idea: instead of attacking `x` directly, first fix a small
bound `B` and let the **factor base** be the set of primes `≤ B`. An integer
is `B`-smooth if all of its prime factors are `≤ B`.

1. Repeatedly try random exponents `t` and check whether `g^t mod p` is
   `B`-smooth. Each hit gives a linear relation between `t` and the (unknown)
   discrete logarithms of the factor-base primes, modulo `p - 1`.
2. Once enough relations are collected, solve the resulting linear system to
   recover `log_g(p_i)` for every prime `p_i` in the factor base.
3. Find one more smooth relation involving the actual target `b`, which turns
   directly into `log_g(b) = x`.

Classical index calculus therefore needs discrete logarithms for the whole
factor base — `k + 1` of them in total, where `k` is the factor-base size —
before it can answer the original question.

## The "double" idea

The paper this project implements changes the shape of the computation. Run
**two** relation-collection processes at once:

- **g-side:** sample random `t_i`, keep the ones where `g^{t_i} mod p` is
  smooth, and use them to build a linear system whose solution gives
  `log_g(p_i)` for some subset of factor-base primes.
- **b-side:** simultaneously, sample random `t̄_i`, keep the ones where
  `b^{t̄_i} mod p` is smooth, and solve for `log_b(p̄_i)` — logarithms of
  (possibly different) factor-base primes, but taken **base `b`** instead of
  base `g`.

Both sides can be partially solved as soon as their linear system reaches
full rank, even before every prime in the factor base has a known log on
that side. The moment **any single prime** `p̄` has a known logarithm on
*both* sides — `α = log_g(p̄)` and `β = log_b(p̄)` — the answer falls out
algebraically:

```
g^α ≡ p̄ ≡ b^β (mod p)          both sides equal the same prime
     ≡ (g^x)^β (mod p)          since b = g^x, by definition of x
⟹   α ≡ x·β (mod p - 1)
⟹   x ≡ α · β⁻¹ (mod p - 1)     the answer, with no further work
```

By the pigeonhole principle, once `(g-side logs found) + (b-side logs found)
≥ k + 1`, an overlap between the two sides is *guaranteed* — but in practice
one shows up far sooner (with 5 logs found on each side, the paper shows the
overlap probability already exceeds 94%). So the double algorithm typically
needs noticeably fewer than `k + 1` total discrete logarithms, and the two
independent channels can run in parallel. It also has a generality advantage:
classical index calculus requires `g` to actually generate the group, or
some factor-base logs simply won't exist; the double algorithm still works
with a non-generator base, since it only ever needs a partial, overlapping
set of logs.

Two supporting ideas make the "collect relations" step fast enough to be
practical at all:

- **Dickman's ρ function** estimates how common `B`-smooth numbers are near
  a given size, which tells you how many random samples to expect per
  smooth hit — turning "try random exponents until one works" into a
  plannable batch size instead of an open-ended guess.
- **Batch smoothness testing** (Bernstein's product-tree method) tests an
  entire batch of candidates against the whole factor base at once, instead
  of trial-dividing each candidate one at a time — the difference between
  testing thousands of candidates in roughly `O(batch · polylog)` time
  versus `O(batch × |factor base|)`.

## Data representation

- **Fixed-width integers (`u128`, an unsigned 128-bit type):** the values
  that show up on the arithmetic hot path — candidates being smoothness
  tested, moduli, exponents — fit comfortably in 128 bits at this project's
  target scale, so they're kept as native machine words rather than
  arbitrary-precision numbers wherever possible. This keeps the innermost
  loops (modular multiplication, exponentiation, gcd) working with plain
  CPU arithmetic and a hand-written Montgomery-multiplication backend,
  rather than paying for general-purpose bignum overhead on every step.
- **Arbitrary precision (GMP's `mpz_class`):** used where a value can
  genuinely outgrow 128 bits (products accumulated across a whole factor
  base, intermediate results in the CRT/Hensel-lifting step of the linear
  solve) or where the convenience of a full bignum library matters more than
  raw speed. Conversions between the two representations are a fixed cost
  at the boundary between "fast path" and "needs arbitrary precision."
- **Factor base as a product tree:** rather than a flat list of primes, the
  factor base is organized as a binary tree of products (leaves are the
  primes themselves; each internal node is the product of its two
  children; the root is the product of the entire factor base). This
  structure is what makes batch smoothness testing possible.
- **Relations as sparse rows:** a relation coming out of a smooth
  factorization is stored as a sparse list of `(factor-base index,
  exponent)` pairs rather than a dense row of mostly-zero coefficients —
  the natural shape given that any one smooth number only involves a
  handful of small primes out of a factor base that can have thousands of
  entries. Collections of such rows become the sparse matrix later handed
  to the linear solver.

## Python interface

The C++ core compiles into `smooth._core`, re-exported by the `smooth`
package. The functions currently exposed are organized loosely by the stage
of the pipeline they belong to:

**Number theory primitives**
- `sieve_to(n)` — primes up to `n`, via a sieve of Eratosthenes.
- `is_prime(n)` — primality test (deterministic for smaller `n`, probabilistic Miller–Rabin above that).
- `factorize_naive(n)` — trial-division factorization; returns `(factors, remainder)`, where `remainder` is whatever wasn't factored by small primes alone.
- `squfof(n)` — Shanks' Square Form Factorization; returns one nontrivial factor of `n`.
- `factorize(n)` — full factorization of `n`, combining trial division and SQUFOF.

**Smoothness and the Dickman function**
- `is_smooth(x, y)` — whether `x` is `y`-smooth (all prime factors ≤ `y`).
- `log_dickman(u)` — `ln(ρ(u))`, the log of Dickman's ρ function, used to estimate how common smooth numbers are.
- `mp_ln(x)` — natural log of a (possibly very large) integer `x`.
- `log_mul(x, log_rho)` — computes `x * exp(log_rho)` for large `x`, without losing precision to floating point along the way.
- `psi_approx(x, y)` — an estimate of `Ψ(x, y)`, the count of `y`-smooth integers up to `x`.

**Batch smoothness testing (Bernstein's product-tree method)**
- `build_product_tree(values)` — builds the product-tree levels over a list of values (typically the factor base).
- `smooth_candidates(p_levels, X)` — given a product tree and a batch of candidates `X`, returns the indices of the candidates that are smooth over the tree's base set.
- `tree_factorize(p_levels, d)` — factors a known-smooth `d` over the product tree, returning its sparse `(index, exponent)` relation row.

A typical exploratory session looks like:

```python
import smooth

smooth.is_prime(2**61 - 1)                 # True
smooth.factorize(360)                      # [(2, 3), (3, 2), (5, 1)]

primes = smooth.sieve_to(50)               # factor base candidates
tree = smooth.build_product_tree(primes)

candidates = [2*3*5*7*11, 41*43*47, 59*61]
smooth_idx = smooth.smooth_candidates(tree, candidates)   # which candidates are smooth
smooth.tree_factorize(tree, candidates[smooth_idx[0]])    # its relation row
```

The pieces that turn these primitives into an actual end-to-end discrete-log
solve — generator search, the g-side/b-side relation collectors, the
overlap-finder, and the final combination step — are still under active
development and not yet part of the public interface.
