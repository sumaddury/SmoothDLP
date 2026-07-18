#pragma once
#include <vector>
#include <utility>
#include <cstdint>
#include "types.h"

namespace gauss {

/**
  All primes up to MAX_SMALL_PRIME (1,000,000), precomputed once at load
  time via sieveTo -- the small-prime trial-division base used by
  factorize_naive/factorize, and the exact-lookup table isPrime uses below
  MAX_SMALL_PRIME.
*/
extern std::vector<uint32_t> full_primes_array;

/**
  Every prime <= n, in increasing order (sieve of Eratosthenes).
*/
std::vector<uint32_t> sieveTo(uint64_t n);

/**
  Primality test for n. Exact (bitset lookup against full_primes_array) for
  n <= MAX_SMALL_PRIME; above that, Miller-Rabin via Montgomery arithmetic
  (mont::mulmod) -- a fixed deterministic witness set below a ~2^64-ish
  bound, RANDOM_MR_ROUNDS random witnesses above it (so it becomes
  probabilistic, not exact, for very large n).
*/
bool isPrime(u128 n);

/**
  Trial-divides n by every prime in full_primes_array (up to sqrt(n), or up
  to full_primes_array's cap for n beyond MAX_NAIVE), returning
  (factors found so far, remaining unfactored cofactor). The remainder is
  1 iff n was fully factored by small primes alone -- otherwise it's n's
  largest cofactor with no small prime factors, still possibly composite.
*/
std::pair<FactorList, u128> factorize_naive(u128 n);

/**
  Shanks' Square Form Factorization (SQUFOF): returns a single nontrivial
  factor of composite n (not necessarily prime, not necessarily the
  smallest) by searching for a square form among several multiplier
  candidates (the table M). Returns 1 if no factor is found.
*/
u128 squfof(u128 n);

/**
  Full factorization of n into (prime, exponent) pairs, whose
  reconstruction (product of prime^exponent) equals n exactly. Combines
  factorize_naive (small primes), perfect square/cube stripping, and
  SQUFOF (for cofactors with no small factors), recursing until every
  remaining piece is confirmed prime via isPrime.
*/
FactorList factorize(u128 n);

} // namespace gauss
