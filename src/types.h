#pragma once
#include <cstdint>
#include <vector>
#include <utility>
#include <gmpxx.h>

using u128 = unsigned __int128;
using U128Vector = std::vector<u128>;
using MpzVector = std::vector<mpz_class>;
using SparseList = std::vector<std::pair<size_t, uint32_t>>;    // (column index, exponent/coefficient)
using RelationMatrix = std::vector<SparseList>;                  // one SparseList per row
using FactorList = std::vector<std::pair<u128, uint32_t>>;       // (prime, exponent) pairs

/**
 * Converts a u128 to the equivalent arbitrary-precision mpz_class, by
 * splitting it into 64-bit high/low halves (GMP has no native 128-bit
 * constructor) and recombining as hi*2^64 + lo.
 */
inline mpz_class u128_to_mpz(u128 v) {
    mpz_class hi((unsigned long)(uint64_t)(v >> 64));
    mpz_class lo((unsigned long)(uint64_t)(v));
    return (hi << 64) + lo;
}

/**
 * Converts an mpz_class back to a u128, taking the low 128 bits (via a
 * 64-bit high/low split, the inverse of u128_to_mpz). Undefined/truncating
 * if m doesn't fit in 128 bits or is negative.
 */
inline u128 mpz_to_u128(const mpz_class& m) {
    mpz_class tmp = m;
    uint64_t lo = mpz_get_ui(tmp.get_mpz_t());
    mpz_class shifted = tmp >> 64;
    uint64_t hi = mpz_get_ui(shifted.get_mpz_t());
    return ((u128)hi << 64) | (u128)lo;
}

/**
 * Count of leading zero bits in x, over the full 128-bit width (clz128(0)
 * == 128). Equivalently, 128 - bit_length(x).
 */
inline int clz128(u128 x) {
    uint64_t hi = (uint64_t)(x >> 64);
    if (hi != 0) return __builtin_clzll(hi);
    uint64_t lo = (uint64_t)x;
    return lo != 0 ? 64 + __builtin_clzll(lo) : 128;
}

/**
 * Count of trailing zero bits in x (ctz128(0) == 128). Equivalently, the
 * largest k such that 2^k divides x -- e.g. x >> ctz128(x) is x's odd part.
 */
inline int ctz128(u128 x) {
    uint64_t lo = (uint64_t)x;
    if (lo != 0) return __builtin_ctzll(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return (hi != 0) ? 64 + __builtin_ctzll(hi) : 128;
}

/**
 * Number of set (1) bits in x, over the full 128-bit width.
 */
inline int popcount128(u128 x) {
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);

    return __builtin_popcountll(lo) + __builtin_popcountll(hi);
}
