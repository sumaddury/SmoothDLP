#pragma once
#include <cstdint>
#include <vector>
#include <utility>
#include <gmpxx.h>

using u128 = unsigned __int128;
using U128Vector = std::vector<u128>;
using MpzVector = std::vector<mpz_class>;
using SparseList = std::vector<std::pair<size_t, uint32_t>>;
using RelationMatrix = std::vector<SparseList>;
using FactorList = std::vector<std::pair<u128, uint32_t>>;

inline mpz_class u128_to_mpz(u128 v) {
    mpz_class hi((unsigned long)(uint64_t)(v >> 64));
    mpz_class lo((unsigned long)(uint64_t)(v));
    return (hi << 64) + lo;
}

inline u128 mpz_to_u128(const mpz_class& m) {
    mpz_class tmp = m;
    uint64_t lo = mpz_get_ui(tmp.get_mpz_t());
    mpz_class shifted = tmp >> 64;
    uint64_t hi = mpz_get_ui(shifted.get_mpz_t());
    return ((u128)hi << 64) | (u128)lo;
}

inline int clz128(u128 x) {
    uint64_t hi = (uint64_t)(x >> 64);
    if (hi != 0) return __builtin_clzll(hi);
    uint64_t lo = (uint64_t)x;
    return lo != 0 ? 64 + __builtin_clzll(lo) : 128;
}

inline int ctz128(u128 x) {
    uint64_t lo = (uint64_t)x;
    if (lo != 0) return __builtin_ctzll(lo);
    uint64_t hi = (uint64_t)(x >> 64);
    return (hi != 0) ? 64 + __builtin_ctzll(hi) : 128;
}

inline int popcount128(u128 x) {
    uint64_t lo = (uint64_t)x;
    uint64_t hi = (uint64_t)(x >> 64);

    return __builtin_popcountll(lo) + __builtin_popcountll(hi);
}
