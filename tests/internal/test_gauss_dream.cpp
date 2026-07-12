// Standalone correctness tests for the complete, headered functions in
// src/gauss_dream.{h,cpp}: sieveTo, isPrime, factorize_naive, squfof,
// factorize. Depends on montgomery.h (isPrime's Miller-Rabin runs through
// mont::mulmod), so links montgomery.cpp too -- GMP only, otherwise, used
// here purely as an independent oracle (mpz_probab_prime_p, trial division),
// never as part of the code under test. See the Makefile in this directory
// for the exact build command (`make -C tests/internal test`).

#include "gauss_dream.h"
#include "types.h"

#include <gmpxx.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

int g_checks = 0;
int g_failures = 0;

void report_failure(const char* file, int line, const std::string& msg) {
    ++g_failures;
    std::cerr << "FAIL " << file << ":" << line << ": " << msg << "\n";
}

#define CHECK(cond) \
    do { \
        ++g_checks; \
        if (!(cond)) report_failure(__FILE__, __LINE__, #cond); \
    } while (0)

std::mt19937_64 rng(0x600DF00DULL);

u128 random_u128() {
    uint64_t hi = rng(), lo = rng();
    return ((u128)hi << 64) | (u128)lo;
}

// Random value with exactly `bits` significant bits (top bit set), odd.
u128 random_bitlength(int bits) {
    u128 v = random_u128();
    if (bits < 128) v &= (((u128)1 << bits) - 1);
    v |= (u128)1 << (bits - 1);
    return v | 1;
}

bool reference_is_prime(u128 n) {
    if (n < 2) return false;
    mpz_class m = u128_to_mpz(n);
    return mpz_probab_prime_p(m.get_mpz_t(), 25) != 0; // 1=probable, 2=definite
}

u128 next_probable_prime(u128 n) {
    n |= 1;
    while (!reference_is_prime(n)) n += 2;
    return n;
}

// ---- gauss_dream::sieveTo ------------------------------------------------
//
// Contract: sieveTo(n) returns every prime <= n, in increasing order.

void test_sieve_to_known_small_case() {
    auto primes = sieveTo(30);
    std::vector<uint32_t> expected = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    CHECK(primes == expected);
}

void test_sieve_to_matches_reference() {
    for (uint64_t n : {10ULL, 97ULL, 1000ULL, 50000ULL}) {
        std::vector<bool> composite(n + 1, false);
        std::vector<uint32_t> expected;
        for (uint64_t i = 2; i <= n; ++i) {
            if (composite[i]) continue;
            expected.push_back((uint32_t)i);
            for (uint64_t j = i * i; j <= n; j += i) composite[j] = true;
        }
        auto got = sieveTo(n);
        CHECK(got == expected);
    }
}

// ---- gauss_dream::isPrime -------------------------------------------------
//
// Contract: isPrime(n) is a primality test -- exact for n small enough to
// consult the precomputed sieve table, Miller-Rabin (deterministic witness
// set below a fixed 64-bit-ish bound, randomized above it) otherwise.

void test_is_prime_known_primes_and_composites() {
    for (u128 p : {2u, 3u, 5u, 7u, 997u, 7919u, 104729u}) CHECK(isPrime(p));
    for (u128 c : {0u, 1u, 4u, 6u, 100u, 999u, 1000000u}) CHECK(!isPrime(c));
}

void test_is_prime_exhaustive_small_range_vs_reference() {
    for (u128 n = 0; n < 200000; ++n) {
        bool got = isPrime(n);
        bool want = reference_is_prime(n);
        if (got != want) {
            report_failure(__FILE__, __LINE__,
                "isPrime mismatch at n=" + u128_to_mpz(n).get_str());
        }
        ++g_checks;
    }
}

void test_is_prime_known_large_mersenne_primes() {
    CHECK(isPrime(((u128)1 << 61) - 1));
    CHECK(isPrime(((u128)1 << 127) - 1));
}

void test_is_prime_known_large_composite() {
    CHECK(!isPrime((((u128)1 << 61) - 1) * 3));
}

void test_is_prime_randomized_vs_reference() {
    for (int bits : {40, 64, 96, 127}) {
        for (int i = 0; i < 100; ++i) {
            u128 n = random_bitlength(bits);
            bool got = isPrime(n);
            bool want = reference_is_prime(n);
            if (got != want) {
                report_failure(__FILE__, __LINE__,
                    "isPrime randomized mismatch at n=" + u128_to_mpz(n).get_str() +
                    " (bits=" + std::to_string(bits) + ")");
            }
            ++g_checks;
        }
    }
}

// ---- gauss_dream::factorize_naive -----------------------------------------
//
// Contract: factorize_naive(n) trial-divides by small primes only, returning
// (factors found so far, remaining unfactored part). remainder == 1 iff n
// was fully factored by trial division alone.

mpz_class reconstruct(const std::vector<std::pair<u128, uint32_t>>& factors) {
    mpz_class acc = 1;
    for (auto& [p, e] : factors)
        for (uint32_t i = 0; i < e; ++i) acc *= u128_to_mpz(p);
    return acc;
}

void test_factorize_naive_known_small_composite() {
    auto [factors, remainder] = factorize_naive(360);
    CHECK(remainder == 1);
    CHECK(reconstruct(factors) == 360);
}

void test_factorize_naive_small_times_large_prime() {
    u128 mersenne61 = ((u128)1 << 61) - 1;
    auto [factors, remainder] = factorize_naive(2 * mersenne61);
    CHECK(remainder == 1);
    CHECK(reconstruct(factors) == u128_to_mpz(2 * mersenne61));
}

void test_factorize_naive_large_semiprime_no_small_factors() {
    u128 p = 1099511628211ULL, q = 2199023255579ULL;
    CHECK(reference_is_prime(p) && reference_is_prime(q));
    u128 n = p * q;
    auto [factors, remainder] = factorize_naive(n);
    CHECK(factors.empty());
    CHECK(remainder == n);
}

void test_factorize_naive_randomized_reconstruction() {
    for (int iter = 0; iter < 30; ++iter) {
        u128 n = 1;
        int primes_used = 2 + (int)(rng() % 3);
        for (int i = 0; i < primes_used; ++i) {
            u128 candidate = random_bitlength(12 + (int)(rng() % 10));
            u128 p = next_probable_prime(candidate);
            // Keep the product within u128 range.
            mpz_class check = u128_to_mpz(n) * u128_to_mpz(p);
            if (check >= (mpz_class(1) << 127)) break;
            n *= p;
        }
        auto [factors, remainder] = factorize_naive(n);
        mpz_class reconstructed = reconstruct(factors) * u128_to_mpz(remainder);
        if (reconstructed != u128_to_mpz(n)) {
            report_failure(__FILE__, __LINE__,
                "factorize_naive randomized reconstruction mismatch for n=" +
                u128_to_mpz(n).get_str());
        }
        ++g_checks;
    }
}

// ---- gauss_dream::squfof --------------------------------------------------
//
// Contract: squfof(n) (Shanks' Square Form Factorization) returns a single
// nontrivial factor of composite n (not necessarily prime, not necessarily
// the smallest).

void test_squfof_returns_nontrivial_factor_randomized() {
    for (int iter = 0; iter < 20; ++iter) {
        u128 p = next_probable_prime(random_bitlength(20 + (int)(rng() % 15)));
        u128 q = next_probable_prime(random_bitlength(20 + (int)(rng() % 15)));
        if (p == q) continue;
        mpz_class n_mp = u128_to_mpz(p) * u128_to_mpz(q);
        if (n_mp >= (mpz_class(1) << 127)) continue;
        u128 n = p * q;

        u128 f = squfof(n);
        bool ok = (f != 1) && (f != n) && (n % f == 0);
        if (!ok) {
            report_failure(__FILE__, __LINE__,
                "squfof failed to find a nontrivial factor of n=" +
                u128_to_mpz(n).get_str() + " (got f=" + u128_to_mpz(f).get_str() + ")");
        }
        ++g_checks;
    }
}

// ---- gauss_dream::factorize ------------------------------------------------
//
// Contract: factorize(n) fully factors n (trial division + SQUFOF for
// larger cofactors), returning (prime, exponent) pairs whose reconstruction
// equals n exactly, with every reported prime genuinely prime.

void test_factorize_known_small_composite() {
    auto factors = factorize(360);
    std::sort(factors.begin(), factors.end());
    std::vector<std::pair<u128, uint32_t>> expected = {{2, 3}, {3, 2}, {5, 1}};
    CHECK(factors == expected);
}

void test_factorize_known_large_semiprime() {
    u128 p = 1099511628211ULL, q = 2199023255579ULL;
    u128 n = p * q;
    auto factors = factorize(n);
    CHECK(reconstruct(factors) == u128_to_mpz(n));
    std::vector<u128> primes_found;
    for (auto& [pr, e] : factors) {
        primes_found.push_back(pr);
        CHECK(e == 1);
    }
    std::sort(primes_found.begin(), primes_found.end());
    std::vector<u128> expected = {p, q};
    std::sort(expected.begin(), expected.end());
    CHECK(primes_found == expected);
}

void test_factorize_randomized_reconstruction_and_primality() {
    for (int iter = 0; iter < 30; ++iter) {
        u128 n = 1;
        int primes_used = 2 + (int)(rng() % 3);
        for (int i = 0; i < primes_used; ++i) {
            u128 candidate = random_bitlength(12 + (int)(rng() % 10));
            u128 p = next_probable_prime(candidate);
            mpz_class check = u128_to_mpz(n) * u128_to_mpz(p);
            if (check >= (mpz_class(1) << 127)) break;
            n *= p;
        }
        auto factors = factorize(n);
        if (reconstruct(factors) != u128_to_mpz(n)) {
            report_failure(__FILE__, __LINE__,
                "factorize reconstruction mismatch for n=" + u128_to_mpz(n).get_str());
        }
        ++g_checks;
        for (auto& [p, e] : factors) {
            (void)e;
            CHECK(reference_is_prime(p));
        }
    }
}

struct TestCase {
    const char* name;
    void (*fn)();
};

} // namespace

int main() {
    TestCase tests[] = {
        {"sieve_to_known_small_case", test_sieve_to_known_small_case},
        {"sieve_to_matches_reference", test_sieve_to_matches_reference},
        {"is_prime_known_primes_and_composites", test_is_prime_known_primes_and_composites},
        {"is_prime_exhaustive_small_range_vs_reference", test_is_prime_exhaustive_small_range_vs_reference},
        {"is_prime_known_large_mersenne_primes", test_is_prime_known_large_mersenne_primes},
        {"is_prime_known_large_composite", test_is_prime_known_large_composite},
        {"is_prime_randomized_vs_reference", test_is_prime_randomized_vs_reference},
        {"factorize_naive_known_small_composite", test_factorize_naive_known_small_composite},
        {"factorize_naive_small_times_large_prime", test_factorize_naive_small_times_large_prime},
        {"factorize_naive_large_semiprime_no_small_factors", test_factorize_naive_large_semiprime_no_small_factors},
        {"factorize_naive_randomized_reconstruction", test_factorize_naive_randomized_reconstruction},
        {"squfof_returns_nontrivial_factor_randomized", test_squfof_returns_nontrivial_factor_randomized},
        {"factorize_known_small_composite", test_factorize_known_small_composite},
        {"factorize_known_large_semiprime", test_factorize_known_large_semiprime},
        {"factorize_randomized_reconstruction_and_primality", test_factorize_randomized_reconstruction_and_primality},
    };

    for (auto& t : tests) {
        int before = g_failures;
        t.fn();
        std::cout << (g_failures == before ? "[PASS] " : "[FAIL] ") << t.name << "\n";
    }

    std::cout << "\n" << (g_checks - g_failures) << "/" << g_checks
               << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
