// Standalone correctness tests for src/montgomery.{h,cpp}.
//
// montgomery has no dependencies on the rest of the project (gauss_dream,
// infra, LinBox, ...), so this test suite only compiles/links montgomery.cpp
// against GMP (already pulled in by types.h for the u128<->mpz_class helpers,
// which double as this suite's oracle). See the Makefile in this directory
// for the exact build command.
//
// All mont:: functions declared in montgomery.h are implemented and tested
// here: mul128, gcd, mulmod_any, inverse, r_squared_mod_n, reduce256, mulmod,
// pow2mod_odd, pow2mod, powmod_odd.
//
// Note the current shape of the pow2mod family (there is no separate
// "pow2mod_any" anymore -- it was folded into pow2mod itself):
//   - pow2mod_odd(base, bits, n, n_prime, r2) computes base^(2^bits) mod n
//     for odd n via Montgomery-form repeated squaring.
//   - pow2mod(base, bits, n, m_prime, r2_m) computes base^(2^bits) mod n for
//     ANY n (odd or even) by splitting n = 2^e * m (m odd, e = ctz(n)),
//     solving the 2^e part with a cheap bitmask (no division), solving the m
//     part by calling pow2mod_odd, and recombining with CRT (Garner). Its
//     m_prime/r2_m parameters are relative to *m* (n's odd part), not to n.
// In both cases, n_prime/m_prime = mont::inverse(n or m) and
// r2/r2_m = mont::r_squared_mod_n(n or m) are precomputed by the caller and
// passed in, rather than recomputed internally on every call -- this avoids
// repeated, expensive mulmod_any-based r_squared_mod_n calls when the same
// modulus is reused across several pow2mod*/mulmod calls (mirroring how
// gauss_dream.cpp's isPrime already hoists n_prime/r2 out of its
// per-witness Miller-Rabin loop).

#include "montgomery.h"
#include "types.h"

#include <gmpxx.h>

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

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

#define CHECK_EQ_U128(actual, expected, ctx) \
    do { \
        ++g_checks; \
        u128 _a = (actual); \
        u128 _e = (expected); \
        if (_a != _e) { \
            report_failure(__FILE__, __LINE__, \
                std::string(ctx) + ": got " + u128_to_mpz(_a).get_str() + \
                ", expected " + u128_to_mpz(_e).get_str()); \
        } \
    } while (0)

constexpr u128 U128_MAX = ~(u128)0;

std::mt19937_64 rng(0xC0FFEEULL);

u128 random_u128() {
    uint64_t hi = rng();
    uint64_t lo = rng();
    return ((u128)hi << 64) | (u128)lo;
}

// Random value with at most `bits` significant bits (bits in [1, 128]).
u128 random_u128_bits(int bits) {
    u128 v = random_u128();
    if (bits < 128) v &= (((u128)1 << bits) - 1);
    return v;
}

// ---- GMP oracles -----------------------------------------------------

void mul128_oracle(u128 a, u128 b, u128& hi, u128& lo) {
    mpz_class A = u128_to_mpz(a);
    mpz_class B = u128_to_mpz(b);
    mpz_class P = A * B;
    mpz_class mask = (mpz_class(1) << 128) - 1;
    mpz_class LO = P & mask;
    mpz_class HI = P >> 128;
    lo = mpz_to_u128(LO);
    hi = mpz_to_u128(HI);
}

u128 gcd_oracle(u128 a, u128 b) {
    mpz_class A = u128_to_mpz(a);
    mpz_class B = u128_to_mpz(b);
    mpz_class G;
    mpz_gcd(G.get_mpz_t(), A.get_mpz_t(), B.get_mpz_t());
    return mpz_to_u128(G);
}

u128 mulmod_any_oracle(u128 a, u128 b, u128 n) {
    mpz_class A = u128_to_mpz(a);
    mpz_class B = u128_to_mpz(b);
    mpz_class N = u128_to_mpz(n);
    mpz_class P = (A * B) % N;
    return mpz_to_u128(P);
}

u128 r_squared_mod_n_oracle(u128 n) {
    mpz_class N = u128_to_mpz(n);
    mpz_class R = mpz_class(1) << 128;
    mpz_class result = (R * R) % N;
    return mpz_to_u128(result);
}

u128 reduce256_oracle(u128 hi, u128 lo, u128 n) {
    mpz_class HI = u128_to_mpz(hi);
    mpz_class LO = u128_to_mpz(lo);
    mpz_class N = u128_to_mpz(n);
    mpz_class value = (HI << 128) + LO;
    mpz_class result = value % N;
    return mpz_to_u128(result);
}

// ---- mont::mul128 ------------------------------------------------------

void test_mul128_edge_cases() {
    struct Case { u128 a, b; };
    Case cases[] = {
        {0, 0},
        {0, 12345},
        {12345, 0},
        {1, 1},
        {1, U128_MAX},
        {U128_MAX, 1},
        {U128_MAX, U128_MAX},
        {(u128)1 << 127, 2},
        {(u128)1 << 127, (u128)1 << 127},
        {((u128)1 << 64) - 1, ((u128)1 << 64) - 1},
        {(u128)1 << 64, (u128)1 << 64},
    };
    for (auto& c : cases) {
        u128 hi, lo, ehi, elo;
        mont::mul128(c.a, c.b, hi, lo);
        mul128_oracle(c.a, c.b, ehi, elo);
        CHECK_EQ_U128(lo, elo, "mul128 lo edge case");
        CHECK_EQ_U128(hi, ehi, "mul128 hi edge case");
    }
}

void test_mul128_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        u128 a = random_u128();
        u128 b = random_u128();
        u128 hi, lo, ehi, elo;
        mont::mul128(a, b, hi, lo);
        mul128_oracle(a, b, ehi, elo);
        CHECK_EQ_U128(lo, elo, "mul128 lo random");
        CHECK_EQ_U128(hi, ehi, "mul128 hi random");
    }
}

// mul128 is commutative: a*b == b*a as a 256-bit product.
void test_mul128_commutative() {
    constexpr int kIters = 5000;
    for (int i = 0; i < kIters; ++i) {
        u128 a = random_u128();
        u128 b = random_u128();
        u128 hi1, lo1, hi2, lo2;
        mont::mul128(a, b, hi1, lo1);
        mont::mul128(b, a, hi2, lo2);
        CHECK_EQ_U128(lo1, lo2, "mul128 commutative lo");
        CHECK_EQ_U128(hi1, hi2, "mul128 commutative hi");
    }
}

// ---- mont::gcd -----------------------------------------------------------

void test_gcd_edge_cases() {
    CHECK_EQ_U128(mont::gcd(0, 0), 0, "gcd(0,0)");
    CHECK_EQ_U128(mont::gcd(0, 5), 5, "gcd(0,5)");
    CHECK_EQ_U128(mont::gcd(5, 0), 5, "gcd(5,0)");
    CHECK_EQ_U128(mont::gcd(1, U128_MAX), 1, "gcd(1,MAX)");
    CHECK_EQ_U128(mont::gcd(17, 13), 1, "gcd(17,13) coprime");
    CHECK_EQ_U128(mont::gcd((u128)1 << 100, (u128)1 << 50), (u128)1 << 50,
                  "gcd(2^100,2^50)");
    CHECK_EQ_U128(mont::gcd(U128_MAX, U128_MAX), U128_MAX, "gcd(MAX,MAX)");

    for (u128 a : {(u128)6, (u128)123456789, (u128)1 << 90}) {
        CHECK_EQ_U128(mont::gcd(a, a), a, "gcd(a,a)");
    }
}

void test_gcd_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        // Bias toward smaller bit-widths sometimes so shared factors of 2
        // actually show up, not just uniformly-random near-coprime pairs.
        int bits_a = 1 + (int)(rng() % 128);
        int bits_b = 1 + (int)(rng() % 128);
        u128 a = random_u128_bits(bits_a);
        u128 b = random_u128_bits(bits_b);
        u128 got = mont::gcd(a, b);
        u128 want = gcd_oracle(a, b);
        CHECK_EQ_U128(got, want, "gcd random");
    }
}

// ---- mont::mulmod_any -----------------------------------------------------

void test_mulmod_any_edge_cases() {
    CHECK_EQ_U128(mont::mulmod_any(12345, 6789, 1), 0, "mulmod_any n=1");
    CHECK_EQ_U128(mont::mulmod_any(0, 6789, 97), 0, "mulmod_any a=0");
    CHECK_EQ_U128(mont::mulmod_any(6789, 0, 97), 0, "mulmod_any b=0");

    // n fits in 64 bits (exercises the single-limb reduction path).
    CHECK_EQ_U128(mont::mulmod_any(123456789, 987654321, 1000000007),
                  mulmod_any_oracle(123456789, 987654321, 1000000007),
                  "mulmod_any small odd n");
    CHECK_EQ_U128(mont::mulmod_any(123456789, 987654321, 100000000),
                  mulmod_any_oracle(123456789, 987654321, 100000000),
                  "mulmod_any small even n");

    // n > 2^64, top bit already set (no normalization shift needed, d == 0).
    u128 n_hi_bit = ((u128)1 << 127) | 12345;
    CHECK_EQ_U128(mont::mulmod_any(U128_MAX, U128_MAX - 1, n_hi_bit),
                  mulmod_any_oracle(U128_MAX, U128_MAX - 1, n_hi_bit),
                  "mulmod_any n with high bit set");

    // n > 2^64, top bit clear (exercises the normalization-shift branch).
    u128 n_shifted = (u128)1 << 100;
    CHECK_EQ_U128(mont::mulmod_any(U128_MAX, U128_MAX, n_shifted),
                  mulmod_any_oracle(U128_MAX, U128_MAX, n_shifted),
                  "mulmod_any n needing normalization");

    // Even modulus > 2^64.
    u128 n_even_big = ((u128)1 << 90) + 6;
    CHECK_EQ_U128(mont::mulmod_any(U128_MAX, 42, n_even_big),
                  mulmod_any_oracle(U128_MAX, 42, n_even_big),
                  "mulmod_any large even n");

    // a, b unreduced (larger than n).
    u128 n_small = 97;
    CHECK_EQ_U128(mont::mulmod_any(U128_MAX, U128_MAX, n_small),
                  mulmod_any_oracle(U128_MAX, U128_MAX, n_small),
                  "mulmod_any unreduced operands");
}

void test_mulmod_any_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        u128 a = random_u128();
        u128 b = random_u128();
        int n_bits = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(n_bits);
        if (n == 0) n = 1;
        u128 got = mont::mulmod_any(a, b, n);
        u128 want = mulmod_any_oracle(a, b, n);
        CHECK_EQ_U128(got, want, "mulmod_any random");
    }
}

// ---- mont::reduce256 -------------------------------------------------
//
// Contract, as currently implemented: reduce256(hi, lo, n) == (hi:lo, a
// 256-bit value) mod n, for any n > 0 -- this is a plain full reduction
// (mulmod_any is now just mul128 + reduce256), *not* Montgomery REDC (no
// n_prime involved). Tested directly with arbitrary hi/lo pairs, not just
// ones reachable as a product of two u128 factors, since e.g. hi close to
// 2^128-1 combined with a large lo exceeds what mul128(a,b) can ever
// produce (max product is (2^128-1)^2 < 2^256 - 2^129).

void test_reduce256_edge_cases() {
    CHECK_EQ_U128(mont::reduce256(0, 0, 97), 0, "reduce256(0,0,n)");
    CHECK_EQ_U128(mont::reduce256(0, 41, 97), 41, "reduce256 already-reduced lo");
    CHECK_EQ_U128(mont::reduce256(0, 96, 97), 96, "reduce256 lo == n-1");
    CHECK_EQ_U128(mont::reduce256(0, 97, 97), 0, "reduce256 lo == n");
    CHECK_EQ_U128(mont::reduce256(U128_MAX, U128_MAX, 1), 0, "reduce256 n=1");

    // n fits in 64 bits (single-limb reduction path).
    CHECK_EQ_U128(mont::reduce256(123456789, 987654321, 1000000007),
                  reduce256_oracle(123456789, 987654321, 1000000007),
                  "reduce256 small odd n");
    CHECK_EQ_U128(mont::reduce256(123456789, 987654321, 100000000),
                  reduce256_oracle(123456789, 987654321, 100000000),
                  "reduce256 small even n");

    // n > 2^64, top bit already set (d == 0, no normalization shift).
    u128 n_hi_bit = ((u128)1 << 127) | 12345;
    CHECK_EQ_U128(mont::reduce256(U128_MAX, U128_MAX - 1, n_hi_bit),
                  reduce256_oracle(U128_MAX, U128_MAX - 1, n_hi_bit),
                  "reduce256 n with high bit set");

    // n > 2^64, top bit clear (exercises the normalization-shift branch).
    u128 n_shifted = (u128)1 << 100;
    CHECK_EQ_U128(mont::reduce256(U128_MAX, U128_MAX, n_shifted),
                  reduce256_oracle(U128_MAX, U128_MAX, n_shifted),
                  "reduce256 n needing normalization");

    // Even modulus > 2^64.
    u128 n_even_big = ((u128)1 << 90) + 6;
    CHECK_EQ_U128(mont::reduce256(U128_MAX, 42, n_even_big),
                  reduce256_oracle(U128_MAX, 42, n_even_big),
                  "reduce256 large even n");

    // hi alone already exceeds anything mul128 could ever produce for this n.
    u128 n_small = 97;
    CHECK_EQ_U128(mont::reduce256(U128_MAX, U128_MAX, n_small),
                  reduce256_oracle(U128_MAX, U128_MAX, n_small),
                  "reduce256 hi beyond mul128 range");

    // Exact multiple of n (remainder should land on exactly 0), built by
    // splitting k*n across the hi:lo boundary.
    {
        mpz_class N = u128_to_mpz((u128)((((u128)1) << 100) + 3));
        mpz_class K = u128_to_mpz(U128_MAX - 12345);
        mpz_class value = N * K;
        mpz_class mask = (mpz_class(1) << 128) - 1;
        u128 lo = mpz_to_u128(value & mask);
        u128 hi = mpz_to_u128(value >> 128);
        CHECK_EQ_U128(mont::reduce256(hi, lo, mpz_to_u128(N)), 0,
                      "reduce256 exact multiple of n");
    }
}

void test_reduce256_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        u128 hi = random_u128();
        u128 lo = random_u128();
        int n_bits = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(n_bits);
        if (n == 0) n = 1;
        u128 got = mont::reduce256(hi, lo, n);
        u128 want = reduce256_oracle(hi, lo, n);
        CHECK_EQ_U128(got, want, "reduce256 random");
    }
}

// ---- mont::inverse ---------------------------------------------------
//
// Contract: for odd n, x = mont::inverse(n) satisfies n*x == 1 (mod 2^128).
// u128 multiplication already wraps mod 2^128 natively, so this needs no
// GMP oracle -- just check the defining property directly. (Note this is
// the *positive* inverse, n*x == 1, not the negated n' with n*n' == -1
// mod 2^128 that REDC conventionally calls "n_prime" -- worth double
// checking against whichever convention reduce256/mulmod end up using.)

void test_inverse_edge_cases() {
    u128 odd_values[] = {
        1,
        3,
        5,
        U128_MAX,           // all-ones, odd
        U128_MAX - 2,       // odd, near max
        ((u128)1 << 127) + 1,
        ((u128)1 << 64) + 1,
        ((u128)1 << 64) - 1,
    };
    for (u128 n : odd_values) {
        u128 inv = mont::inverse(n);
        CHECK_EQ_U128(n * inv, 1, "inverse edge case: n * inverse(n) == 1");
    }
}

void test_inverse_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        u128 n = random_u128() | 1; // force odd
        u128 inv = mont::inverse(n);
        CHECK_EQ_U128(n * inv, 1, "inverse random: n * inverse(n) == 1");
    }
}

// inverse(n) should itself be odd (it's a unit mod 2^128, and the only
// units mod a power of two are odd numbers).
void test_inverse_is_odd() {
    constexpr int kIters = 5000;
    for (int i = 0; i < kIters; ++i) {
        u128 n = random_u128() | 1;
        u128 inv = mont::inverse(n);
        CHECK((inv & 1) == 1);
    }
}

// ---- mont::r_squared_mod_n ---------------------------------------------
//
// Contract: r_squared_mod_n(n) == R^2 mod n, where R = 2^128 -- the value
// needed to convert an integer into Montgomery form. Depends on
// mulmod_any (already verified above), so this cross-checks the whole
// thing against GMP computing R^2 mod n directly with a real 2^128
// literal (no overflow tricks needed on the oracle side).

void test_r_squared_mod_n_edge_cases() {
    u128 values[] = {
        1,
        2,
        3,
        U128_MAX,               // R == n + 1, so R^2 mod n == 1
        U128_MAX - 1,
        (u128)1 << 64,
        ((u128)1 << 64) + 1,
        (u128)1 << 127,
        ((u128)1 << 127) + 1,
    };
    for (u128 n : values) {
        u128 got = mont::r_squared_mod_n(n);
        u128 want = r_squared_mod_n_oracle(n);
        CHECK_EQ_U128(got, want, "r_squared_mod_n edge case");
    }
    CHECK_EQ_U128(mont::r_squared_mod_n(U128_MAX), 1, "r_squared_mod_n(2^128-1) == 1");
}

void test_r_squared_mod_n_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        int bits = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits);
        if (n == 0) n = 1;
        u128 got = mont::r_squared_mod_n(n);
        u128 want = r_squared_mod_n_oracle(n);
        CHECK_EQ_U128(got, want, "r_squared_mod_n random");
    }
}

// ---- mont::mulmod (Montgomery REDC) -----------------------------------
//
// Contract: mulmod(a_bar, b_bar, n, n_prime) computes (a_bar * b_bar * R^-1)
// mod n, where R = 2^128 and n_prime = mont::inverse(n) (the *positive*
// inverse, n * n_prime == 1 mod R -- mulmod negates it internally to get the
// REDC constant n' with n * n' == -1 mod R, so callers always pass
// mont::inverse(n) directly, never a pre-negated value). Requires n odd, and
// a_bar, b_bar already reduced (< n) -- that's how every call site in this
// codebase (gauss_dream.cpp's Miller-Rabin, pow2mod_odd below) uses it, and
// it's also the precondition classic REDC needs (T = a_bar*b_bar < n*R) to
// guarantee the single conditional subtraction fully reduces the result.

u128 mulmod_oracle(u128 a_bar, u128 b_bar, u128 n) {
    mpz_class A = u128_to_mpz(a_bar);
    mpz_class B = u128_to_mpz(b_bar);
    mpz_class N = u128_to_mpz(n);
    mpz_class R = mpz_class(1) << 128;
    mpz_class R_inv;
    mpz_invert(R_inv.get_mpz_t(), R.get_mpz_t(), N.get_mpz_t());
    mpz_class P = ((A * B) % N) * R_inv % N;
    return mpz_to_u128(P);
}

void test_mulmod_edge_cases() {
    // n = 1: everything reduces to 0.
    {
        u128 n = 1;
        u128 n_prime = mont::inverse(n);
        CHECK_EQ_U128(mont::mulmod(0, 0, n, n_prime), 0, "mulmod n=1");
    }
    // a_bar = 0 / b_bar = 0.
    {
        u128 n = 97;
        u128 n_prime = mont::inverse(n);
        CHECK_EQ_U128(mont::mulmod(0, 42, n, n_prime), 0, "mulmod a_bar=0");
        CHECK_EQ_U128(mont::mulmod(42, 0, n, n_prime), 0, "mulmod b_bar=0");
    }
    // n odd, near 2^128.
    {
        u128 n = U128_MAX;
        u128 n_prime = mont::inverse(n);
        u128 a = n - 1, b = n - 2;
        CHECK_EQ_U128(mont::mulmod(a, b, n, n_prime), mulmod_oracle(a, b, n),
                      "mulmod large odd n near 2^128");
    }
    // Smallest odd n > 1.
    {
        u128 n = 3;
        u128 n_prime = mont::inverse(n);
        CHECK_EQ_U128(mont::mulmod(2, 2, n, n_prime), mulmod_oracle(2, 2, n),
                      "mulmod n=3");
    }
}

void test_mulmod_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        int bits = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits) | 1; // force odd
        u128 n_prime = mont::inverse(n);
        u128 a = random_u128() % n;
        u128 b = random_u128() % n;
        u128 got = mont::mulmod(a, b, n, n_prime);
        u128 want = mulmod_oracle(a, b, n);
        CHECK_EQ_U128(got, want, "mulmod random");
    }
}

void test_mulmod_result_in_range() {
    constexpr int kIters = 5000;
    for (int i = 0; i < kIters; ++i) {
        int bits = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits) | 1;
        u128 n_prime = mont::inverse(n);
        u128 a = random_u128() % n;
        u128 b = random_u128() % n;
        u128 got = mont::mulmod(a, b, n, n_prime);
        CHECK(got < n);
    }
}

// ---- mont::pow2mod_odd -------------------------------------------------
//
// Contract: pow2mod_odd(base, bits, n, n_prime, r2) computes
// base^(2^bits) mod n for odd n, via `bits` Montgomery-form squarings (not
// square-and-multiply over an arbitrary exponent -- the exponent here is
// always exactly a power of two, so plain repeated squaring is exact).
// n_prime = mont::inverse(n), r2 = mont::r_squared_mod_n(n): both computed
// once by the caller and passed in (see file-header note above).

u128 pow2mod_oracle(u128 base, unsigned bits, u128 n) {
    mpz_class B = u128_to_mpz(base) % u128_to_mpz(n);
    mpz_class N = u128_to_mpz(n);
    mpz_class E = mpz_class(1) << bits;
    mpz_class P;
    mpz_powm(P.get_mpz_t(), B.get_mpz_t(), E.get_mpz_t(), N.get_mpz_t());
    return mpz_to_u128(P);
}

void test_pow2mod_odd_edge_cases() {
    struct Case { u128 base; unsigned bits; u128 n; };
    Case cases[] = {
        {5, 0, 97},                     // bits=0: base^1 mod n, no squarings
        {0, 10, 97},                    // base = 0
        {96, 5, 97},                    // base = n-1
        {2, 1, 3},                      // smallest odd n > 1
        {U128_MAX - 2, 20, U128_MAX},   // n odd, near 2^128
        {7, 127, 97},                   // large bits
    };
    for (auto& c : cases) {
        u128 n_prime = mont::inverse(c.n);
        u128 r2 = mont::r_squared_mod_n(c.n);
        u128 got = mont::pow2mod_odd(c.base, c.bits, c.n, n_prime, r2);
        u128 want = pow2mod_oracle(c.base, c.bits, c.n);
        CHECK_EQ_U128(got, want, "pow2mod_odd edge case");
    }
}

void test_pow2mod_odd_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n) | 1;
        u128 n_prime = mont::inverse(n);
        u128 r2 = mont::r_squared_mod_n(n);
        u128 base = random_u128();
        unsigned bits = (unsigned)(rng() % 40);
        u128 got = mont::pow2mod_odd(base, bits, n, n_prime, r2);
        u128 want = pow2mod_oracle(base, bits, n);
        CHECK_EQ_U128(got, want, "pow2mod_odd random");
    }
}

// ---- mont::pow2mod (general n, via CRT split n = 2^e * m) --------------
//
// Contract: pow2mod(base, bits, n, m_prime, r2_m) computes
// base^(2^bits) mod n for ANY n >= 1. Internally splits n = 2^e * m
// (e = ctz(n), m odd), solves mod 2^e with a bitmask, solves mod m via
// pow2mod_odd, and recombines with CRT/Garner. Crucially, m_prime/r2_m are
// relative to *m* (n's odd part), not to n -- pow2mod_params_for below
// mirrors exactly what a caller must compute.

struct PowModParams { u128 m_prime, r2_m; };

PowModParams pow2mod_params_for(u128 n) {
    int e = ctz128(n);
    u128 m = n >> e;
    return {mont::inverse(m), mont::r_squared_mod_n(m)};
}

void test_pow2mod_edge_cases() {
    CHECK_EQ_U128(mont::pow2mod(12345, 10, 1, 0, 0), 0, "pow2mod n=1");

    struct Case { u128 base; unsigned bits; u128 n; };
    Case cases[] = {
        {5, 10, 2},                  // n=2: pure power of two (m=1)
        {5, 10, 1024},               // n=2^10: pure power of two (m=1)
        {5, 10, 97},                 // n odd: reduces to pow2mod_odd (e=0)
        {5, 10, U128_MAX},           // n odd, near 2^128
        {123, 8, 96},                // n = 2^5 * 3: mixed
        {123, 8, (u128)3 << 100},    // large mixed n
        {0, 5, 50},                  // base = 0
    };
    for (auto& c : cases) {
        auto p = pow2mod_params_for(c.n);
        u128 got = mont::pow2mod(c.base, c.bits, c.n, p.m_prime, p.r2_m);
        u128 want = pow2mod_oracle(c.base, c.bits, c.n);
        CHECK_EQ_U128(got, want, "pow2mod edge case");
    }
}

void test_pow2mod_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n);
        if (n == 0) n = 1;
        u128 base = random_u128();
        unsigned bits = (unsigned)(rng() % 40);

        auto p = pow2mod_params_for(n);
        u128 got = mont::pow2mod(base, bits, n, p.m_prime, p.r2_m);
        u128 want = pow2mod_oracle(base, bits, n);
        CHECK_EQ_U128(got, want, "pow2mod random");
    }
}

// pow2mod must agree with pow2mod_odd whenever n happens to be odd (e=0,
// the CRT split degenerates to just the odd cofactor).
void test_pow2mod_matches_pow2mod_odd_for_odd_n() {
    constexpr int kIters = 2000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n) | 1;
        u128 base = random_u128();
        unsigned bits = (unsigned)(rng() % 40);

        u128 n_prime = mont::inverse(n);
        u128 r2 = mont::r_squared_mod_n(n);
        u128 want = mont::pow2mod_odd(base, bits, n, n_prime, r2);

        auto p = pow2mod_params_for(n);
        u128 got = mont::pow2mod(base, bits, n, p.m_prime, p.r2_m);
        CHECK_EQ_U128(got, want, "pow2mod matches pow2mod_odd for odd n");
    }
}

// ---- mont::powmod_odd (arbitrary-exponent square-and-multiply) --------
//
// Contract: powmod_odd(base, e, n, n_prime, r2) computes base^e mod n for
// odd n and an arbitrary exponent e (unlike pow2mod_odd, e need not be a
// power of two -- this does real square-and-multiply, not just repeated
// squaring). n_prime = mont::inverse(n), r2 = mont::r_squared_mod_n(n),
// same hoisting convention as pow2mod_odd. Note `e` is a plain `unsigned`
// (32 bits), not u128 like every modulus/base in this file -- so this
// currently cannot represent exponents needing more than 32 bits.

u128 powmod_oracle(u128 base, unsigned e, u128 n) {
    mpz_class B = u128_to_mpz(base) % u128_to_mpz(n);
    mpz_class N = u128_to_mpz(n);
    mpz_class E = e;
    mpz_class P;
    mpz_powm(P.get_mpz_t(), B.get_mpz_t(), E.get_mpz_t(), N.get_mpz_t());
    return mpz_to_u128(P);
}

void test_powmod_odd_edge_cases() {
    struct Case { u128 base; unsigned e; u128 n; };
    Case cases[] = {
        {5, 0, 97},                        // e=0: base^0 == 1 mod n
        {0, 0, 97},                        // base=0, e=0: 0^0 == 1 by convention
        {0, 5, 97},                        // base=0, e>0: 0^e == 0
        {5, 1, 97},                        // e=1: base mod n
        {96, 5, 97},                       // base = n-1 (i.e. -1 mod n), odd e
        {96, 6, 97},                       // base = n-1, even e
        {2, 10, 1024 * 1024 * 1024 + 7},   // not a power-of-two exponent
        {2, 1, 3},                         // smallest odd n > 1
        {U128_MAX - 2, 12345, U128_MAX},   // n odd, near 2^128
        {7, 0xFFFFFFFFu, 97},              // e at the top of unsigned's range
    };
    for (auto& c : cases) {
        u128 n_prime = mont::inverse(c.n);
        u128 r2 = mont::r_squared_mod_n(c.n);
        u128 got = mont::powmod_odd(c.base, c.e, c.n, n_prime, r2);
        u128 want = powmod_oracle(c.base, c.e, c.n);
        CHECK_EQ_U128(got, want, "powmod_odd edge case");
    }
}

void test_powmod_odd_random() {
    constexpr int kIters = 20000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n) | 1;
        u128 n_prime = mont::inverse(n);
        u128 r2 = mont::r_squared_mod_n(n);
        u128 base = random_u128();
        unsigned e = (unsigned)rng(); // exercise the full unsigned range
        u128 got = mont::powmod_odd(base, e, n, n_prime, r2);
        u128 want = powmod_oracle(base, e, n);
        CHECK_EQ_U128(got, want, "powmod_odd random");
    }
}

void test_powmod_odd_result_in_range() {
    constexpr int kIters = 5000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n) | 1;
        u128 n_prime = mont::inverse(n);
        u128 r2 = mont::r_squared_mod_n(n);
        u128 base = random_u128();
        unsigned e = (unsigned)rng();
        u128 got = mont::powmod_odd(base, e, n, n_prime, r2);
        CHECK(got < n);
    }
}

// powmod_odd must agree with pow2mod_odd whenever e is exactly a power of
// two -- two independently-shaped algorithms (square-and-multiply over an
// arbitrary exponent vs. plain repeated squaring over a bits count) that
// should land on exactly the same value.
void test_powmod_odd_matches_pow2mod_odd_for_power_of_two_e() {
    constexpr int kIters = 2000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n) | 1;
        u128 n_prime = mont::inverse(n);
        u128 r2 = mont::r_squared_mod_n(n);
        u128 base = random_u128();
        unsigned bits = (unsigned)(rng() % 32); // keep 2^bits within unsigned range

        u128 want = mont::pow2mod_odd(base, bits, n, n_prime, r2);
        u128 got = mont::powmod_odd(base, (unsigned)1 << bits, n, n_prime, r2);
        CHECK_EQ_U128(got, want, "powmod_odd matches pow2mod_odd for e = 2^bits");
    }
}

// Exponent law: base^(e1) * base^(e2) == base^(e1+e2) mod n. A property
// check independent of the GMP oracle, using mulmod to combine two
// powmod_odd results in Montgomery-ish ordinary form.
void test_powmod_odd_exponent_law() {
    constexpr int kIters = 5000;
    for (int i = 0; i < kIters; ++i) {
        int bits_n = 1 + (int)(rng() % 128);
        u128 n = random_u128_bits(bits_n) | 1;
        u128 n_prime = mont::inverse(n);
        u128 r2 = mont::r_squared_mod_n(n);
        u128 base = random_u128();

        // Keep e1+e2 within unsigned range to avoid wraparound changing the
        // mathematical exponent being compared.
        unsigned e1 = (unsigned)(rng() % (1u << 20));
        unsigned e2 = (unsigned)(rng() % (1u << 20));

        u128 r1 = mont::powmod_odd(base, e1, n, n_prime, r2);
        u128 r2_val = mont::powmod_odd(base, e2, n, n_prime, r2);
        u128 combined = mont::mulmod_any(r1, r2_val, n);

        u128 want = mont::powmod_odd(base, e1 + e2, n, n_prime, r2);
        CHECK_EQ_U128(combined, want, "powmod_odd exponent law: base^e1 * base^e2 == base^(e1+e2)");
    }
}

struct TestCase {
    const char* name;
    void (*fn)();
};

} // namespace

int main() {
    TestCase tests[] = {
        {"mul128_edge_cases", test_mul128_edge_cases},
        {"mul128_random", test_mul128_random},
        {"mul128_commutative", test_mul128_commutative},
        {"gcd_edge_cases", test_gcd_edge_cases},
        {"gcd_random", test_gcd_random},
        {"mulmod_any_edge_cases", test_mulmod_any_edge_cases},
        {"mulmod_any_random", test_mulmod_any_random},
        {"reduce256_edge_cases", test_reduce256_edge_cases},
        {"reduce256_random", test_reduce256_random},
        {"inverse_edge_cases", test_inverse_edge_cases},
        {"inverse_random", test_inverse_random},
        {"inverse_is_odd", test_inverse_is_odd},
        {"r_squared_mod_n_edge_cases", test_r_squared_mod_n_edge_cases},
        {"r_squared_mod_n_random", test_r_squared_mod_n_random},
        {"mulmod_edge_cases", test_mulmod_edge_cases},
        {"mulmod_random", test_mulmod_random},
        {"mulmod_result_in_range", test_mulmod_result_in_range},
        {"pow2mod_odd_edge_cases", test_pow2mod_odd_edge_cases},
        {"pow2mod_odd_random", test_pow2mod_odd_random},
        {"pow2mod_edge_cases", test_pow2mod_edge_cases},
        {"pow2mod_random", test_pow2mod_random},
        {"pow2mod_matches_pow2mod_odd_for_odd_n", test_pow2mod_matches_pow2mod_odd_for_odd_n},
        {"powmod_odd_edge_cases", test_powmod_odd_edge_cases},
        {"powmod_odd_random", test_powmod_odd_random},
        {"powmod_odd_result_in_range", test_powmod_odd_result_in_range},
        {"powmod_odd_matches_pow2mod_odd_for_power_of_two_e", test_powmod_odd_matches_pow2mod_odd_for_power_of_two_e},
        {"powmod_odd_exponent_law", test_powmod_odd_exponent_law},
    };

    for (auto& t : tests) {
        int before = g_failures;
        t.fn();
        std::cout << (g_failures == before ? "[PASS] " : "[FAIL] ") << t.name
                   << "\n";
    }

    std::cout << "\n" << (g_checks - g_failures) << "/" << g_checks
               << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
