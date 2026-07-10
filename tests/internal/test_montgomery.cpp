// Standalone correctness tests for src/montgomery.{h,cpp}.
//
// montgomery has no dependencies on the rest of the project (gauss_dream,
// infra, LinBox, ...), so this test suite only compiles/links montgomery.cpp
// against GMP (already pulled in by types.h for the u128<->mpz_class helpers,
// which double as this suite's oracle). See the Makefile in this directory
// for the exact build command.
//
// Only functions that are actually implemented in montgomery.cpp are tested
// here: mul128, gcd, mulmod_any, inverse, r_squared_mod_n, reduce256. Add
// sections for the remaining mont:: functions (mulmod, pow2mod_odd,
// pow2mod_any, pow2mod) as they get implemented -- testing them before then
// will fail to link, since montgomery.h declares them but no definition
// exists yet.

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
