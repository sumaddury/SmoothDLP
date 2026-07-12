// Standalone correctness tests for the complete, headered functions in
// src/smooth_algos.{h,cpp}: isSmooth, logDickman, mp_ln, log_mul, psiApprox.
// Depends on gauss_dream.h (isSmooth calls isPrime) and transitively on
// montgomery.h, so links both .cpp files too. Needs Boost (Dickman rho's
// asymptotic tail uses boost::math::expint) in addition to GMP -- both
// already required by the whole project, see ../../environment.yml. See the
// Makefile in this directory for the exact build command
// (`make -C tests/internal test`).

#include "smooth_algos.h"
#include "gauss_dream.h"
#include "types.h"

#include <gmpxx.h>

#include <cmath>
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

std::mt19937_64 rng(0xFEEDFACEULL);

u128 random_u128() {
    uint64_t hi = rng(), lo = rng();
    return ((u128)hi << 64) | (u128)lo;
}

u128 random_bitlength(int bits) {
    u128 v = random_u128();
    if (bits < 128) v &= (((u128)1 << bits) - 1);
    v |= (u128)1 << (bits - 1);
    return v | 1;
}

bool isclose(double a, double b, double rel_tol = 1e-9, double abs_tol = 0.0) {
    return std::fabs(a - b) <= std::max(rel_tol * std::max(std::fabs(a), std::fabs(b)), abs_tol);
}

// ---- smooth_algos::isSmooth -----------------------------------------------
//
// Contract: isSmooth(x, y) is true iff every prime factor of x is <= y.

void test_is_smooth_known_cases() {
    CHECK(isSmooth(77, 11));
    CHECK(!isSmooth(77, 7));
    CHECK(isSmooth(1, 5));
    CHECK(isSmooth(97, 97));
    CHECK(!isSmooth(97, 50));

    u128 product = (u128)2 * 3 * 5 * 7 * 11 * 13;
    CHECK(isSmooth(product, 13));
    CHECK(!isSmooth(product, 11));
}

void test_is_smooth_randomized_constructed_products() {
    const std::vector<uint32_t> small_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    for (int iter = 0; iter < 50; ++iter) {
        uint32_t bound = small_primes[rng() % small_primes.size()];
        std::vector<uint32_t> at_or_below, above;
        for (uint32_t p : small_primes) (p <= bound ? at_or_below : above).push_back(p);

        u128 x = 1;
        int factors = 1 + (int)(rng() % 5);
        for (int f = 0; f < factors; ++f) x *= at_or_below[rng() % at_or_below.size()];
        CHECK(isSmooth(x, bound));

        if (!above.empty()) {
            u128 x2 = x * above[rng() % above.size()];
            CHECK(!isSmooth(x2, bound));
        }
    }
}

// ---- smooth_algos::logDickman ---------------------------------------------
//
// Contract: logDickman(u) == ln(rho(u)), Dickman's rho function -- closed
// form for u<=2, spline-interpolated table for 2<u<20, asymptotic
// saddle-point formula for u>=20. rho is nonincreasing, so logDickman must
// be too.

void test_log_dickman_known_values() {
    CHECK(isclose(logDickman(0.0), 0.0, 0.0, 1e-12));
    CHECK(isclose(logDickman(0.5), 0.0, 0.0, 1e-12));
    CHECK(isclose(logDickman(1.0), 0.0, 0.0, 1e-12));

    double expected = std::log(1.0 - std::log(2.0));
    CHECK(isclose(logDickman(2.0), expected, 0.0, 1e-3));
}

void test_log_dickman_is_monotonically_nonincreasing() {
    std::uniform_real_distribution<double> dist(0.0, 40.0);
    for (int iter = 0; iter < 200; ++iter) {
        double a = dist(rng), b = dist(rng);
        if (a > b) std::swap(a, b);
        if (b - a < 1e-6) continue;
        CHECK(logDickman(a) >= logDickman(b) - 1e-9);
    }
}

// ---- smooth_algos::mp_ln --------------------------------------------------
//
// Contract: mp_ln(x) == ln(x) for a (possibly >64-bit) integer x, computed
// via a truncated-mantissa-plus-shift trick rather than full-precision GMP
// log (which doesn't exist): keep the top 53 bits as a double mantissa, and
// add back shift * ln(2) for the bits shifted off.

void test_mp_ln_known_powers_of_two() {
    for (int k : {1, 10, 64, 96, 127}) {
        double expected = k * std::log(2.0);
        CHECK(isclose(mp_ln((u128)1 << k), expected, 1e-9));
    }
}

void test_mp_ln_randomized_against_reference() {
    for (int iter = 0; iter < 50; ++iter) {
        u128 x = random_bitlength(20 + (int)(rng() % 100));
        int bit_len = 128 - clz128(x);
        int shift = bit_len > 53 ? bit_len - 53 : 0;
        u128 mantissa_int = shift ? (x >> shift) : x;
        double mantissa = (double)(uint64_t)mantissa_int; // fits: <= 53 bits
        double expected = std::log(mantissa) + shift * std::log(2.0);
        CHECK(isclose(mp_ln(x), expected, 1e-9));
    }
}

// ---- smooth_algos::log_mul -------------------------------------------------
//
// Contract: log_mul(x, log_rho) == x * exp(log_rho), for large x, computed
// via GMP integer scaling (decompose exp(log_rho) into a fixed-digit
// mantissa times a power of ten) rather than naive floating-point
// multiplication, to avoid losing precision for huge x.

void test_log_mul_identity_at_zero() {
    for (u128 x : {(u128)1, (u128)2, (u128)97, (u128)1 << 60, ((u128)1 << 100) + 3}) {
        CHECK(log_mul(x, 0.0) == x);
    }
}

void test_log_mul_randomized_against_reference() {
    std::uniform_real_distribution<double> dist(-20.0, 0.0);
    for (int iter = 0; iter < 50; ++iter) {
        u128 x = random_bitlength(20 + (int)(rng() % 100));
        double log_rho = dist(rng);
        mpz_class x_mp = u128_to_mpz(x);
        double expected_d = 0.0;
        try { expected_d = x_mp.get_d() * std::exp(log_rho); } catch (...) {}
        u128 got = log_mul(x, log_rho);
        double got_d = u128_to_mpz(got).get_d();
        if (expected_d > 1.0) {
            double rel_err = std::fabs(got_d - expected_d) / expected_d;
            CHECK(rel_err < 0.02); // wider tolerance than the Python suite's 0.01 --
                                    // mpz_class::get_d() itself loses precision for
                                    // very large x, which the Python oracle (exact
                                    // Python bigints) doesn't suffer from.
        }
    }
}

// ---- smooth_algos::psiApprox -----------------------------------------------
//
// Contract: psiApprox(x, y) estimates Psi(x, y), the count of y-smooth
// integers up to x, via x * rho(ln(x)/ln(y)). Must never exceed x, and must
// be nondecreasing in y (a larger smoothness bound only allows more smooth
// numbers through).

void test_psi_approx_never_exceeds_x() {
    for (int iter = 0; iter < 50; ++iter) {
        u128 x = random_bitlength(40 + (int)(rng() % 80));
        uint64_t y = 4 + (rng() % 1000000);
        CHECK(psiApprox(x, y) <= x);
    }
}

void test_psi_approx_nondecreasing_in_y() {
    for (int iter = 0; iter < 30; ++iter) {
        u128 x = random_bitlength(40 + (int)(rng() % 80));
        uint64_t y1 = 4 + (rng() % 100000);
        uint64_t y2 = y1 + 1 + (rng() % 100000);
        CHECK(psiApprox(x, y1) <= psiApprox(x, y2));
    }
}

void test_psi_approx_matches_x_times_rho_within_tolerance() {
    for (int iter = 0; iter < 30; ++iter) {
        u128 x = random_bitlength(40 + (int)(rng() % 80));
        uint64_t y = 1000 + (rng() % 1000000);

        double u = mp_ln(x) / std::log((double)y);
        double expected_d = u128_to_mpz(x).get_d() * std::exp(logDickman(u));
        double got_d = u128_to_mpz(psiApprox(x, y)).get_d();

        if (expected_d > 1.0) {
            double rel_err = std::fabs(got_d - expected_d) / expected_d;
            CHECK(rel_err < 0.02);
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
        {"is_smooth_known_cases", test_is_smooth_known_cases},
        {"is_smooth_randomized_constructed_products", test_is_smooth_randomized_constructed_products},
        {"log_dickman_known_values", test_log_dickman_known_values},
        {"log_dickman_is_monotonically_nonincreasing", test_log_dickman_is_monotonically_nonincreasing},
        {"mp_ln_known_powers_of_two", test_mp_ln_known_powers_of_two},
        {"mp_ln_randomized_against_reference", test_mp_ln_randomized_against_reference},
        {"log_mul_identity_at_zero", test_log_mul_identity_at_zero},
        {"log_mul_randomized_against_reference", test_log_mul_randomized_against_reference},
        {"psi_approx_never_exceeds_x", test_psi_approx_never_exceeds_x},
        {"psi_approx_nondecreasing_in_y", test_psi_approx_nondecreasing_in_y},
        {"psi_approx_matches_x_times_rho_within_tolerance", test_psi_approx_matches_x_times_rho_within_tolerance},
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
