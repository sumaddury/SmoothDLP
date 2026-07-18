// Standalone correctness tests for the complete, headered functions in
// src/infra.{h,cpp}: buildProductTree, smoothCandidates, treeFactorize.
//
// (linSolve/crtSolve/rank_relation_gf2 are disabled in infra.cpp -- the
// LinBox/Givaro integration they depend on isn't trusted yet, and they were
// never declared in infra.h to begin with. Nothing here touches them.)
//
// infra.{h,cpp} depends on montgomery.h (smoothCandidates calls mont::pow2mod
// /mont::gcd), so this suite links montgomery.cpp too, but -- now that the
// LinBox-based code is disabled -- needs nothing beyond GMP, same as
// test_montgomery.cpp. See the Makefile in this directory for the exact
// build command (`make -C tests/internal test`).

#include "infra.h"
#include "montgomery.h"
#include "types.h"

#include <gmpxx.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace infra;

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

std::mt19937_64 rng(0x51DEBEEFULL);

u128 random_u128() {
    uint64_t hi = rng(), lo = rng();
    return ((u128)hi << 64) | (u128)lo;
}

// Random value with at most `bits` significant bits (bits in [1, 128]).
u128 random_u128_bits(int bits) {
    u128 v = random_u128();
    if (bits < 128) v &= (((u128)1 << bits) - 1);
    return v;
}

const std::vector<uint32_t> kSmallFactorBase = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};

MpzVector factor_base_mpz() {
    MpzVector v;
    for (uint32_t p : kSmallFactorBase) v.emplace_back(p);
    return v;
}

// ---- infra::buildProductTree -------------------------------------------
//
// Contract: buildProductTree(level) returns the levels of a binary product
// tree bottom-up -- level 0 is the input values themselves; each level above
// pairs up adjacent entries and multiplies them (odd one out carries
// through unchanged); the top level is a single value: the product of
// everything in the input.

void test_build_product_tree_known_case() {
    MpzVector values = {2, 3, 5, 7};
    auto levels = buildProductTree(values);
    CHECK(levels.front() == values);
    CHECK(levels.back().size() == 1);
    CHECK(levels.back()[0] == 210); // 2*3*5*7
}

void test_build_product_tree_odd_count_carries_through() {
    MpzVector values = {2, 3, 5}; // odd count -> last level-0 entry carries up unpaired
    auto levels = buildProductTree(values);
    CHECK(levels.back()[0] == 30); // 2*3*5, regardless of pairing order
}

void test_build_product_tree_single_value() {
    MpzVector values = {97};
    auto levels = buildProductTree(values);
    CHECK(levels.size() == 1);
    CHECK(levels[0][0] == 97);
}

void test_build_product_tree_root_matches_reference_product_random() {
    constexpr int kIters = 200;
    for (int iter = 0; iter < kIters; ++iter) {
        std::size_t n = 1 + (rng() % 40);
        MpzVector values;
        values.reserve(n);
        mpz_class expected = 1;
        for (std::size_t i = 0; i < n; ++i) {
            mpz_class v = u128_to_mpz(random_u128_bits(1 + (int)(rng() % 20)) | 1);
            values.push_back(v);
            expected *= v;
        }
        auto levels = buildProductTree(values);
        CHECK(levels.back()[0] == expected);
    }
}

// ---- infra::smoothCandidates --------------------------------------------
//
// Contract: smoothCandidates(p_levels, X) returns the indices of X whose
// value is fully "smooth" over the factor base p_levels was built from --
// i.e. every prime factor of X[i] appears somewhere in the factor base.
// Internally: reduce the factor-base product mod each candidate, repeatedly
// square (Bernstein's trick, via mont::pow2mod) until the exponent exceeds
// the candidate's bit length, then gcd against the candidate -- a full gcd
// match means the candidate's entire value divides out over the base.

void test_smooth_candidates_known_smooth_and_nonsmooth() {
    auto levels = buildProductTree(factor_base_mpz());
    U128Vector X = {
        (u128)2 * 3 * 5 * 7 * 11,   // smooth
        (u128)2 * 3 * 53,           // 53 outside the base -> not smooth
        (u128)41 * 43 * 47,         // smooth
        (u128)59 * 61,              // both outside the base -> not smooth
    };
    auto idx = smoothCandidates(levels, X);
    std::sort(idx.begin(), idx.end());
    std::vector<size_t> expected = {0, 2};
    CHECK(idx == expected);
}

void test_smooth_candidates_randomized_constructed_batch() {
    auto levels = buildProductTree(factor_base_mpz());
    const std::vector<uint32_t> outside_base = {53, 59, 61, 67, 71, 73, 79};

    U128Vector X;
    std::vector<bool> expected_smooth;
    constexpr int kBatch = 100;
    for (int i = 0; i < kBatch; ++i) {
        bool make_smooth = (rng() % 2) == 0;
        u128 x = 1;
        int factors = 1 + (int)(rng() % 4);
        for (int f = 0; f < factors; ++f) x *= kSmallFactorBase[rng() % kSmallFactorBase.size()];
        if (!make_smooth) x *= outside_base[rng() % outside_base.size()];
        X.push_back(x);
        expected_smooth.push_back(make_smooth);
    }

    auto idx = smoothCandidates(levels, X);
    std::vector<bool> got_smooth(X.size(), false);
    for (size_t i : idx) got_smooth[i] = true;

    for (size_t i = 0; i < X.size(); ++i) {
        CHECK(got_smooth[i] == expected_smooth[i]);
    }
}

void test_smooth_candidates_single_prime_candidates() {
    // Every factor-base prime alone should be reported smooth against its
    // own base (a direct exercise of the pow2mod call site for many
    // distinct, small odd/even moduli in one batch).
    auto levels = buildProductTree(factor_base_mpz());
    U128Vector X;
    for (uint32_t p : kSmallFactorBase) X.push_back(p);
    auto idx = smoothCandidates(levels, X);
    CHECK(idx.size() == X.size());
}

// ---- infra::treeFactorize ------------------------------------------------
//
// Contract: treeFactorize(p_levels, d) fully factors a known-smooth `d` over
// the product tree's base primes, returning a sparse list of
// (factor-base index, exponent) pairs whose reconstruction equals d exactly.

u128 reconstruct(const SparseList& factors) {
    u128 acc = 1;
    for (auto& [idx, exp] : factors) {
        for (uint32_t e = 0; e < exp; ++e) acc *= kSmallFactorBase[idx];
    }
    return acc;
}

void test_tree_factorize_known_composite() {
    auto levels = buildProductTree(factor_base_mpz());
    u128 d = (u128)8 * 3 * 11; // 2^3 * 3 * 11
    auto factors = treeFactorize(levels, d);
    CHECK(reconstruct(factors) == d);
}

void test_tree_factorize_randomized_reconstruction() {
    auto levels = buildProductTree(factor_base_mpz());
    constexpr int kIters = 200;
    for (int iter = 0; iter < kIters; ++iter) {
        u128 d = 1;
        int factors = 1 + (int)(rng() % 6);
        for (int f = 0; f < factors; ++f) d *= kSmallFactorBase[rng() % kSmallFactorBase.size()];

        auto result = treeFactorize(levels, d);
        u128 got = reconstruct(result);
        if (got != d) {
            report_failure(__FILE__, __LINE__,
                "treeFactorize reconstruction mismatch for d=" +
                u128_to_mpz(d).get_str());
        }
        ++g_checks;
    }
}

// ---- infra::addRelations -------------------------------------------------
//
// Contract: addRelations(base, params, M, X) grows M/X by relation batches
// (random t in [0, p-2], base^t mod p, kept iff B-smooth over
// params.p_levels) until at least k = params.p_levels[0].size() new rows
// have been added. Each new row of M is that relation's sparse
// factorization; the matching entry appended to X is the exponent t that
// produced it -- NOT the smooth value itself -- since X is the right-hand
// side of the M*L === X (mod p-1) system the rest of the pipeline solves
// for base's discrete logs. That's the specific thing these tests verify:
// every (M[i], X[i]) pair must satisfy reconstruct(M[i]) == base^X[i] mod p.

ProblemParams make_params(u128 p, u128 mask, int mask_bitlength, u128 B,
                           double smooth_density, const std::vector<MpzVector>& levels) {
    return ProblemParams{p, mont::inverse(p), mont::r_squared_mod_n(p), mask, B,
                          smooth_density, levels, mask_bitlength};
}

void test_add_relations_reaches_factor_base_size_and_reconstructs() {
    // p = 47: every residue mod p is < 47, so it's automatically fully
    // smooth over a base containing every prime <= 47 -- a single batch of
    // exactly k relations always suffices, deterministically (no flakiness
    // from addRelations' internal random_device-seeded draws).
    auto levels = buildProductTree(factor_base_mpz());
    const size_t k = kSmallFactorBase.size();
    u128 p = 47, base = 5;
    ProblemParams params = make_params(p, 63, 6, 47, 1.0, levels);

    std::vector<SparseList> M;
    U128Vector X;
    addRelations(base, params, M, X);

    CHECK(M.size() == k);
    CHECK(X.size() == k);
    for (size_t i = 0; i < M.size(); ++i) {
        u128 want = mont::powmod_odd(base, X[i], params.p, params.p_prime, params.r2);
        CHECK(reconstruct(M[i]) == want);
    }
}

void test_add_relations_partial_density_still_reconstructs() {
    // A realistic case where not every candidate is smooth (p is far larger
    // than the factor base's top prime), and smooth_density is deliberately
    // set well above the true rate, so batches under-deliver and multiple
    // rounds are the likely path -- exactly the loop iteration this routine
    // needs to get right. Regardless of how many rounds it actually takes,
    // every relation produced must still be internally consistent.
    auto levels = buildProductTree(factor_base_mpz());
    const size_t k = kSmallFactorBase.size();
    u128 p = 1009, base = 2;
    ProblemParams params = make_params(p, 1023, 10, 47, 0.05, levels);

    std::vector<SparseList> M;
    U128Vector X;
    addRelations(base, params, M, X);

    CHECK(M.size() == X.size());
    CHECK(M.size() >= k);
    for (size_t i = 0; i < M.size(); ++i) {
        CHECK(X[i] <= params.p - 2);
        u128 want = mont::powmod_odd(base, X[i], params.p, params.p_prime, params.r2);
        CHECK(reconstruct(M[i]) == want);
    }
}

struct TestCase {
    const char* name;
    void (*fn)();
};

} // namespace

int main() {
    TestCase tests[] = {
        {"build_product_tree_known_case", test_build_product_tree_known_case},
        {"build_product_tree_odd_count_carries_through", test_build_product_tree_odd_count_carries_through},
        {"build_product_tree_single_value", test_build_product_tree_single_value},
        {"build_product_tree_root_matches_reference_product_random", test_build_product_tree_root_matches_reference_product_random},
        {"smooth_candidates_known_smooth_and_nonsmooth", test_smooth_candidates_known_smooth_and_nonsmooth},
        {"smooth_candidates_randomized_constructed_batch", test_smooth_candidates_randomized_constructed_batch},
        {"smooth_candidates_single_prime_candidates", test_smooth_candidates_single_prime_candidates},
        {"tree_factorize_known_composite", test_tree_factorize_known_composite},
        {"tree_factorize_randomized_reconstruction", test_tree_factorize_randomized_reconstruction},
        {"add_relations_reaches_factor_base_size_and_reconstructs", test_add_relations_reaches_factor_base_size_and_reconstructs},
        {"add_relations_partial_density_still_reconstructs", test_add_relations_partial_density_still_reconstructs},
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
