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
#include "gauss_dream.h"
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

// ---- infra::ProblemParams -------------------------------------------------
//
// Contract: ProblemParams(p) self-computes everything addRelations needs for
// that p: Montgomery constants for p; mask_bitlength/mask, a rejection-
// sampling range covering exactly [0, p-2]; a smoothness bound B; the
// product tree built over the factor base up to B; p's own factorization;
// and smooth_density, which must be an actual probability in (0, 1] (not
// logDickman's raw log-space output, which is <= 0 and would make
// addRelations divide by a non-positive number).

u128 find_prime_at_least(u128 start) {
    u128 candidate = start | 1;
    for (int i = 0; i < 1'000'000; ++i) {
        if (gauss::isPrime(candidate)) return candidate;
        candidate += 2;
    }
    throw std::runtime_error("find_prime_at_least: no prime found in range");
}

void check_problem_params_invariants(u128 p) {
    ProblemParams params(p);

    // mask_bitlength/mask describe exactly the range [0, p-2]: mask is the
    // smallest all-ones bitmask that covers every set bit of p-2.
    CHECK(params.mask_bitlength == 128 - clz128(p - 2));
    CHECK((params.mask & (p - 2)) == (p - 2));
    CHECK(params.mask_bitlength == 128 || (params.mask >> params.mask_bitlength) == 0);

    // B is a usable smoothness bound, and every factor-base prime is <= B.
    CHECK(params.B >= 2);
    for (const mpz_class& prime : params.p_levels[0]) {
        CHECK(mpz_to_u128(prime) <= params.B);
    }

    // The product tree's root is exactly the product of the base primes.
    mpz_class expected_root = 1;
    for (const mpz_class& prime : params.p_levels[0]) expected_root *= prime;
    CHECK(params.p_levels.back()[0] == expected_root);

    // p_factorization is (p-1)'s factorization (needed downstream for
    // CRT/Hensel lifting over the group order), not p's -- it must
    // reconstruct to p-1 exactly.
    u128 reconstructed = 1;
    for (auto& [prime, exp] : params.p_factorization)
        for (uint32_t e = 0; e < exp; ++e) reconstructed *= prime;
    CHECK(reconstructed == p - 1);

    CHECK(params.smooth_density > 0.0);
    CHECK(params.smooth_density <= 1.0);
}

void test_problem_params_invariants_across_prime_sizes() {
    // A spread across the paper's actual target range (30-75 bits), plus a
    // couple of small primes below it -- see CLAUDE.md's scope note that
    // this project targets 30-75 bit primes, not the full u128 range.
    const u128 primes[] = {
        (u128)47,
        (u128)1009,
        (u128)8191,                        // 2^13 - 1, Mersenne prime
        (u128)524287,                       // 2^19 - 1, Mersenne prime
        (u128)2147483647ULL,                 // 2^31 - 1, Mersenne prime
        (u128)2305843009213693951ULL,        // 2^61 - 1, Mersenne prime
        find_prime_at_least((u128)1 << 74),  // ~75-bit, the paper's top end
    };
    for (u128 p : primes) check_problem_params_invariants(p);
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
// every (M[i], X[i]) pair must satisfy reconstruct(M[i]) == base^X[i] mod p,
// against ProblemParams's own self-computed factor base (not a hand-picked
// one), since ProblemParams no longer accepts one from the caller.

u128 reconstruct_from_params(const ProblemParams& params, const SparseList& factors) {
    u128 acc = 1;
    for (auto& [idx, exp] : factors) {
        u128 prime = mpz_to_u128(params.p_levels[0][idx]);
        for (uint32_t e = 0; e < exp; ++e) acc *= prime;
    }
    return acc;
}

void test_add_relations_reaches_factor_base_size_and_reconstructs() {
    // p = 47: ProblemParams's own heuristic picks a tiny factor base for a
    // prime this small, so a handful of relations always suffices quickly.
    u128 p = 47, base = 5;
    ProblemParams params(p);
    const size_t k = params.p_levels[0].size();

    RelationMatrix M;
    U128Vector X;
    addRelations(base, params, M, X);

    CHECK(M.size() >= k);
    CHECK(X.size() == M.size());
    for (size_t i = 0; i < M.size(); ++i) {
        u128 want = mont::powmod_odd(base, X[i], params.p, params.p_prime, params.r2);
        CHECK(reconstruct_from_params(params, M[i]) == want);
    }
}

void test_add_relations_larger_prime_still_reconstructs() {
    // A realistic case where not every candidate is smooth (p is far larger
    // than the factor base's top prime), so batches may need multiple
    // rounds -- exactly the loop iteration this routine needs to get right.
    // Regardless of how many rounds it actually takes, every relation
    // produced must still be internally consistent.
    u128 p = 1009, base = 2;
    ProblemParams params(p);
    const size_t k = params.p_levels[0].size();

    RelationMatrix M;
    U128Vector X;
    addRelations(base, params, M, X);

    CHECK(M.size() == X.size());
    CHECK(M.size() >= k);
    for (size_t i = 0; i < M.size(); ++i) {
        CHECK(X[i] <= params.p - 2);
        u128 want = mont::powmod_odd(base, X[i], params.p, params.p_prime, params.r2);
        CHECK(reconstruct_from_params(params, M[i]) == want);
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
        {"problem_params_invariants_across_prime_sizes", test_problem_params_invariants_across_prime_sizes},
        {"add_relations_reaches_factor_base_size_and_reconstructs", test_add_relations_reaches_factor_base_size_and_reconstructs},
        {"add_relations_larger_prime_still_reconstructs", test_add_relations_larger_prime_still_reconstructs},
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
