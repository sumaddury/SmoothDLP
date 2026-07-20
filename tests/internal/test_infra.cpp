// Standalone correctness tests for the complete, headered functions in
// src/infra.{h,cpp}: buildProductTree, smoothCandidates, treeFactorize,
// ProblemParams, addRelations, crtSolve, and filterDetermined.
//
// infra.{h,cpp} depends on montgomery.h (smoothCandidates calls mont::pow2mod
// /mont::gcd) and, since crtSolve calls into lin_alg::henselLift, on LinBox/
// Givaro/NTL/FFLAS-FFPACK too. This file unity-includes infra.cpp directly
// (not just infra.h) for the same reason infra.cpp itself unity-includes
// lin_alg.cpp: some tests here call lin_alg::rankModPrime directly (to
// characterize crtSolve's failure mode against rank), and infra.h doesn't
// expose lin_alg's types/declarations -- getting them via a normal
// "#include lin_alg.h" here, alongside infra.cpp's own unity-include of it,
// would be a second translation unit pulling in the same LinBox commentator
// machinery, hitting the ~60 duplicate-symbol errors documented in
// lin_alg.cpp's top-of-file comment. Unity-including infra.cpp keeps
// everything in one TU instead. See the Makefile in this directory for the
// exact build command (`make -C tests/internal test`).

#include "infra.cpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
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

// ---- infra::crtSolve -------------------------------------------------------
//
// Contract: crtSolve(params, M, X) solves M*L === X (mod p-1) by Hensel
// lifting each of p-1's prime-power factors independently (lin_alg::
// henselLift) and CRT-recombining (Garner's algorithm) the results into one
// solution mod p-1. As with lin_alg::solveModPrime, rank deficiency is not
// an error -- the returned L must satisfy every row, and must equal the
// planted solution exactly on every coordinate outside whatever kernel a
// factor's solve had.

// Rows of M*C === X (mod p-1) that C fails to satisfy -- the infra-level
// analogue of lin_alg's own unsatisfiedRows test helper, checked against
// the composite modulus p-1 (never itself passed to a field solve, but
// still a well-defined thing to check a candidate solution against).
std::vector<size_t> unsatisfied_rows_mod_p_minus_1(const RelationMatrix& M, const U128Vector& X,
                                                    const U128Vector& C, u128 p_minus_1) {
    std::vector<size_t> bad;
    if (X.size() != M.size()) {
        for (size_t i = 0; i < M.size(); ++i) bad.push_back(i);
        return bad;
    }
    for (size_t i = 0; i < M.size(); ++i) {
        u128 acc = 0;
        bool in_range = true;
        for (const auto& entry : M[i]) {
            if (entry.first >= C.size()) { in_range = false; break; }
            u128 coeff = (u128)entry.second % p_minus_1;
            acc = (acc + (coeff * C[entry.first]) % p_minus_1) % p_minus_1;
        }
        if (!in_range || acc != X[i] % p_minus_1) bad.push_back(i);
    }
    return bad;
}

u128 crt_random_below(u128 q) {
    uint64_t hi = rng(), lo = rng();
    return ((((u128)hi << 64) | (u128)lo) % (q - 1)) + 1;
}

// Builds an (n_cols + extra_rows) x n_cols sparse system over params' own
// factor base, with a planted L_true mod p-1 and X set to the reduced dot
// product -- the shape crtSolve's tests need, since it solves mod the
// composite p-1, not a single prime q the way lin_alg's own buildSystem
// does.
struct CrtSystem {
    RelationMatrix M;
    U128Vector X;
    U128Vector L_true;
};

CrtSystem build_crt_system(const ProblemParams& params, size_t extra_rows, int weight) {
    const size_t n = params.p_levels[0].size();
    const u128 p_minus_1 = params.p - 1;
    const size_t m = n + extra_rows;

    CrtSystem s;
    s.M.assign(m, SparseList{});
    for (size_t i = 0; i < m; ++i) {
        s.M[i].push_back({(i * 7) % n, (uint32_t)(1 + rng() % 1000)});
        for (int k = 0; k < weight; ++k)
            s.M[i].push_back({rng() % n, (uint32_t)(1 + rng() % 1000)});
    }
    s.L_true.resize(n);
    for (size_t j = 0; j < n; ++j) s.L_true[j] = crt_random_below(p_minus_1);
    s.X.assign(m, 0);
    for (size_t i = 0; i < m; ++i) {
        u128 acc = 0;
        for (const auto& e : s.M[i])
            acc = (acc + ((u128)e.second % p_minus_1) * s.L_true[e.first]) % p_minus_1;
        s.X[i] = acc;
    }
    return s;
}

// crtSolve now returns an EMPTY vector when henselLift can't reach full
// precision for some factor (see infra.h's doc) -- a real, data-dependent
// possibility, not a bug. These tests assert success for their specific,
// already-observed-to-work inputs (a regression pin), but check emptiness
// explicitly rather than running row checks against an empty C. The
// dedicated statistics test below is what actually characterizes how often
// empty happens and whether more relations fixes it -- these are not the
// tests for that question.

// crt_solve_failure_vs_rank_by_relation_multiple (below) measures this
// directly: failure rate reaches a stable 0% at m=3k (three times the
// factor-base size) across every prime tested, and stays there through
// 30k -- never creeps back up. These regression tests build their systems
// at 3k+ accordingly, not at an arbitrary absolute row count (an earlier
// version of this file used absolute extra_rows values, which turned out
// to mean wildly different things for a k=6 factor base vs a k=64 one --
// exactly the flaw the k-relative statistics test below was written to
// stop hiding).

void test_crt_solve_recovers_planted_solution_small_prime() {
    // p = 1009: p-1 = 1008 = 2^4 * 3^2 * 7 -- three distinct prime-power
    // factors, two of them with exponent > 1, so this exercises both the
    // multi-factor CRT recombination and multi-level Hensel lifting in one
    // hand-verifiable case.
    ProblemParams params(1009);
    const size_t k = params.p_levels[0].size();
    CrtSystem s = build_crt_system(params, /*extra_rows=*/2 * k, /*weight=*/3);

    U128Vector C = crtSolve(params, s.M, s.X);
    CHECK(!C.empty());
    CHECK(unsatisfied_rows_mod_p_minus_1(s.M, s.X, C, params.p - 1).empty());
    CHECK(C == s.L_true);
}

void test_crt_solve_matches_planted_solution_across_prime_sizes() {
    // A spread across the paper's target range, matching
    // test_problem_params_invariants_across_prime_sizes' choices -- p-1's
    // factorization shape (how many distinct primes, how large the
    // exponents) varies a lot across these, which is exactly what exercises
    // different paths through the CRT loop and Hensel lift depths.
    const u128 primes[] = {
        (u128)1009,
        (u128)8191,                          // p-1 = 8190 = 2*3^2*5*7*13
        (u128)524287,                        // p-1 = 2*3*7*11*31*41
        (u128)2147483647ULL,                 // p-1 = 2*3*7*11*31*151*331
    };
    for (u128 p : primes) {
        ProblemParams params(p);
        const size_t k = params.p_levels[0].size();
        CrtSystem s = build_crt_system(params, /*extra_rows=*/2 * k, /*weight=*/3);

        U128Vector C = crtSolve(params, s.M, s.X);
        CHECK(!C.empty());
        CHECK(unsatisfied_rows_mod_p_minus_1(s.M, s.X, C, params.p - 1).empty());
        CHECK(C == s.L_true);
    }
}

void test_crt_solve_handles_rectangular_overdetermined() {
    ProblemParams params(1009);
    const size_t k = params.p_levels[0].size();
    CrtSystem s = build_crt_system(params, /*extra_rows=*/3 * k, /*weight=*/5);

    U128Vector C = crtSolve(params, s.M, s.X);
    CHECK(!C.empty());
    CHECK(unsatisfied_rows_mod_p_minus_1(s.M, s.X, C, params.p - 1).empty());
    CHECK(C == s.L_true);
}

// End-to-end: real relations from addRelations (not a synthetic planted
// system), solved by crtSolve. There's no independently-known L_true here
// -- the check is that whatever crtSolve returns satisfies every relation
// addRelations actually produced, mod the real p-1, which is the only thing
// that matters once relations come from genuine smooth values rather than a
// hand-planted system.
//
// addRelations itself only ever collects the bare minimum (>= k rows) --
// exactly the m=1k regime the statistics test below shows failing 73-93% of
// the time. Getting a reliable relation set is this test's caller's job,
// same as it would be for any real driver built on top of crtSolve: call
// addRelations repeatedly until comfortably past 3k, not once.
//
// base MUST be a genuine primitive root mod p for this to ever converge --
// addRelations/crtSolve do not check this (nor should they: it's the
// caller's precondition to satisfy, the same way std::lower_bound doesn't
// verify its range is sorted). base=2 is NOT a primitive root mod 1009
// (confirmed by direct order check) and produces a PERMANENT rank ceiling
// no amount of retrying clears -- found by exactly this test, originally
// written with base=2 by mistake. base=11 below is a confirmed primitive
// root mod 1009. k=6 is still tiny enough that a single draw can land
// rank-deficient by ordinary chance even with a real generator (the
// statistics test below averages over many trials to characterize that;
// a single real run doesn't get that luxury), which is what the retry loop
// here is actually for now.
void test_crt_solve_on_real_relations_from_add_relations() {
    u128 p = 1009, base = 11;
    ProblemParams params(p);
    const size_t k = params.p_levels[0].size();
    RelationMatrix M;
    U128Vector X;
    while (M.size() < 3 * k) addRelations(base, params, M, X);

    U128Vector C = crtSolve(params, M, X);
    int retries = 0;
    while (C.empty() && retries < 10) {
        addRelations(base, params, M, X);
        C = crtSolve(params, M, X);
        ++retries;
    }
    std::cout << "  [diag] resolved after " << retries << " retr" << (retries == 1 ? "y" : "ies")
              << ", final M.size()=" << M.size() << "\n";
    if (C.empty()) {
        for (auto& [q, e] : params.p_factorization) {
            if (e < 2) continue;
            size_t rank = lin_alg::rankModPrime(M, k, q);
            std::cout << "  [diag] STILL failing: factor q=" << (uint64_t)q << " e=" << e
                      << " rank=" << rank << "/" << k << "\n";
        }
    }
    CHECK(!C.empty());
    CHECK(unsatisfied_rows_mod_p_minus_1(M, X, C, params.p - 1).empty());
}

// ---- crtSolve: failure-rate characterization -------------------------------
//
// crtSolve returns empty when henselLift can't reach full precision for some
// factor -- a rank-deficient base solve whose kernel direction has raw
// coefficients divisible by q but not q^2 (see lin_alg.h's henselLift doc).
// This is NOT hypothetical: it was found via exactly this kind of testing,
// on real synthetic data derived from an actual ProblemParams(524287). Two
// things need to be true for "return empty, let the caller add relations"
// to be an acceptable fix rather than a silent-failure-shaped one:
//   1. it must be genuinely uncommon, not the common case -- if double
//      index calculus's whole premise is relation sets barely past the
//      minimum k, and this failure hit routinely, "just add relations"
//      wouldn't be a real fix, it'd be a euphemism for "doesn't work".
//   2. it must actually resolve when more relations are added -- if it
//      didn't, "add relations and retry" would be advice that never
//      terminates.
// Both are checked here, empirically, across a spread of primes and factor
// shapes, not assumed.

// Relation count expressed in units of k (the factor-base size), since k
// itself varies enormously across primes (single digits for p=1009, in the
// hundreds for p=2147483647) -- an absolute row count like "extra_rows=30"
// means something completely different for each, which is what made the
// earlier version of this test meaningless to compare across primes.
//
// Also tracks rank at each step, for every factor with e>=2 (only those can
// ever trigger this failure -- an e=1 factor never enters henselLift's
// correction loop at all), to answer directly: does the failure rate hit
// zero exactly when those factors reach full column rank, or does it clear
// up earlier?
void test_crt_solve_failure_vs_rank_by_relation_multiple() {
    const u128 primes[] = {(u128)1009, (u128)524287, (u128)2147483647ULL};
    const int multipliers[] = {1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30};
    const int trials = 15;

    for (u128 p : primes) {
        ProblemParams params(p);
        const size_t k = params.p_levels[0].size();

        std::vector<std::pair<u128, int>> multi_factors;
        for (auto& [q, e] : params.p_factorization)
            if (e >= 2) multi_factors.push_back({q, e});

        std::cout << "\np=" << (uint64_t)p << " k=" << k << " factors with e>=2: ";
        for (auto& [q, e] : multi_factors) std::cout << "(" << (uint64_t)q << "," << e << ") ";
        std::cout << "\n";

        bool reached_zero_failure = false;
        for (int mult : multipliers) {
            const size_t m = (size_t)mult * k;
            const size_t extra_rows = m - k;

            int failures = 0, full_rank_trials = 0;
            std::vector<double> avg_deficiency(multi_factors.size(), 0.0);

            for (int t = 0; t < trials; ++t) {
                CrtSystem s = build_crt_system(params, extra_rows, 3);
                U128Vector C = crtSolve(params, s.M, s.X);
                if (C.empty()) ++failures;

                bool all_full = true;
                for (size_t fi = 0; fi < multi_factors.size(); ++fi) {
                    size_t rank = lin_alg::rankModPrime(s.M, k, multi_factors[fi].first);
                    avg_deficiency[fi] += (double)(k - rank);
                    if (rank != k) all_full = false;
                }
                if (all_full) ++full_rank_trials;
            }

            double failure_rate = 100.0 * failures / trials;
            double full_rank_rate = 100.0 * full_rank_trials / trials;
            std::cout << "  m=" << mult << "k: failure=" << failure_rate
                      << "%  all-e>=2-factors-full-rank=" << full_rank_rate << "%  avg deficiency ";
            for (size_t fi = 0; fi < multi_factors.size(); ++fi)
                std::cout << "(q=" << (uint64_t)multi_factors[fi].first << ": "
                          << (avg_deficiency[fi] / trials) << ") ";
            std::cout << "\n";

            // Once fully clear, a further increase in relations must not
            // bring failure back.
            if (failure_rate == 0.0) reached_zero_failure = true;
            else CHECK(!reached_zero_failure);
        }
        // every prime must reach 0% failure somewhere in this sweep --
        // otherwise "add relations" isn't actually a terminating strategy
        // for it within a reasonable relation count.
        CHECK(reached_zero_failure);
    }
}

// ---- infra::filterDetermined -----------------------------------------------
//
// Contract: filterDetermined(base, params, L) returns exactly the (prime,
// log) pairs from L that verify directly against real modular
// exponentiation -- base^L[i] mod p == params.factor_base[i] -- without
// relying on any rank/kernel analysis of how crtSolve arrived at L. This is
// the paper's Omega_g/Omega_b construction (Algorithm 1's second phase).

// Brute-force base^t mod p for every t in [0, p-2], independent of
// mont::powmod_odd (which filterDetermined itself calls) -- a real oracle,
// not a restatement of the implementation under test. Only practical for
// the small primes used below.
std::unordered_map<u128, u128> brute_force_discrete_logs(u128 base, u128 p) {
    std::unordered_map<u128, u128> table;
    u128 val = 1;
    for (u128 t = 0; t < p - 1; ++t) {
        table[val] = t;
        val = (val * base) % p;
    }
    return table;
}

void test_filter_determined_keeps_only_directly_verified_entries() {
    // base=11 is a confirmed primitive root mod 1009 (see the crtSolve
    // tests above); 3k relations is comfortably past the point this
    // module's own statistics test shows failure rate hits a stable 0%.
    u128 p = 1009, base = 11;
    ProblemParams params(p);
    const size_t k = params.p_levels[0].size();

    RelationMatrix M;
    U128Vector X;
    while (M.size() < 3 * k) addRelations(base, params, M, X);

    U128Vector C = crtSolve(params, M, X);
    CHECK(!C.empty());

    const auto oracle = brute_force_discrete_logs(base, p);
    const auto determined = filterDetermined(base, params, C);

    // At 3k relations for a genuine primitive root, C should already be
    // fully correct, so nothing should be dropped -- but the real property
    // under test is that every pair returned is independently verifiable
    // against the oracle, regardless of how much of C happened to survive.
    CHECK(!determined.empty());
    CHECK(determined.size() <= C.size());
    for (auto& [prime, log] : determined) {
        auto it = oracle.find(prime);
        CHECK(it != oracle.end());
        if (it != oracle.end()) CHECK(it->second == log);
    }
}

void test_filter_determined_excludes_a_corrupted_coordinate_and_keeps_the_rest() {
    // Deterministic, not dependent on ever actually hitting a rank-deficient
    // crtSolve output: take a real, verified-correct C and corrupt exactly
    // one coordinate by hand. filterDetermined must drop that one pair and
    // keep every other -- an inverted or mismatched predicate would get
    // this backwards (drop everything correct, keep the wrong one), so this
    // is a direct regression check on that failure mode.
    u128 p = 1009, base = 11;
    ProblemParams params(p);
    const size_t k = params.p_levels[0].size();

    RelationMatrix M;
    U128Vector X;
    while (M.size() < 3 * k) addRelations(base, params, M, X);

    U128Vector C = crtSolve(params, M, X);
    CHECK(!C.empty());
    CHECK(C.size() >= 2);

    const size_t corrupt_idx = 0;
    U128Vector C_corrupt = C;
    C_corrupt[corrupt_idx] = (C_corrupt[corrupt_idx] + 1) % (p - 1); // now genuinely wrong

    const auto determined = filterDetermined(base, params, C_corrupt);

    const uint32_t corrupted_prime = params.factor_base[corrupt_idx];
    for (auto& [prime, log] : determined) CHECK(prime != corrupted_prime);

    for (size_t i = 1; i < C.size(); ++i) {
        const uint32_t prime = params.factor_base[i];
        bool found = false;
        for (auto& [p2, log] : determined) {
            if (p2 == prime) { CHECK(log == C[i]); found = true; break; }
        }
        CHECK(found);
    }
}

void test_filter_determined_on_empty_input_returns_empty() {
    ProblemParams params(1009);
    U128Vector empty;
    CHECK(filterDetermined((u128)11, params, empty).empty());
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
        {"crt_solve_recovers_planted_solution_small_prime", test_crt_solve_recovers_planted_solution_small_prime},
        {"crt_solve_matches_planted_solution_across_prime_sizes", test_crt_solve_matches_planted_solution_across_prime_sizes},
        {"crt_solve_handles_rectangular_overdetermined", test_crt_solve_handles_rectangular_overdetermined},
        {"crt_solve_on_real_relations_from_add_relations", test_crt_solve_on_real_relations_from_add_relations},
        {"crt_solve_failure_vs_rank_by_relation_multiple", test_crt_solve_failure_vs_rank_by_relation_multiple},
        {"filter_determined_keeps_only_directly_verified_entries", test_filter_determined_keeps_only_directly_verified_entries},
        {"filter_determined_excludes_a_corrupted_coordinate_and_keeps_the_rest", test_filter_determined_excludes_a_corrupted_coordinate_and_keeps_the_rest},
        {"filter_determined_on_empty_input_returns_empty", test_filter_determined_on_empty_input_returns_empty},
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
