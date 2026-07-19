// Standalone correctness tests for src/lin_alg.{h,cpp}: solveModPrime and
// rankModPrime, the LinBox-backed sparse solve of M*L === X (mod q).
//
// Unlike the other suites in this directory, this one DOES link LinBox,
// Givaro, NTL and FFLAS-FFPACK -- it is the only part of the project that
// touches them. See the Makefile in this directory for the exact build
// command (`make -C tests/internal test`).
//
// Systems are built by planting a known solution L_true, forming
// X = M*L_true (mod q), and then checking what comes back. Two distinct
// properties are checked throughout, and they are NOT the same thing:
//   * the returned L satisfies every row (result.unsatisfied_rows empty), and
//   * the returned L equals L_true coordinate-for-coordinate.
// The first always holds for a consistent system; the second holds only for
// coordinates that are uniquely determined, which is exactly the
// partial-rank behaviour the double-index-calculus driver depends on.

#include "lin_alg.h"
#include "types.h"

#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace lin_alg;

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

std::mt19937_64 rng(0xC0FFEEULL);

// A few primes spanning the supported range. 2^32 matters because the
// field this module deliberately avoids (Modular<uint64_t,__uint128_t>)
// starts silently returning wrong answers just above it -- everything from
// P_2_33 up is a regression pin against that.
const u128 P_SMALL  = 1000003ULL;
const u128 P_2_31   = 2147483647ULL;
const u128 P_2_33   = 8589934583ULL;
const u128 P_2_40   = 1099511627791ULL;
const u128 P_2_63   = 9223372036854775783ULL;
const u128 P_2_64   = 18446744073709551557ULL;   // largest prime < 2^64

u128 randomBelow(u128 q) {
    u128 hi = rng(), lo = rng();
    return (((hi << 64) | lo) % (q - 1)) + 1;
}

/**
  Builds an m x n sparse system with a planted solution: each row gets
  `weight` random entries plus one on a rotating diagonal (so no row is
  empty), L_true is random mod q, and X is set to M*L_true mod q.
*/
struct System {
    RelationMatrix M;
    U128Vector X;
    U128Vector L_true;
};

System buildSystem(size_t m, size_t n, int weight, u128 q) {
    System s;
    s.M.assign(m, SparseList{});
    for (size_t i = 0; i < m; ++i) {
        s.M[i].push_back({(i * 7) % n, (uint32_t)(1 + rng() % 1000)});
        for (int k = 0; k < weight; ++k)
            s.M[i].push_back({rng() % n, (uint32_t)(1 + rng() % 1000)});
    }
    s.L_true.resize(n);
    for (size_t j = 0; j < n; ++j) s.L_true[j] = randomBelow(q);
    s.X.assign(m, 0);
    for (size_t i = 0; i < m; ++i) {
        u128 acc = 0;
        for (const auto& e : s.M[i])
            acc = (acc + ((u128)e.second % q) * s.L_true[e.first]) % q;
        s.X[i] = acc;
    }
    return s;
}

bool matchesPlanted(const SolveResult& r, const U128Vector& L_true) {
    if (r.L.size() != L_true.size()) return false;
    for (size_t j = 0; j < L_true.size(); ++j)
        if (r.L[j] != L_true[j]) return false;
    return true;
}

// ---- solveModPrime: basic correctness ---------------------------------------
//
// Contract: for a consistent, full-rank system the returned L both satisfies
// every row and equals the planted solution exactly.

void test_solve_recovers_planted_solution_small_prime() {
    System s = buildSystem(40, 20, 3, P_SMALL);
    SolveResult r = solveModPrime(s.M, s.X, 20, P_SMALL);
    CHECK(r.status == SolveStatus::OK);
    CHECK(r.unsatisfied_rows.empty());
    CHECK(matchesPlanted(r, s.L_true));
}

// The whole point of the field choice in lin_alg.cpp. Every one of these
// primes above 2^32 is a size at which the naive Modular<uint64_t,__uint128_t>
// field silently produces a wrong vector; all of them must round-trip here.
void test_solve_is_correct_across_prime_sizes() {
    const u128 primes[] = {P_SMALL, P_2_31, P_2_33, P_2_40, P_2_63, P_2_64};
    for (u128 q : primes) {
        System s = buildSystem(60, 30, 3, q);
        SolveResult r = solveModPrime(s.M, s.X, 30, q);
        CHECK(r.status == SolveStatus::OK);
        CHECK(r.unsatisfied_rows.empty());
        CHECK(matchesPlanted(r, s.L_true));
    }
}

void test_solve_randomized_full_rank_systems() {
    for (int iter = 0; iter < 25; ++iter) {
        size_t n = 10 + rng() % 40;
        size_t m = n + 10 + rng() % 30;          // overdetermined
        System s = buildSystem(m, n, 3, P_2_64);
        SolveResult r = solveModPrime(s.M, s.X, n, P_2_64);
        CHECK(r.status == SolveStatus::OK);
        CHECK(r.unsatisfied_rows.empty());
    }
}

// Relation collection produces far more rows than columns; that shape must
// work, and with enough independent rows the solution is pinned exactly.
void test_solve_handles_rectangular_overdetermined() {
    System s = buildSystem(300, 100, 4, P_2_64);
    SolveResult r = solveModPrime(s.M, s.X, 100, P_2_64);
    CHECK(r.status == SolveStatus::OK);
    CHECK(r.unsatisfied_rows.empty());
    CHECK(matchesPlanted(r, s.L_true));
}

// X entries arriving unreduced (they are exponents mod p-1, not mod q)
// must give the same answer as pre-reduced ones.
void test_solve_reduces_unreduced_rhs() {
    System s = buildSystem(40, 20, 3, P_2_40);
    SolveResult a = solveModPrime(s.M, s.X, 20, P_2_40);
    for (size_t i = 0; i < s.X.size(); ++i) s.X[i] += P_2_40 * (1 + (i % 3));
    SolveResult b = solveModPrime(s.M, s.X, 20, P_2_40);
    CHECK(a.status == SolveStatus::OK);
    CHECK(b.status == SolveStatus::OK);
    CHECK(a.L == b.L);
}

// Repeated column indices in one row must be summed, not overwritten.
// LinBox's setEntry overwrites, so without explicit coalescing the matrix
// LinBox holds silently differs from the one the caller described and the
// system reads as unsolvable. buildSystem generates repeats naturally (its
// random columns collide), so several suites above depend on this too --
// this test pins it directly rather than leaving it implicit.
void test_duplicate_column_entries_are_summed() {
    const u128 q = P_2_64;
    RelationMatrix M(3);
    M[0] = {{0, 5}, {1, 7}, {0, 4}};       // column 0 appears twice: 5 + 4 == 9
    M[1] = {{0, 2}, {1, 1}};
    M[2] = {{1, 3}, {0, 6}};
    U128Vector L_true = {11, 23};
    U128Vector X(3);
    for (size_t i = 0; i < 3; ++i) {
        u128 acc = 0;
        for (const auto& e : M[i]) acc = (acc + ((u128)e.second % q) * L_true[e.first]) % q;
        X[i] = acc;
    }
    SolveResult r = solveModPrime(M, X, 2, q);
    CHECK(r.status == SolveStatus::OK);
    CHECK(r.unsatisfied_rows.empty());
    CHECK(r.L.size() == 2 && r.L[0] == 11 && r.L[1] == 23);

    // An explicitly pre-summed row must give the identical answer.
    RelationMatrix M2 = M;
    M2[0] = {{0, 9}, {1, 7}};
    SolveResult r2 = solveModPrime(M2, X, 2, q);
    CHECK(r2.status == SolveStatus::OK);
    CHECK(r2.L == r.L);
}

// ---- solveModPrime: partial rank --------------------------------------------
//
// Contract: rank deficiency is not an error. With a controlled kernel, the
// returned L must still satisfy every row, and must be exactly right on
// every coordinate outside the kernel -- only kernel coordinates may differ.

/**
  Builds a system whose column `dup` is an exact copy of column `src`, so
  the kernel is exactly span(e_src - e_dup): coordinates src and dup are
  the only ones not uniquely determined.
*/
System buildDuplicateColumnSystem(size_t m, size_t n, size_t src, size_t dup, u128 q) {
    System s;
    s.M.assign(m, SparseList{});
    for (size_t i = 0; i < m; ++i) {
        uint32_t shared = (uint32_t)(1 + rng() % 1000);
        s.M[i].push_back({src, shared});
        s.M[i].push_back({dup, shared});               // identical column
        for (int k = 0; k < 3; ++k) {
            size_t j = rng() % n;
            if (j != src && j != dup) s.M[i].push_back({j, (uint32_t)(1 + rng() % 1000)});
        }
    }
    s.L_true.resize(n);
    for (size_t j = 0; j < n; ++j) s.L_true[j] = randomBelow(q);
    s.X.assign(m, 0);
    for (size_t i = 0; i < m; ++i) {
        u128 acc = 0;
        for (const auto& e : s.M[i])
            acc = (acc + ((u128)e.second % q) * s.L_true[e.first]) % q;
        s.X[i] = acc;
    }
    return s;
}

void test_solve_partial_rank_locks_determined_coordinates() {
    const size_t n = 20, src = 3, dup = 17;
    System s = buildDuplicateColumnSystem(60, n, src, dup, P_2_64);

    CHECK(rankModPrime(s.M, n, P_2_64) == n - 1);

    SolveResult r = solveModPrime(s.M, s.X, n, P_2_64);
    CHECK(r.status == SolveStatus::OK);
    CHECK(r.unsatisfied_rows.empty());       // still a genuine solution

    // Every coordinate outside the kernel must be exactly the planted value;
    // only src/dup are permitted to differ.
    for (size_t j = 0; j < n; ++j) {
        if (j == src || j == dup) continue;
        CHECK(r.L[j] == s.L_true[j]);
    }
    // ...and the kernel direction must be respected: L[src]+L[dup] is the
    // invariant, since the two columns are identical.
    u128 planted_sum = (s.L_true[src] + s.L_true[dup]) % P_2_64;
    u128 got_sum     = (r.L[src] + r.L[dup]) % P_2_64;
    CHECK(planted_sum == got_sum);
}

// ---- rankModPrime -----------------------------------------------------------
//
// Contract: rank is a property of M mod q, not of M. The same integer matrix
// can be full rank modulo one prime and deficient modulo another -- which is
// why rank must be recomputed per CRT component and never reused.

void test_rank_is_full_for_generic_system() {
    System s = buildSystem(60, 25, 4, P_2_64);
    CHECK(rankModPrime(s.M, 25, P_2_64) == 25);
}

void test_rank_depends_on_the_modulus() {
    // Column `dup` equals column `src` plus a multiple of q1, so the two
    // columns coincide mod q1 but not mod q2. Entries stay inside uint32_t.
    const u128 q1 = 65521, q2 = 1000003;
    const size_t m = 60, n = 20, src = 3, dup = 11;
    System s;
    s.M.assign(m, SparseList{});
    for (size_t i = 0; i < m; ++i) {
        uint32_t base = (uint32_t)(1 + rng() % 50000);
        s.M[i].push_back({src, base});
        s.M[i].push_back({dup, (uint32_t)(base + (uint32_t)q1 * (1 + rng() % 5))});
        for (int k = 0; k < 3; ++k) {
            size_t j = rng() % n;
            if (j != src && j != dup) s.M[i].push_back({j, (uint32_t)(1 + rng() % 50000)});
        }
    }
    CHECK(rankModPrime(s.M, n, q1) == n - 1);   // columns collide mod q1
    CHECK(rankModPrime(s.M, n, q2) == n);       // and separate mod q2
}

// ---- Wiedemann --------------------------------------------------------------
//
// Contract: on a square system the Wiedemann path must produce a genuine
// solution too (it is the method of choice for large systems, where
// elimination's fill-in dominates).

void test_wiedemann_solves_square_systems() {
    System s = buildSystem(60, 60, 4, P_2_64);
    SolveResult r = solveModPrime(s.M, s.X, 60, P_2_64, Method::Wiedemann);
    if (r.status == SolveStatus::OK) {
        CHECK(r.unsatisfied_rows.empty());
    } else {
        // Wiedemann is randomized and may report failure on an unlucky
        // draw; what it must never do is return OK with a bad vector.
        CHECK(r.status == SolveStatus::FAILED || r.status == SolveStatus::INCONSISTENT);
    }
}

void test_wiedemann_and_elimination_agree_on_row_satisfaction() {
    System s = buildSystem(50, 50, 4, P_2_40);
    SolveResult a = solveModPrime(s.M, s.X, 50, P_2_40, Method::SparseElimination);
    SolveResult b = solveModPrime(s.M, s.X, 50, P_2_40, Method::Wiedemann);
    CHECK(a.status == SolveStatus::OK);
    CHECK(a.unsatisfied_rows.empty());
    if (b.status == SolveStatus::OK) CHECK(b.unsatisfied_rows.empty());
}

// ---- Inconsistent systems ---------------------------------------------------
//
// Contract: an unsolvable system must never come back as OK with an empty
// residual. Either it is flagged INCONSISTENT, or it is OK with the
// offending rows listed -- both are honest answers, a silent false OK is not.

void test_inconsistent_system_is_not_silently_accepted() {
    System s = buildSystem(30, 15, 3, P_2_64);
    s.M.push_back(s.M[0]);                       // duplicate a row...
    s.X.push_back((s.X[0] + 1) % P_2_64);        // ...with a different RHS
    SolveResult r = solveModPrime(s.M, s.X, 15, P_2_64);
    if (r.status == SolveStatus::OK) {
        CHECK(!r.unsatisfied_rows.empty());
    } else {
        CHECK(r.status == SolveStatus::INCONSISTENT || r.status == SolveStatus::FAILED);
    }
}

// ---- Input validation -------------------------------------------------------
//
// Contract: malformed input is rejected before any LinBox call, with a
// specific status rather than a crash.

void test_rejects_modulus_above_supported_range() {
    System s = buildSystem(10, 5, 2, P_SMALL);
    u128 too_big = MAX_MODULUS + 1;              // 2^64, one past the limit
    SolveResult r = solveModPrime(s.M, s.X, 5, too_big);
    CHECK(r.status == SolveStatus::MODULUS_TOO_LARGE);
    CHECK(rankModPrime(s.M, 5, too_big) == SIZE_MAX);
}

void test_accepts_modulus_at_the_top_of_the_range() {
    // 2^64-59 is the largest prime that fits; it must be accepted, so the
    // bound check is off-by-one-free at the end that matters.
    System s = buildSystem(20, 10, 3, P_2_64);
    SolveResult r = solveModPrime(s.M, s.X, 10, P_2_64);
    CHECK(r.status == SolveStatus::OK);
}

void test_rejects_out_of_range_column_index() {
    System s = buildSystem(10, 5, 2, P_SMALL);
    s.M[0].push_back({99, 3});                   // column 99 with n_cols == 5
    SolveResult r = solveModPrime(s.M, s.X, 5, P_SMALL);
    CHECK(r.status == SolveStatus::BAD_DIMENSIONS);
    CHECK(rankModPrime(s.M, 5, P_SMALL) == SIZE_MAX);
}

void test_rejects_rhs_length_mismatch() {
    System s = buildSystem(10, 5, 2, P_SMALL);
    s.X.pop_back();
    SolveResult r = solveModPrime(s.M, s.X, 5, P_SMALL);
    CHECK(r.status == SolveStatus::BAD_DIMENSIONS);
}

void test_rejects_empty_inputs() {
    RelationMatrix empty_M;
    U128Vector empty_X;
    CHECK(solveModPrime(empty_M, empty_X, 5, P_SMALL).status == SolveStatus::BAD_DIMENSIONS);
    System s = buildSystem(10, 5, 2, P_SMALL);
    CHECK(solveModPrime(s.M, s.X, 0, P_SMALL).status == SolveStatus::BAD_DIMENSIONS);
}

// ---- Non-mutation -----------------------------------------------------------
//
// Contract: LinBox's solve consumes the matrix it is given, so lin_alg must
// copy internally and leave the caller's inputs untouched.

void test_caller_inputs_are_not_modified() {
    System s = buildSystem(40, 20, 3, P_2_64);
    RelationMatrix M_before = s.M;
    U128Vector X_before = s.X;
    solveModPrime(s.M, s.X, 20, P_2_64);
    rankModPrime(s.M, 20, P_2_64);
    CHECK(s.M == M_before);
    CHECK(s.X == X_before);
}

struct TestCase {
    const char* name;
    void (*fn)();
};

} // namespace

int main() {
    TestCase tests[] = {
        {"solve_recovers_planted_solution_small_prime", test_solve_recovers_planted_solution_small_prime},
        {"solve_is_correct_across_prime_sizes", test_solve_is_correct_across_prime_sizes},
        {"solve_randomized_full_rank_systems", test_solve_randomized_full_rank_systems},
        {"solve_handles_rectangular_overdetermined", test_solve_handles_rectangular_overdetermined},
        {"solve_reduces_unreduced_rhs", test_solve_reduces_unreduced_rhs},
        {"duplicate_column_entries_are_summed", test_duplicate_column_entries_are_summed},
        {"solve_partial_rank_locks_determined_coordinates", test_solve_partial_rank_locks_determined_coordinates},
        {"rank_is_full_for_generic_system", test_rank_is_full_for_generic_system},
        {"rank_depends_on_the_modulus", test_rank_depends_on_the_modulus},
        {"wiedemann_solves_square_systems", test_wiedemann_solves_square_systems},
        {"wiedemann_and_elimination_agree_on_row_satisfaction", test_wiedemann_and_elimination_agree_on_row_satisfaction},
        {"inconsistent_system_is_not_silently_accepted", test_inconsistent_system_is_not_silently_accepted},
        {"rejects_modulus_above_supported_range", test_rejects_modulus_above_supported_range},
        {"accepts_modulus_at_the_top_of_the_range", test_accepts_modulus_at_the_top_of_the_range},
        {"rejects_out_of_range_column_index", test_rejects_out_of_range_column_index},
        {"rejects_rhs_length_mismatch", test_rejects_rhs_length_mismatch},
        {"rejects_empty_inputs", test_rejects_empty_inputs},
        {"caller_inputs_are_not_modified", test_caller_inputs_are_not_modified},
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
