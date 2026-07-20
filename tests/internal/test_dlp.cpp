// Standalone correctness tests for src/dlp.{h,cpp}: the DLP class, which
// drives infra.{h,cpp}'s building blocks (addRelations, crtSolve,
// filterDetermined, modInv) through the paper's Algorithm 1 end to end.
//
// Unlike test_infra.cpp, this file does NOT unity-include dlp.cpp or
// infra.cpp -- dlp.cpp itself only includes dlp.h (declarations), not
// infra.cpp, so DLP::solve's definition and infra::*'s definitions live in
// genuinely separate translation units here, exactly the way a real
// consumer of this library would link against them. (This distinction is
// exactly what surfaced two real bugs during development: infra::modInv was
// briefly marked `inline` with its only body in infra.cpp -- invisible from
// test_infra.cpp's unity build, but a hard link failure the moment
// something links against infra.cpp as an ordinary object file, which is
// what this file does.) See this directory's Makefile for the exact source
// list and LinBox-aware build flags this needs (the same ones test_infra.cpp
// needs, since infra.cpp still unity-includes lin_alg.cpp internally).

#include "dlp.h"
#include "gauss_dream.h"

#include <iostream>
#include <random>
#include <string>
#include <unordered_map>

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

// Real order check via gauss::factorize(p-1) -- not a restatement of
// anything DLP/infra does internally, since neither exposes a primitive-
// root test. candidate is a primitive root mod p iff candidate^((p-1)/q) is
// never 1 (mod p) for any prime q dividing p-1.
bool is_primitive_root(u128 candidate, u128 p) {
    const u128 order = p - 1;
    const FactorList f = gauss::factorize(order);
    for (auto& [q, e] : f) {
        (void)e;
        const u128 exp = order / q;
        u128 result = 1, base = candidate % p, e2 = exp;
        while (e2) {
            if (e2 & 1) result = mont::mulmod_any(result, base, p);
            base = mont::mulmod_any(base, base, p);
            e2 >>= 1;
        }
        if (result == 1) return false;
    }
    return true;
}

u128 find_primitive_root(u128 p, u128 skip) {
    for (u128 cand = 2; cand < p; ++cand)
        if (cand != skip && is_primitive_root(cand, p)) return cand;
    return 0;
}

// Brute-force log_g(b) mod p by direct enumeration -- a real, independent
// oracle (does not call powmod_odd or anything else DLP::solve itself
// uses), only practical for the smaller primes below.
u128 brute_force_discrete_log(u128 g, u128 b, u128 p) {
    u128 val = 1;
    for (u128 t = 0; t < p - 1; ++t) {
        if (val == b) return t;
        val = mont::mulmod_any(val, g, p);
    }
    return (u128)-1;  // not found -- would mean g isn't actually a generator of b
}

// ---- DLP::solve -------------------------------------------------------
//
// Contract: DLP(p).solve(g, b) returns x in [0, p-2] with g^x === b (mod p),
// given g and b are both genuine primitive roots mod p (dlp.h's documented
// precondition). Every test below either cross-checks against an
// independent brute-force oracle (small/medium primes) or, where brute
// force isn't practical, verifies the defining property g^x === b (mod p)
// directly -- which is itself a complete, rigorous correctness check, not
// a smoke test: mont::mulmod_any-based exponentiation here is independent
// of the Montgomery-form powmod_odd path DLP::solve's own relation
// collection uses internally.

u128 verify_pow(u128 base, u128 e, u128 n) {
    u128 result = 1, b = base % n;
    while (e) {
        if (e & 1) result = mont::mulmod_any(result, b, n);
        b = mont::mulmod_any(b, b, n);
        e >>= 1;
    }
    return result;
}

void test_dlp_solve_matches_brute_force_small_prime() {
    u128 p = 1009;
    u128 g = find_primitive_root(p, 0);
    u128 b = find_primitive_root(p, g);
    CHECK(g != 0 && b != 0 && g != b);

    const u128 want = brute_force_discrete_log(g, b, p);
    CHECK(want != (u128)-1);

    dlp::DLP solver(p);
    const u128 got = solver.solve(g, b);
    CHECK(got == want);
    CHECK(verify_pow(g, got, p) == b);
}

void test_dlp_solve_matches_brute_force_across_prime_sizes() {
    // A spread matching the shapes exercised in test_infra.cpp's own
    // problem-params/crtSolve sweeps -- p-1's factorization varies a lot
    // across these (1008 = 2^4*3^2*7 vs 8190 = 2*3^2*5*7*13, etc.), which
    // is exactly what exercises different CRT/Hensel-lift paths inside
    // crtSolve on both the g-side and b-side solves DLP::solve runs.
    const u128 primes[] = {(u128)1009, (u128)8191, (u128)524287};
    for (u128 p : primes) {
        u128 g = find_primitive_root(p, 0);
        u128 b = find_primitive_root(p, g);
        CHECK(g != 0 && b != 0 && g != b);

        const u128 want = brute_force_discrete_log(g, b, p);
        CHECK(want != (u128)-1);

        dlp::DLP solver(p);
        const u128 got = solver.solve(g, b);
        CHECK(got == want);
        CHECK(verify_pow(g, got, p) == b);
    }
}

void test_dlp_solve_self_consistent_for_large_prime() {
    // p ~ 2^31: brute force (O(p) group operations) is no longer practical,
    // but g^x === b (mod p) is still a complete correctness certificate on
    // its own -- there's exactly one x in [0, p-2] satisfying it, since g
    // is a primitive root, so a matching check here is not weaker evidence
    // than the brute-force cross-check above, just a different (and here,
    // necessary) way of establishing the same fact.
    u128 p = 2147483647ULL;
    u128 g = find_primitive_root(p, 0);
    u128 b = find_primitive_root(p, g);
    CHECK(g != 0 && b != 0 && g != b);

    dlp::DLP solver(p);
    const u128 got = solver.solve(g, b);
    CHECK(got != 0);  // g,b both primitive roots => genuine x is never 0, see dlp.h
    CHECK(verify_pow(g, got, p) == b);
}

void test_dlp_solve_repeated_calls_on_same_instance_are_stable() {
    // The g-side/b-side crtSolve calls inside one solve() run concurrently
    // via std::async, serialized against each other by a mutex specifically
    // because calling crtSolve from two threads at once was found (by
    // direct testing during development) to crash the underlying LinBox/
    // Givaro/NTL/FFLAS-FFPACK stack nondeterministically -- sometimes after
    // several successful calls, not on the first one. A single solve()
    // call is a weak regression test for that: this repeats solve() on one
    // shared DLP instance (reusing its ProblemParams, the way a real caller
    // solving multiple targets against the same p would) enough times that
    // an intermittent crash would very likely have shown up if the mutex
    // fix ever regressed.
    u128 p = 1009;
    u128 g = find_primitive_root(p, 0);
    u128 b = find_primitive_root(p, g);
    CHECK(g != 0 && b != 0 && g != b);

    dlp::DLP solver(p);
    for (int i = 0; i < 15; ++i) {
        const u128 got = solver.solve(g, b);
        CHECK(verify_pow(g, got, p) == b);
    }
}

void test_dlp_solve_inverse_relationship_when_swapping_roles() {
    // If b = g^x then g = b^y with x*y === 1 (mod p-1) -- true here
    // specifically because b is required to be a primitive root too (see
    // dlp.h's precondition), which forces gcd(x, p-1) = 1 (a value is a
    // primitive root's image under exponentiation by x iff x is itself a
    // unit mod p-1). Solving both directions and checking this cross-
    // relationship is a correctness check that doesn't depend on any
    // independent oracle at all -- an internally-consistent-but-wrong pair
    // of answers (e.g. both off by the same corrupting factor) would still
    // be caught here, unlike a check that only ever verifies one direction.
    u128 p = 8191;
    u128 g = find_primitive_root(p, 0);
    u128 b = find_primitive_root(p, g);
    CHECK(g != 0 && b != 0 && g != b);

    dlp::DLP solver(p);
    const u128 x = solver.solve(g, b);
    const u128 y = solver.solve(b, g);

    CHECK(verify_pow(g, x, p) == b);
    CHECK(verify_pow(b, y, p) == g);
    CHECK(mont::mulmod_any(x, y, p - 1) == 1);
}

struct TestCase {
    const char* name;
    void (*fn)();
};

}  // namespace

int main() {
    TestCase tests[] = {
        {"dlp_solve_matches_brute_force_small_prime", test_dlp_solve_matches_brute_force_small_prime},
        {"dlp_solve_matches_brute_force_across_prime_sizes", test_dlp_solve_matches_brute_force_across_prime_sizes},
        {"dlp_solve_self_consistent_for_large_prime", test_dlp_solve_self_consistent_for_large_prime},
        {"dlp_solve_repeated_calls_on_same_instance_are_stable", test_dlp_solve_repeated_calls_on_same_instance_are_stable},
        {"dlp_solve_inverse_relationship_when_swapping_roles", test_dlp_solve_inverse_relationship_when_swapping_roles},
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
