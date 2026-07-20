#pragma once
#include "infra.h"
#include "montgomery.h"

namespace dlp {

/**
  Top-level driver for the paper's Double Index Calculus Algorithm (Huang,
  Zhang, Zhao, Peng, Liao, Wang, 2024, Algorithm 1): given a prime p, finds
  x such that g^x === b (mod p), by collecting smooth relations for g and b
  *simultaneously* against a shared factor base and combining as soon as a
  factor-base prime's discrete log is known on both sides -- the building
  blocks (infra::addRelations, infra::crtSolve, infra::filterDetermined,
  infra::modInv) are described individually in infra.h; this class is just
  the loop that ties them into the paper's algorithm end to end.

  A DLP instance owns one p's infra::ProblemParams (factor base, product
  tree, p-1's factorization -- all expensive to build) so it can be reused
  across multiple solve() calls against the same p without rebuilding them.
*/
class DLP {
private:
  infra::ProblemParams params;

public:
  DLP(u128 p);

  /**
    Returns x in [0, p-2] such that g^x === b (mod p).

    PRECONDITION, not checked here: BOTH g and b must be genuine primitive
    roots mod p (order exactly p-1) -- the same precondition
    infra::addRelations documents on its own `base` argument, required here
    twice over since g and b are each passed to addRelations in that role,
    once per side. If either isn't, that side's relation-collection process
    hits addRelations' documented permanent-rank-ceiling failure mode, and
    solve() never terminates within its relation-count cap (it will exhaust
    every doubling step and return 0 -- see below -- not throw or hang).

    Relation collection proceeds in doubling rounds -- target row counts
    k, 2k, 4k, ..., 16k (k = params.factor_base.size()) -- collecting more
    for both g and b and re-attempting crtSolve each round rather than
    committing to a single relation count up front, since (per infra.h's
    own characterization) a relation set right at k is unreliable but one
    a few multiples past it essentially always succeeds. Each round, the
    g-side and b-side pipelines (addRelations + crtSolve + filterDetermined)
    run concurrently via std::async -- addRelations is safe to call this
    way (pure GMP/Montgomery arithmetic on disjoint per-side state), but
    crtSolve is NOT: the underlying LinBox/Givaro/NTL/FFLAS-FFPACK stack
    has been confirmed (by direct testing -- two threads calling crtSolve
    concurrently crash within a handful of trials, reliably, regardless of
    input) to have thread-unsafe global state somewhere in that chain. So
    the two crtSolve calls are serialized against each other through a
    file-local mutex (dlp.cpp) even though everything else about the two
    sides runs in parallel. Do not remove that mutex to "simplify" this --
    it is load-bearing, not a leftover.

    Once both sides solve in a given round, this implements the paper's
    second phase directly: walk params.factor_base's directly-verified
    Omega_g/Omega_b pairs (see infra::filterDetermined) in the sorted order
    filterDetermined already returns them in, two-pointer-merge to find
    overlapping primes, and for the first overlap whose b-side log (beta)
    is actually invertible mod p-1 (infra::modInv -- p-1 is always even, so
    an arbitrary overlap prime's beta is not guaranteed invertible, and one
    that isn't is simply skipped in favor of the next overlap), return
    x = alpha * beta^-1 (mod p-1).

    Returns 0 if no round up to 16k produces both a solve and an invertible
    overlap. This is a safe, unambiguous "not found" sentinel, not a
    possible genuine answer: since b is required to be a primitive root
    (see the precondition above), x = 0 would mean b === g^0 === 1 (mod p),
    which is not a primitive root (order 1, not p-1) -- so a genuine x is
    never 0 under this function's own precondition.
  */
  u128 solve(u128 g, u128 b);
};

}
