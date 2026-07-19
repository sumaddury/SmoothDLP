#pragma once
#include <cstdint>
#include <vector>
#include "types.h"

namespace infra {

/**
  Builds a binary product tree over `level` (typically the factor base):
  returns each level bottom-up, level[0] == the input values themselves,
  each level above pairs up adjacent entries and multiplies them (an odd
  entry out carries through unpaired), and the last level is a single
  value -- the product of everything in the input. This structure is what
  makes batch smoothness testing (smoothCandidates/treeFactorize) fast.
*/
std::vector<MpzVector> buildProductTree(MpzVector level);

/**
  Given a product tree p_levels (see buildProductTree) built over a factor
  base, and a batch of candidates X, returns the indices into X of every
  candidate that is fully smooth over that factor base (i.e. every prime
  factor of X[i] appears in the base). Implements Bernstein's batch
  smoothness test: reduce the factor-base product mod each candidate,
  repeatedly square (mont::pow2mod) enough times to exceed the candidate's
  bit length, then gcd against the candidate -- a full match (gcd == x)
  means x's value divides out completely over the base.
*/
std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const U128Vector& X);

/**
  Given a product tree p_levels and a value d already known to be smooth
  over its base (see smoothCandidates), returns d's full factorization as
  sparse (factor-base index, exponent) pairs, by walking down the tree with
  gcd tests to find which primes divide d and with what multiplicity.
*/
SparseList treeFactorize(const std::vector<MpzVector>& p_levels, u128 d);

/**
  Fixed, precomputed inputs shared by every addRelations call for one side
  (g-side or b-side) of one discrete-log instance: p is the prime modulus;
  p_prime/r2 are p's Montgomery constants (mont::inverse(p)/
  mont::r_squared_mod_n(p)), hoisted once for reuse across every powmod_odd
  call rather than recomputed per relation; mask_bitlength/mask describe
  the random-exponent rejection-sampling range (mask_bitlength is p-2's bit
  count, mask the smallest all-ones bitmask covering [0, p-2]); B is the
  smoothness bound the factor base is built with; smooth_density is the
  estimated probability that a random residue mod p is B-smooth (see
  smooth_algos.h's logDickman -- this is the actual probability, not its
  log, since addRelations divides by it directly), used to size relation-
  collection batches; p_levels is the product tree built over the factor
  base (buildProductTree) and p_factorization is (p-1)'s factorization, not
  p's -- it's the group order's prime-power factors that crtSolve's CRT/
  Hensel decomposition actually runs over -- both computed once here and
  owned by value: nothing outside ProblemParams holds this data, so a
  reference member would dangle the instant the constructor returned.
*/
struct ProblemParams {
  const u128 p, p_prime, r2;
  const int mask_bitlength;
  const u128 mask, B;
  const double smooth_density;
  const std::vector<MpzVector> p_levels;
  const FactorList p_factorization;

  ProblemParams(u128 p);
};

/**
  Given a product tree p_levels over the factor base, and base g or b, adds
  smooth relations to M and X until M has grown by at least k rows, where k
  is the factor base size (params.p_levels[0].size()). Each new row is a
  sparse (factor-base index, exponent) factorization of some B-smooth
  base^t mod p, and the matching entry appended to X is the exponent t
  itself (not the smooth value) -- X is the right-hand side of the linear
  system M*L === X (mod p-1) that recovers base's discrete logs, so it must
  hold the t's, which is what that system is actually solving for.

  PRECONDITION, not checked here: base must be a genuine primitive root mod
  p (order exactly p-1), for either the g or b role. If it isn't, base^t
  only ever ranges over the proper subgroup base generates, which caps what
  M*L === X (mod p-1) can ever learn about factor-base primes outside that
  subgroup's structure -- confirmed directly: for such a base, crtSolve can
  hit a rank ceiling on some factor of p-1 that persists no matter how many
  relations are added (100+ relations past the minimum, still deficient),
  where a confirmed primitive root on the identical (p, factor base) clears
  it within a small multiple of k every time. This module does not verify
  the precondition (the same way std::lower_bound does not verify its range
  is sorted) -- an invalid base fails silently as persistent crtSolve
  failure, not as a diagnosable error here.
*/
void addRelations(
    u128 base,
    const ProblemParams& params,
    RelationMatrix& M,
    U128Vector& X);

/**
  Solves M*L === X (mod p-1) for L, given a matrix/right-hand-side pair
  already accumulated by addRelations (so M's columns index params.p_levels'
  factor base, and X holds the relation exponents t). p-1 is composite, so
  this isn't a single field solve: params.p_factorization gives p-1's
  distinct prime-power factors q_i^e_i, each solved independently via
  lin_alg::henselLift (a base GF(q_i) solve, Hensel-lifted up to q_i^e_i),
  and the per-factor results are then recombined into one solution mod p-1
  via Garner's CRT algorithm.

  Returns an EMPTY vector if henselLift fails for any single factor --
  i.e. if it cannot reach that factor's full q_i^e_i precision from this
  relation set (see lin_alg.h's henselLift doc for the mechanism: a
  rank-deficient base solve whose kernel direction has raw, unreduced
  coefficients that are divisible by q_i but not q_i^2, so it looks free at
  the base level but is not actually free once lifted). This is a real,
  data-dependent possibility, not a bug to route around silently -- there is
  no partial/best-effort L worth returning in that case.

  Two different things can cause the rank deficiency behind that failure,
  and they resolve completely differently:
    - Ordinary statistical rank deficiency (a relation set that just hasn't
      accumulated enough independent rows yet): "collect more relations for
      this base and try again" is a real, converging fix. Measured directly
      in this module's test suite: 0% failure, stable, by roughly 3x the
      factor-base size k, across every prime and factor shape tested, and
      it never regresses back to failing once clear.
    - base not being a genuine primitive root (see addRelations' doc): this
      produces a PERMANENT rank ceiling that retrying can never clear, no
      matter how many relations are added -- also measured directly (100+
      relations past the minimum, still deficient on the same factor every
      time). If crtSolve keeps returning empty well past a few multiples of
      k, this is the thing to check, not "collect even more relations".

  As with lin_alg::solveModPrime/henselLift, rank deficiency mod any q_i
  that henselLift DOES successfully lift through is not an error --
  coordinates outside every factor's kernel come back correct, others come
  back arbitrary. Whether a specific coordinate (e.g. a candidate overlap
  prime's discrete log) is trustworthy mod the full p-1 is the caller's job
  to check, the same way it always has been upstream of this function --
  CRT-combining does not manufacture certainty an individual factor's solve
  didn't have.
*/
U128Vector crtSolve(
    const ProblemParams& params,
    const RelationMatrix& M,
    const U128Vector& X);

} // namespace infra
