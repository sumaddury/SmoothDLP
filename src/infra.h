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
  base (buildProductTree) and p_factorization is p's own factorization
  (needed downstream for CRT over p-1) -- both computed once here and owned
  by value: nothing outside ProblemParams holds this data, so a reference
  member would dangle the instant the constructor returned.
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
*/
void addRelations(
    u128 base,
    const ProblemParams& params,
    RelationMatrix& M,
    U128Vector& X);
  

} // namespace infra
