#pragma once
#include <cstdint>
#include <vector>
#include "types.h"

/**
 * Builds a binary product tree over `level` (typically the factor base):
 * returns each level bottom-up, level[0] == the input values themselves,
 * each level above pairs up adjacent entries and multiplies them (an odd
 * entry out carries through unpaired), and the last level is a single
 * value -- the product of everything in the input. This structure is what
 * makes batch smoothness testing (smoothCandidates/treeFactorize) fast.
 */
std::vector<MpzVector> buildProductTree(MpzVector level);

/**
 * Given a product tree p_levels (see buildProductTree) built over a factor
 * base, and a batch of candidates X, returns the indices into X of every
 * candidate that is fully smooth over that factor base (i.e. every prime
 * factor of X[i] appears in the base). Implements Bernstein's batch
 * smoothness test: reduce the factor-base product mod each candidate,
 * repeatedly square (mont::pow2mod) enough times to exceed the candidate's
 * bit length, then gcd against the candidate -- a full match (gcd == x)
 * means x's value divides out completely over the base.
 */
std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const U128Vector& X);

/**
 * Given a product tree p_levels and a value d already known to be smooth
 * over its base (see smoothCandidates), returns d's full factorization as
 * sparse (factor-base index, exponent) pairs, by walking down the tree with
 * gcd tests to find which primes divide d and with what multiplicity.
 */
SparseList treeFactorize(const std::vector<MpzVector>& p_levels, u128 d);
