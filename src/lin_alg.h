#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include "types.h"

namespace lin_alg {

/**
  Largest modulus solveModPrime/rankModPrime accept (2^64 - 1). The
  underlying field is Givaro::Modular<RecInt::ruint<7>>, a fixed-width
  128-bit representation whose maxCardinality is 2^64, so every residue
  mod q fits in a uint64_t and all intermediate products stay exact.

  This bound is a hard correctness boundary, not a performance one: the
  obvious-looking alternative Givaro::Modular<uint64_t, __uint128_t>
  advertises maxCardinality 2^64-1 and has correct field arithmetic, but
  LinBox's sparse elimination silently returns WRONG solutions with it for
  q above roughly 2^32 -- no exception, no warning, just a bad vector.
  Modular<int64_t, __uint128_t> fails identically. Do not "optimize" this
  file by switching to either of them; see tests/internal/test_lin_alg.cpp,
  which pins the working field against a full 64-bit prime.
*/
extern const u128 MAX_MODULUS;

/**
  Which LinBox algorithm solveModPrime/rankModPrime dispatch to.

  SparseElimination is structured Gaussian elimination (LinBox's
  GaussDomain): faster for the matrix sizes this project reaches (measured
  ~2x faster at n<=2000), tolerates rectangular systems, but suffers
  fill-in so its cost grows faster than Wiedemann's.

  Wiedemann is the Krylov/black-box method: O(n) working memory and cost
  driven by matrix-vector products rather than fill-in, so it overtakes
  elimination on larger systems (measured crossover near n=4000). It
  expects a square system -- pass a square M if you select it.
*/
enum class Method { SparseElimination, Wiedemann };

/**
  Outcome of a solve attempt. OK means LinBox returned a vector (which may
  still be only partially meaningful -- see SolveResult's own doc and the
  rank discussion there; checking which coordinates are actually
  determined is the caller's job, via rankModPrime, not this library's).
  INCONSISTENT means the system has no solution at all, which for genuine
  relation data means a corrupt relation rather than an unlucky matrix.
  MODULUS_TOO_LARGE and BAD_DIMENSIONS are caller-input rejections, checked
  before any LinBox call. FAILED is any other LinBox-reported failure.
*/
enum class SolveStatus { OK, INCONSISTENT, MODULUS_TOO_LARGE, BAD_DIMENSIONS, FAILED };

/**
  Result of solveModPrime. L has length n_cols and is only populated when
  status == OK; it holds one residue mod q per factor-base column.

  When status is OK, L solves the system. What that does NOT mean is that
  every coordinate of L is uniquely determined: if M is rank-deficient mod
  q, the solution space is a coset of a nonzero kernel, and coordinates
  lying in that kernel are arbitrary while still satisfying every row.
  Sorting the determined coordinates from the arbitrary ones happens
  downstream by re-exponentiation (g^L_j == p_j mod p) after CRT
  recombination, which cannot be done here -- it needs L mod (p-1), not L
  mod q.
*/
struct SolveResult {
  SolveStatus status;
  U128Vector L;
};

/**
  Solves the sparse system M*L === X (mod q) for L, over the prime field
  GF(q). M is given by rows, each a SparseList of (column index, exponent)
  pairs with column indices < n_cols; X is the dense right-hand side, one
  entry per row of M. Entries of both are reduced mod q internally, so the
  caller may pass raw unreduced values (X's entries are the relation
  exponents t, which are reduced mod p-1, not mod q).

  M may be rectangular -- in relation collection it is overdetermined
  (more relations than factor-base columns), which is the expected shape.
  Rank deficiency is not an error: LinBox still returns a valid solution
  of the system, and empirically every uniquely-determined coordinate of
  that solution is exactly correct, with only kernel coordinates arbitrary.
  That is what makes partial solves usable here.

  q must be an odd prime <= MAX_MODULUS. Primality is a precondition and is
  NOT checked (checking it would pull the whole factoring stack into this
  translation unit); a composite q yields meaningless output.

  The caller's M and X are not modified. LinBox's solve consumes the matrix
  it is handed, so this builds its own copy internally.
*/
SolveResult solveModPrime(
    const RelationMatrix& M,
    const U128Vector& X,
    size_t n_cols,
    u128 q,
    Method method = Method::SparseElimination);

/**
  Rank of M over GF(q), or SIZE_MAX if the inputs are rejected (q above
  MAX_MODULUS, or a column index >= n_cols). Diagnostic companion to
  solveModPrime: comparing the rank against n_cols tells you the dimension
  of the kernel, hence how many coordinates of a returned solution are
  arbitrary.

  Rank is a property of M mod q, not of M -- the same integer matrix can be
  full rank modulo one prime and deficient modulo another, so this must be
  called per CRT component and its answer never reused across components.
*/
size_t rankModPrime(const RelationMatrix& M, size_t n_cols, u128 q);

} // namespace lin_alg
