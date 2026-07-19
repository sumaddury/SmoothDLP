#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include "types.h"

// Include order below is load-bearing: <linbox/ring/modular.h> must precede
// the matrix/vector headers so LinBox's MVProductDomain specializations are
// visible at the point Matrix/Vector below are instantiated. Get it wrong
// and the failure is a wall of undefined symbols at LINK time, not a compile
// error. Every translation unit that includes this header inherits this
// ordering automatically -- don't reorder these.
#include <linbox/linbox-config.h>
#include <linbox/ring/modular.h>
#include <givaro/modular-ruint.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/vector/vector.h>

namespace lin_alg {

/**
  Largest modulus this module's solves accept (2^64 - 1). The underlying
  field is Givaro::Modular<RecInt::ruint<7>>, a fixed-width 128-bit
  representation whose maxCardinality is 2^64, so every residue mod q fits
  in a uint64_t and all intermediate products stay exact.

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
  The field/matrix/vector types every solve in this module operates over.
  Exposed here (not kept .cpp-local) so a caller can build its own Matrix/
  Vector once -- e.g. a Hensel-lift driver that converts a RelationMatrix to
  a Matrix a single time via its own conversion step, then calls
  solveModPrime repeatedly against the same A with a fresh b per lift level
  -- without duplicating this module's field setup or its reasoning above.
*/
using Field = Givaro::Modular<RecInt::ruint<7>>;
using Matrix = LinBox::SparseMatrix<Field, LinBox::SparseMatrixFormat::SparseSeq>;
using Vector = LinBox::DenseVector<Field>;

/**
  Constructs the GF(q) field object described above. Does not itself check
  q against MAX_MODULUS -- solveModPrime/rankModPrime check internally
  before ever touching a Field; a caller building its own Matrix ahead of a
  solveModPrime call should check q <= MAX_MODULUS first.
*/
Field makeField(u128 q);

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
  Outcome of a solve attempt. OK means LinBox wrote a solution into the
  caller's x (which may still be only partially meaningful -- checking
  which coordinates are actually determined is the caller's job, via
  rankModPrime, not this library's). INCONSISTENT means the system has no
  solution at all, which for genuine relation data means a corrupt relation
  rather than an unlucky matrix. MODULUS_TOO_LARGE and BAD_DIMENSIONS are
  caller-input rejections, checked before any LinBox call. FAILED is any
  other LinBox-reported failure.
*/
enum class SolveStatus { OK, INCONSISTENT, MODULUS_TOO_LARGE, BAD_DIMENSIONS, FAILED };

/**
  Solves A*x === b (mod q) for x, over the prime field GF(q). This is a raw,
  no-conversion layer: A and b must already be built (over a Field from
  makeField(q), or one solveModPrime built internally in a prior call with
  the same q), and x must already be sized to A.coldim() -- nothing here
  builds a Matrix/Vector from a RelationMatrix/U128Vector or reduces raw
  exponents mod q, unlike this module's previous all-in-one form. x is only
  written when status == OK.

  A is not modified: LinBox's solve() dispatch (for both methods below)
  takes A by const reference and copies it internally before any destructive
  elimination step, verified directly against the installed LinBox rather
  than assumed from its docs.

  q must be an odd prime <= MAX_MODULUS; passed explicitly only for the
  cheap bounds check below; A/x/b already carry the actual Field. Primality
  of q is a precondition and is NOT checked (checking it would pull the
  whole factoring stack into this module); a composite q yields meaningless
  output.

  This is meant to be called repeatedly against the same A with a different
  b each time -- e.g. once per Hensel-lift level -- without repaying the
  RelationMatrix-to-Matrix conversion cost on every call; that conversion is
  the caller's responsibility, done once, not this function's.
*/
SolveStatus solveModPrime(
    const Matrix& A,
    Vector& x,
    const Vector& b,
    u128 q,
    Method method = Method::SparseElimination);

/**
  Rank of M over GF(q), or SIZE_MAX if the inputs are rejected (q above
  MAX_MODULUS, or a column index >= n_cols). Diagnostic companion to
  solveModPrime: comparing the rank against n_cols tells you the dimension
  of the kernel, hence how many coordinates of a returned solution are
  arbitrary.

  Unlike solveModPrime above, this still takes the raw RelationMatrix
  directly and builds its own Matrix internally -- it's a one-shot
  diagnostic call, not part of the repeated-solve Hensel-lift path, so there
  is no matrix to reuse across calls.

  Rank is a property of M mod q, not of M -- the same integer matrix can be
  full rank modulo one prime and deficient modulo another, so this must be
  called per CRT component and its answer never reused across components.
*/
size_t rankModPrime(const RelationMatrix& M, size_t n_cols, u128 q);

/**
  Lifts a base solution of M*L === X (mod q) up to a solution mod q^e, via
  e-1 rounds of Hensel lifting. q must be prime; q^e is intended to be one
  prime-power factor of p-1 in this module's CRT/Hensel decomposition of the
  full mod-(p-1) system. Internally this is e calls to solveModPrime against
  the same q (one field, reused every level, matching solveModPrime's own
  repeated-RHS design) rather than a solve against a growing modulus -- only
  the right-hand side changes per level, never the field.

  X must be the genuine, unreduced integer right-hand side (not reduced mod
  q ahead of time): reducing it early would erase exactly the higher digits
  every level past the first needs to recover.

  Returns L mod q^e as plain integers, one per factor-base column (index
  matching n_cols/M's column indices).
*/
U128Vector henselLift(
    const RelationMatrix& M,
    const U128Vector& X,
    size_t n_cols,
    u128 q,
    int e,
    Method method = Method::SparseElimination);

} // namespace lin_alg
