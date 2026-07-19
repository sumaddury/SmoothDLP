// LinBox-backed sparse solve over GF(q) for a single prime q (see lin_alg.h
// for the interface contract, and for why the field type is what it is).
//
// Include order below is load-bearing: <linbox/ring/modular.h> must precede
// the matrix headers so LinBox's MVProductDomain specializations are visible
// at the point the matrix templates are instantiated. Get it wrong and the
// failure is a wall of undefined symbols at LINK time, not a compile error.

#include "lin_alg.h"

#include <linbox/linbox-config.h>
#include <linbox/ring/modular.h>
#include <givaro/modular-ruint.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/vector/vector.h>
#include <linbox/solutions/solve.h>
#include <linbox/solutions/rank.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/error.h>

#include <algorithm>
#include <exception>

namespace lin_alg {

const u128 MAX_MODULUS = ~(u128)0 >> 64;   // 2^64 - 1

namespace {

using Field = Givaro::Modular<RecInt::ruint<7>>;
using Matrix = LinBox::SparseMatrix<Field, LinBox::SparseMatrixFormat::SparseSeq>;
using Vector = LinBox::DenseVector<Field>;

/**
  Constructs the GF(q) field object. q is known to fit in 64 bits by the
  time this is called, so the ruint<7> is built from a uint64_t. Written
  with a named intermediate rather than Field(RecInt::ruint<7>(q)), which
  the most-vexing-parse rule would read as a function declaration.
*/
Field makeField(u128 q) {
    RecInt::ruint<7> modulus((uint64_t)q);
    return Field(modulus);
}

/**
  True iff every column index in M is < n_cols and X (when given) has one
  entry per row of M -- i.e. the shapes the LinBox calls are about to
  assume actually hold. Checked up front because an out-of-range column
  index would otherwise be an out-of-bounds write inside LinBox rather
  than a diagnosable error.
*/
bool dimensionsValid(const RelationMatrix& M, size_t n_cols, const U128Vector* X) {
    if (X != nullptr && X->size() != M.size()) return false;
    for (const SparseList& row : M)
        for (const auto& entry : row)
            if (entry.first >= n_cols) return false;
    return true;
}

/**
  Fills a LinBox sparse matrix from M, reducing each exponent mod q and
  dropping entries that vanish (an exponent divisible by q is genuinely
  zero in GF(q), so storing it would only cost space).

  Repeated column indices within a row are summed, which is the standard
  reading of a sparse coordinate list and the one the residual check uses.
  This matters because LinBox's setEntry OVERWRITES rather than
  accumulates: feeding it a row with a repeated column directly would
  silently build a different matrix than the caller described, and the
  system would then look inconsistent for no visible reason. Relations
  produced by treeFactorize always have distinct indices, so this is
  defensive rather than load-bearing -- but it is cheap, and the failure
  it prevents is invisible.

  Sorting each row also hands SparseSeq its entries in column order, which
  is the order that format wants.
*/
void fillMatrix(Matrix& A, const Field& F, const RelationMatrix& M, u128 q) {
    SparseList row;
    for (size_t i = 0; i < M.size(); ++i) {
        row = M[i];
        std::sort(row.begin(), row.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        for (size_t k = 0; k < row.size();) {
            size_t col = row[k].first;
            u128 sum = 0;
            while (k < row.size() && row[k].first == col) {
                sum = (sum + (u128)row[k].second) % q;
                ++k;
            }
            if (sum == 0) continue;
            Field::Element e;
            F.init(e, (uint64_t)sum);
            A.setEntry(i, col, e);
        }
    }
    A.finalize();
}

/**
  Rows of M*L === X (mod q) that the candidate L fails to satisfy, computed
  directly from the caller's M rather than from any LinBox structure -- so
  this is an independent check on LinBox's answer, not a restatement of it.

  All arithmetic stays exact in u128: both factors are reduced below
  q <= 2^64-1 before multiplying, so the product is at most (2^64-1)^2,
  which is strictly below 2^128.
*/
std::vector<size_t> residualRows(const RelationMatrix& M, const U128Vector& X,
                                 const U128Vector& L, u128 q) {
    std::vector<size_t> bad;
    for (size_t i = 0; i < M.size(); ++i) {
        u128 acc = 0;
        for (const auto& entry : M[i]) {
            u128 coeff = (u128)entry.second % q;
            u128 term  = (coeff * L[entry.first]) % q;
            acc = (acc + term) % q;
        }
        if (acc != X[i] % q) bad.push_back(i);
    }
    return bad;
}

} // namespace

SolveResult solveModPrime(const RelationMatrix& M, const U128Vector& X,
                          size_t n_cols, u128 q, Method method) {
    SolveResult result;
    result.status = SolveStatus::FAILED;

    if (q < 2 || q > MAX_MODULUS) {
        result.status = SolveStatus::MODULUS_TOO_LARGE;
        return result;
    }
    if (n_cols == 0 || M.empty() || !dimensionsValid(M, n_cols, &X)) {
        result.status = SolveStatus::BAD_DIMENSIONS;
        return result;
    }

    Field F = makeField(q);
    Matrix A(F, M.size(), n_cols);
    fillMatrix(A, F, M, q);

    Vector b(F, M.size());
    for (size_t i = 0; i < M.size(); ++i) {
        Field::Element e;
        F.init(e, (uint64_t)(X[i] % q));
        b.setEntry(i, e);
    }

    Vector x(F, n_cols);
    try {
        if (method == Method::Wiedemann)
            LinBox::solve(x, A, b, LinBox::Method::Wiedemann());
        else
            LinBox::solve(x, A, b, LinBox::Method::SparseElimination());
    } catch (const LinBox::LinboxMathInconsistentSystem&) {
        result.status = SolveStatus::INCONSISTENT;
        return result;
    } catch (const LinBox::LinboxError&) {
        // NOTE: LinboxError does NOT derive from std::exception, so it needs
        // its own handler -- a lone catch(std::exception) silently misses
        // every LinBox failure and lets it escape as an uncaught throw.
        result.status = SolveStatus::FAILED;
        return result;
    } catch (const std::exception&) {
        result.status = SolveStatus::FAILED;
        return result;
    } catch (...) {
        result.status = SolveStatus::FAILED;
        return result;
    }

    result.L.resize(n_cols);
    for (size_t j = 0; j < n_cols; ++j) {
        uint64_t v = 0;
        F.convert(v, x[j]);
        result.L[j] = (u128)v;
    }
    result.status = SolveStatus::OK;
    result.unsatisfied_rows = residualRows(M, X, result.L, q);
    return result;
}

size_t rankModPrime(const RelationMatrix& M, size_t n_cols, u128 q) {
    if (q < 2 || q > MAX_MODULUS) return SIZE_MAX;
    if (n_cols == 0 || M.empty() || !dimensionsValid(M, n_cols, nullptr)) return SIZE_MAX;

    Field F = makeField(q);
    Matrix A(F, M.size(), n_cols);
    fillMatrix(A, F, M, q);

    unsigned long r = 0;
    try {
        LinBox::rank(r, A, LinBox::Method::SparseElimination());
    } catch (const LinBox::LinboxError&) {   // see note in solveModPrime
        return SIZE_MAX;
    } catch (const std::exception&) {
        return SIZE_MAX;
    } catch (...) {
        return SIZE_MAX;
    }
    return (size_t)r;
}

} // namespace lin_alg
