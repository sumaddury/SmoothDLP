// LinBox-backed sparse solve over GF(q) for a single prime q (see lin_alg.h
// for the interface contract, and for why the field type is what it is).

#include "lin_alg.h"

#include <linbox/solutions/solve.h>
#include <linbox/solutions/rank.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/error.h>

#include <algorithm>
#include <exception>

namespace lin_alg {

const u128 MAX_MODULUS = ~static_cast<u128>(0) >> 64;   // 2^64 - 1

Field makeField(u128 q) {
    const RecInt::ruint<7> modulus(static_cast<uint64_t>(q));
    return Field(modulus);
}

namespace {

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
            const size_t col = row[k].first;
            u128 sum = 0;
            while (k < row.size() && row[k].first == col) {
                sum = (sum + static_cast<u128>(row[k].second)) % q;
                ++k;
            }
            if (sum == 0) continue;
            Field::Element e;
            F.init(e, static_cast<uint64_t>(sum));
            A.setEntry(i, col, e);
        }
    }
    A.finalize();
}

/**
  Fills a LinBox dense vector from V, reducing each entry mod q. The
  right-hand-side counterpart to fillMatrix: relation exponents (t values)
  arrive as unreduced u128 integers (they are exponents mod p-1, not mod q),
  so this is where that reduction actually happens for a given q.
*/
void fillVector(Vector& v, const Field& F, const U128Vector& V, u128 q) {
    for (size_t i = 0; i < V.size(); ++i) {
        Field::Element e;
        F.init(e, static_cast<uint64_t>(V[i] % q));
        v.setEntry(i, e);
    }
}

/**
  Reads a LinBox vector back out as plain u128 residues mod q -- the inverse
  of fillVector. Every entry comes back in [0, q), since that is all a
  Field::Element can ever represent; anything beyond a single mod-q digit
  (e.g. accumulating a solution across Hensel-lift levels) is the caller's
  responsibility, not this function's.
*/
U128Vector getVector(const Field& F, const Vector& v) {
    U128Vector out(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        uint64_t val;
        F.convert(val, v[i]);
        out[i] = static_cast<u128>(val);
    }
    return out;
}

} // namespace

SolveStatus solveModPrime(const Matrix& A, Vector& x, const Vector& b, u128 q, Method method) {
    if (q < 2 || q > MAX_MODULUS) return SolveStatus::MODULUS_TOO_LARGE;
    if (A.rowdim() == 0 || A.coldim() == 0 ||
        A.rowdim() != b.size() || A.coldim() != x.size())
        return SolveStatus::BAD_DIMENSIONS;

    try {
        if (method == Method::Wiedemann)
            LinBox::solve(x, A, b, LinBox::Method::Wiedemann());
        else
            LinBox::solve(x, A, b, LinBox::Method::SparseElimination());
    } catch (const LinBox::LinboxMathInconsistentSystem&) {
        return SolveStatus::INCONSISTENT;
    } catch (const LinBox::LinboxError&) {
        // NOTE: LinboxError does NOT derive from std::exception, so it needs
        // its own handler -- a lone catch(std::exception) silently misses
        // every LinBox failure and lets it escape as an uncaught throw.
        return SolveStatus::FAILED;
    } catch (const std::exception&) {
        return SolveStatus::FAILED;
    } catch (...) {
        return SolveStatus::FAILED;
    }

    return SolveStatus::OK;
}

size_t rankModPrime(const RelationMatrix& M, size_t n_cols, u128 q) {
    if (q < 2 || q > MAX_MODULUS) return SIZE_MAX;
    if (n_cols == 0 || M.empty() || !dimensionsValid(M, n_cols, nullptr)) return SIZE_MAX;

    const Field F = makeField(q);
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
    return static_cast<size_t>(r);
}

namespace {

/**
  Exact-integer sparse matrix-vector product M*L, with NO modular reduction
  -- unlike fillMatrix, which reduces mod q for the GF(q) solve itself, this
  is used by henselLift to compute the residual (X - M*L_i) that determines
  the next lift digit, which needs the true integer value, not a reduction
  mod q that has already thrown away the higher digits.
*/
U128Vector modSpMV(const RelationMatrix& M, const U128Vector& L) {
    U128Vector out(M.size());
    for (size_t i = 0; i < M.size(); ++i) {
        u128 accum = 0;
        for (const auto& [j, exponent] : M[i]) accum += L[j] * exponent;
        out[i] = accum;
    }
    return out;
}

} // namespace

U128Vector henselLift(
    const RelationMatrix& M,
    const U128Vector& X,
    size_t n_cols,
    u128 q,
    int e,
    Method method) {

    const Field F = makeField(q);
    Matrix A(F, M.size(), n_cols);
    fillMatrix(A, F, M, q);

    Vector b(F, X.size());
    fillVector(b, F, X, q);

    Vector x(F, n_cols);
    solveModPrime(A, x, b, q, method);

    U128Vector L = getVector(F, x);
    u128 q_k = q;

    for (int k = 1; k < e; k++) {
        U128Vector ML_k = modSpMV(M, L);

        // (X[i] - ML_k[i]) is only guaranteed >= 0 as an unbounded integer,
        // not as u128 arithmetic -- nothing orders X[i] against ML_k[i]
        // directly. Reduce both mod q^(k+1) first (the residual is only
        // ever meaningful at that precision anyway, see lin_alg.h), then
        // subtract with an explicit wraparound instead of letting unsigned
        // subtraction silently underflow.
        const u128 q_k1 = q_k * q;
        U128Vector R_k(M.size());
        for (size_t i = 0; i < M.size(); ++i) {
            const u128 xr = X[i] % q_k1;
            const u128 mr = ML_k[i] % q_k1;
            const u128 diff = (xr >= mr) ? (xr - mr) : (xr + q_k1 - mr);
            R_k[i] = diff / q_k;
        }

        Vector r(F, X.size());
        fillVector(r, F, R_k, q);

        Vector delta(F, n_cols);
        solveModPrime(A, delta, r, q, method);

        U128Vector D = getVector(F, delta);
        for (size_t col = 0; col < n_cols; ++col) L[col] += q_k * D[col];

        q_k = q_k1;
    }

    return L;
}

} // namespace lin_alg
