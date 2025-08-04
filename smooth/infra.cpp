#include <gmpxx.h>
#include "smooth_algos.h"
#include "gauss_dream.h"
#include "infra.h"
#include <random>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>
#include <unordered_set>
#include <cassert>
#include <cstring>
#include <numeric>
#include <linbox/matrix/sparse-matrix.h>

#include <linbox/linbox-config.h>
#include <linbox/algorithms/block-wiedemann.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/solutions/solution-tags.h> 
#include <linbox/solutions/solve.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/error.h>
#include <linbox/solutions/solve/solve-wiedemann.h>
#include <linbox/ring/modular.h>
#include <linbox/solutions/rank.h>
#include <linbox/field/gf2.h>
#include <linbox/blackbox/zero-one.h>
#include <linbox/algorithms/gauss-gf2.h>
#include <linbox/blackbox/permutation.h> 
#include <linbox/algorithms/gauss.h>


using Field64  = Givaro::Modular<int64_t, __uint128_t>;
using Field32  = Givaro::Modular<uint32_t, uint64_t>;
template<typename Field> using Sparse = LinBox::SparseMatrix<Field>;
template<typename Field> using Vec = LinBox::DenseVector<Field>;
struct Part { MpzVector x;  mpz_class mod; };
using Seq   = LinBox::SparseMatrixFormat::SparseSeq;
using Gf2 = LinBox::GF2;


std::vector<MpzVector> buildProductTree(MpzVector level) {
    size_t lvl_size = level.size();
    size_t height = 1 + (lvl_size > 1 ? (32 - __builtin_clz((unsigned)lvl_size - 1)) : 0);
    std::vector<MpzVector> levels;
    levels.reserve(height);

    level.reserve(lvl_size);
    levels.push_back(std::move(level));
    MpzVector next_level;
    while (lvl_size > 1) {
        next_level.clear();
        next_level.reserve((lvl_size + 1) >> 1);

        const auto &cur = levels.back();
        for (size_t i = 0; i < lvl_size; i+= 2) {
            if (i + 1 < lvl_size) {
                next_level.emplace_back(cur[i] * cur[i + 1]);
            }
            else next_level.emplace_back(cur[lvl_size - 1]);
        }
        
        levels.push_back(std::move(next_level));
        lvl_size = levels.back().size();
    }
    return levels;
}

MpzVector batchRemainders(const std::vector<MpzVector>& p_levels,
                                        const MpzVector& X) {
    const mpz_class &Z = p_levels.back()[0];
    auto x_levels = buildProductTree(X);
    size_t h = x_levels.size() - 1;
    const mpz_class& root = x_levels[h][0];

    MpzVector rems, next_rems;
    next_rems.reserve(X.size());
    rems.reserve(X.size());
    rems.push_back(Z % root);
    size_t m;

    for (size_t level = h; level > 0; --level) {
        const auto& node_prods = x_levels[level - 1];
        const auto& prev_rems = rems;
        m = node_prods.size();
        next_rems.clear();
        next_rems.reserve(m);

        for (size_t idx = 0; idx < m; ++idx) {
            next_rems.emplace_back(prev_rems[idx >> 1] % node_prods[idx]);
        }
        rems.swap(next_rems);
    }

    return rems;
}

std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const MpzVector& X) {
    MpzVector rems = batchRemainders(p_levels, X);
    
    auto max_it = std::max_element(X.begin(), X.end());
    mpz_class max_x = *max_it;

    unsigned long bits = mpz_sizeinbase(max_x.get_mpz_t(), 2);
    mpz_class M;
    mpz_ui_pow_ui(M.get_mpz_t(), 2, bits);

    std::vector<size_t> smooth_cands;

    mpz_class y, g;
    for (size_t i = 0; i < X.size(); ++i) {
        const mpz_class& x_mpz = X[i];
        const mpz_class& r = rems[i];
        mpz_powm(y.get_mpz_t(), r.get_mpz_t(), M.get_mpz_t(), x_mpz.get_mpz_t());
        g = gcd(x_mpz, y);

        if (g == x_mpz) smooth_cands.push_back(i);

    }
    return smooth_cands;
}

inline uint32_t find_pow(mpz_class d, const mpz_class& P) {
    uint32_t exp = 0;
    do {
        exp++;
        d /= P;
    } while (d % P == 0);
    return exp;
}

SparseList treeFactorize(const std::vector<MpzVector>& p_levels, const mpz_class &d_mp) {
    SparseList result;
    result.reserve(12);
    std::vector<std::pair<size_t,size_t>> stack;
    stack.reserve(p_levels.size());
    stack.emplace_back(p_levels.size()-1, 0);

    mpz_class g;
    mpz_class d = d_mp;
    size_t level, idx, right, left;
    uint32_t exp;
    while(!stack.empty()) {
        if (d == 1) break;
        std::tie(level, idx) = stack.back();
        stack.pop_back();
        const mpz_class& P = p_levels[level][idx];
        g = gcd(d, P);
        if (g == 1) continue;
        if (level == 0) {
            exp = find_pow(d, P);
            result.emplace_back(idx, exp);
        } else {
            left = 2 * idx;
            right = left + 1;
            stack.push_back(std::make_pair(level - 1, left));
            if (right < p_levels[level - 1].size()) stack.push_back(std::make_pair(level - 1, right));
        }
    }
    return result;
}

template<typename Field> static void linSolveImpl(const RelationMatrix& M_rows, const MpzVector& X_col, unsigned long p_ui, MpzVector& L_out) {
    Field F(static_cast<typename Field::Element>(p_ui));

    const std::size_t m = M_rows.size();
    const std::size_t k = L_out.size();

    Sparse<Field> M(F, m, k);
    Vec<Field> b(F, m), x(F, k), xBackup(F, k);

    for (std::size_t i = 0; i < m; ++i) {
        for (auto [j, e32] : M_rows[i]) {
            typename Field::Element e; F.init(e, e32);
            M.setEntry(i, j, e);
        }
        typename Field::Element rhs; F.init(rhs, mpz_fdiv_ui(X_col[i].get_mpz_t(), p_ui));
        b[i] = rhs;
    }
    M.finalize();

    if constexpr (std::is_same_v<Field, Field64>) LinBox::solve(x, M, b, LinBox::Method::Wiedemann());
    else {
        try {
            LinBox::solve(x, M, b, LinBox::Method::SparseElimination());
        } catch (const LinBox::LinboxError& e) {
            x = xBackup;
            LinBox::solve(x, M, b, LinBox::Method::Wiedemann());
        }
    }

    for (std::size_t j = 0; j < k; ++j) {
        typename Field::Element v = x[j];
        L_out[j] = mpz_class(static_cast<unsigned long>(v));
    }
}

void linSolve(const RelationMatrix& M_rows, const MpzVector& X_col,
              const mpz_class& prime_mpz, MpzVector& L_out)
{
    unsigned long bits = mpz_sizeinbase(prime_mpz.get_mpz_t(), 2);
    unsigned long p_ui = prime_mpz.get_ui();

    if (bits > 32) linSolveImpl<Field64>(M_rows, X_col, p_ui, L_out);
    else linSolveImpl<Field32>(M_rows, X_col, p_ui, L_out);
}

inline mpz_class dot_row(const std::vector<std::pair<size_t,uint32_t>>& row, const std::vector<mpz_class>& vec) {
    mpz_class sum = 0;
    for (auto [col, coeff32] : row)
        sum += mpz_class(coeff32) * vec[col];
    return sum;
}

inline mpz_class inv_mod(const mpz_class& a, const mpz_class& m) {
    mpz_class inv;
    int ok = mpz_invert(inv.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t());
    assert(ok && "prefix and modulus must be coprime");
    return inv;
}

inline void hensel_update(MpzVector& x, const MpzVector& delta, const mpz_class& q_pow) {
    const std::size_t n = x.size();
    for (std::size_t j = 0; j < n; ++j)
        x[j] += q_pow * delta[j];
}

inline void garner_step(MpzVector& X, MpzVector& prefix, const MpzVector& xi, const mpz_class& mi) {
    const mpz_class inv = inv_mod(prefix[0], mi);
    const std::size_t n = X.size();
    for (std::size_t j = 0; j < n; ++j) {
        mpz_class t = ( (xi[j] - X[j]) * inv ) % mi;
        X[j] += prefix[j] * t;
        prefix[j] *= mi;
    }
}

std::size_t rank_relation_gf2(const RelationMatrix& rows, std::size_t k)
{
    Gf2 F;
    Gf2::Element one;
    F.init(one, true);

    LinBox::ZeroOne<Gf2> A(F, rows.size(), k);
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (const auto& cell : rows[i])
            if (cell.second & 1U)
                A.setEntry(i, cell.first, one);

    std::size_t r = 0;
    LinBox::GaussDomain<Gf2> D(F);
    D.rank(r, A);
    return r;
}

MpzVector crtSolve(RelationMatrix& M_rows, const MpzVector& X_col, const FactorList& factorList, const std::size_t k) {    
    std::vector<Part> partial;
    partial.reserve(factorList.size());

    for (const auto& [q_mpz, e] : factorList) {
        const mpz_class q = q_mpz;
        MpzVector x_q(k);
        linSolve(M_rows, X_col, q, x_q);

        mpz_class q_pow = q;
        if (e > 1) {
            MpzVector Ax(M_rows.size());
            for (uint32_t pow = 1; pow < e; ++pow) {
                for (std::size_t row = 0; row < M_rows.size(); ++row) {
                    mpz_class dot = dot_row(M_rows[row], x_q);
                    mpz_class tmp = (X_col[row] - dot) / q_pow;
                    Ax[row] = tmp % q;
                }
                MpzVector delta(k);
                linSolve(M_rows, Ax, q, delta);

                hensel_update(x_q, delta, q_pow);

                q_pow *= q;
            }
        } else q_pow = q;

        partial.push_back({ std::move(x_q),  q_pow });
    }

    const mpz_class modulus = std::accumulate(partial.begin(), partial.end(), mpz_class(1), [](mpz_class a,const Part& p){ return a * p.mod; });

    MpzVector L(k, 0);
    MpzVector prefix(k, 1);

    for (const Part& part : partial) garner_step(L, prefix, part.x, part.mod);
    for (auto& v : L) mpz_mod(v.get_mpz_t(), v.get_mpz_t(), modulus.get_mpz_t());
    return L;
}
