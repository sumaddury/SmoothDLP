#include <gmpxx.h>
#include "smooth_algos.h"
#include "gauss_dream.h"
#include <random>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>
#include <unordered_set>
#include <cassert>
#include <cstring>
#include <linbox/matrix/sparse-matrix.h>
// #include <linbox/vector/blas-vector.h>
#include <linbox/algorithms/block-wiedemann.h>
// #include <givaro/modular-ruint.h>
#include <linbox/matrix/matrix-domain.h>

// extern template class LinBox::MVProductDomain<Givaro::Modular<uint64_t>>;

using F60 = Givaro::Modular<int64_t, __uint128_t>;
struct Part { std::vector<mpz_class> x;  mpz_class mod; };

static mpz_class inv_mod(const mpz_class& a, const mpz_class& m) {
    mpz_class inv;
    int ok = mpz_invert(inv.get_mpz_t(), a.get_mpz_t(), m.get_mpz_t());
    assert(ok && "prefix and modulus must be coprime");
    return inv;
}

std::vector<std::vector<mpz_class>> buildProductTree(std::vector<mpz_class> level) {
    size_t lvl_size = level.size();
    size_t height = 1 + (lvl_size > 1 ? (32 - __builtin_clz((unsigned)lvl_size - 1)) : 0);
    std::vector<std::vector<mpz_class>> levels;
    levels.reserve(height);

    level.reserve(lvl_size);
    levels.push_back(std::move(level));
    std::vector<mpz_class> next_level;
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

std::vector<mpz_class> batchRemainders(const std::vector<std::vector<mpz_class>>& p_levels,
                                        const std::vector<mpz_class>& X) {
    const mpz_class &Z = p_levels.back()[0];
    auto x_levels = buildProductTree(X);
    size_t h = x_levels.size() - 1;
    const mpz_class& root = x_levels[h][0];

    std::vector<mpz_class> rems, next_rems;
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

std::vector<mpz_class> smoothCandidates(const std::vector<std::vector<mpz_class>>& p_levels, const std::vector<mpz_class>& X, const mpz_class& P, const mpz_class& base) {
    std::vector<mpz_class> rems = batchRemainders(p_levels, X);
    
    auto max_it = std::max_element(X.begin(), X.end());
    mpz_class max_x = *max_it;

    unsigned long bits = mpz_sizeinbase(max_x.get_mpz_t(), 2);
    mpz_class M;
    mpz_ui_pow_ui(M.get_mpz_t(), 2, bits);

    std::vector<mpz_class> smooth_cands;

    mpz_class y, g, x_mpz;
    for (size_t i = 0; i < X.size(); ++i) {
        const mpz_class& t_mpz = X[i];
        mpz_powm(x_mpz.get_mpz_t(), base.get_mpz_t(), t_mpz.get_mpz_t(), P.get_mpz_t());
        const mpz_class& r = rems[i];
        mpz_powm(y.get_mpz_t(), r.get_mpz_t(), M.get_mpz_t(), x_mpz.get_mpz_t());
        g = gcd(x_mpz, y);

        if (g == x_mpz) smooth_cands.emplace_back(t_mpz);

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

std::vector<std::pair<size_t, uint32_t>> treeFactorize(const std::vector<std::vector<mpz_class>>& p_levels, const mpz_class &d_mp) {
    std::vector<std::pair<size_t, uint32_t>> result;
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


void linSolve(const std::vector<std::vector<std::pair<size_t,uint32_t>>>& M_rows, const std::vector<mpz_class>& X_col,
                const mpz_class& prime_mpz, std::vector<mpz_class>& L_out) {

    if (!(mpz_sizeinbase(prime_mpz.get_mpz_t(), 2) < 60)) throw std::invalid_argument("linSolve: modulus too large");

    const size_t m = M_rows.size();
    const size_t k = L_out.size();

    int64_t q60 = prime_mpz.get_ui();
    F60 F(q60);

    LinBox::SparseMatrix<F60> M(F, m, k);
    LinBox::DenseVector<F60>  X(F, m), L(F, k);
    LinBox::MatrixDomain<F60> MD(F);
    LinBox::BlockWiedemannSolver<decltype(MD)> solver(MD);

    for (size_t i = 0; i < m; ++i) {
        for (auto const& [j,eij] : M_rows[i]) {
            typename F60::Element e;
            F.init(e, static_cast<int64_t>(eij));
            M.setEntry(i, j, e);
        }
        int64_t xi = mpz_fdiv_ui(X_col[i].get_mpz_t(), q60);
        X[i] = F.init(xi);
    }

    solver.solveNonSingular(L, M, X);

    for (size_t j = 0; j < k; ++j) {
        int64_t v = L[j];
        L_out[j] = mpz_class(static_cast<unsigned long>(v));
    }
}

inline mpz_class dot_row(const std::vector<std::pair<size_t,uint32_t>>& row, const std::vector<mpz_class>& vec) {
    mpz_class sum = 0;
    for (auto [col, coeff32] : row)
        sum += mpz_class(coeff32) * vec[col];
    return sum;
}

inline void hensel_update(std::vector<mpz_class>& x, const std::vector<mpz_class>& delta, const mpz_class& q_pow) {
    const std::size_t n = x.size();
    for (std::size_t j = 0; j < n; ++j)
        x[j] += q_pow * delta[j];
}

inline void garner_step(std::vector<mpz_class>& X, std::vector<mpz_class>& prefix, const std::vector<mpz_class>& xi, const mpz_class& mi) {
    const mpz_class inv = inv_mod(prefix[0], mi);
    const std::size_t n = X.size();
    for (std::size_t j = 0; j < n; ++j) {
        mpz_class t = ( (xi[j] - X[j]) * inv ) % mi;
        X[j] += prefix[j] * t;
        prefix[j] *= mi;
    }
}

std::vector<mpz_class> crtSolve(const std::vector<std::vector<std::pair<size_t,uint32_t>>>& M_rows, const std::vector<mpz_class>& X_col,
                                const std::vector<std::pair<mpz_class,uint32_t>>& factorList ) {
    const std::size_t k = M_rows.front().size();
    std::vector<Part> partial;
    partial.reserve(factorList.size());

    for (const auto& [q_mpz, e] : factorList) {
        const mpz_class q = q_mpz;
        std::vector<mpz_class> x_q(k);
        linSolve(M_rows, X_col, q, x_q);

        mpz_class q_pow = q;
        if (e > 1) {
            std::vector<mpz_class> Ax(M_rows.size());
            for (uint32_t pow = 1; pow < e; ++pow) {
                for (std::size_t row = 0; row < M_rows.size(); ++row) {
                    mpz_class dot = dot_row(M_rows[row], x_q);
                    mpz_class tmp = (X_col[row] - dot) / q_pow;
                    Ax[row] = tmp % q;
                }
                std::vector<mpz_class> delta(k);
                linSolve(M_rows, Ax, q, delta);

                hensel_update(x_q, delta, q_pow);

                q_pow *= q;
            }
        } else q_pow = q;

        partial.push_back({ std::move(x_q),  q_pow });
    }

    const mpz_class modulus = std::accumulate(partial.begin(), partial.end(), mpz_class(1), [](mpz_class a,const Part& p){ return a * p.mod; });

    std::vector<mpz_class> L(k, 0);
    std::vector<mpz_class> prefix(k, 1);

    for (const Part& part : partial) {
        garner_step(L, prefix, part.x, part.mod);
    }
    for (auto& v : L) v %= modulus;
    return L;
}
