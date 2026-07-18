#include <gmpxx.h>
#include "smooth_algos.h"
#include "gauss_dream.h"
#include "infra.h"
#include "montgomery.h"
#include <random>
#include <cstdint>
#include <cmath>
#include <limits>
#include <map>
#include <vector>
#include <utility>
#include <unordered_set>
#include <cassert>
#include <cstring>
#include <numeric>

namespace infra {

// The LinBox/Givaro-based sparse linear solve below (linSolve/linSolveImpl/
// crtSolve/rank_relation_gf2, and the includes/typedefs they need) is
// disabled: we haven't actually figured out how to use Givaro/LinBox
// correctly yet, so no line of code touching either library should be
// trusted. Nothing else in this file depends on any of it -- buildProductTree
// /smoothCandidates/treeFactorize only ever use GMP + montgomery.h. Not
// declared in infra.h, not bound in core.cpp, not covered by any test.
#if 0
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

using Field64 = Givaro::Modular<int64_t, __uint128_t>;
using Field32 = Givaro::Modular<uint32_t, uint64_t>;
template<typename Field> using Sparse = LinBox::SparseMatrix<Field>;
template<typename Field> using Vec = LinBox::DenseVector<Field>;
struct Part { MpzVector x; mpz_class mod; };
using Seq = LinBox::SparseMatrixFormat::SparseSeq;
using Gf2 = LinBox::GF2;
#endif // LinBox/Givaro disabled

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

        const auto& cur = levels.back();
        for (size_t i = 0; i < lvl_size; i += 2) {
            if (i + 1 < lvl_size) {
                next_level.emplace_back(cur[i] * cur[i + 1]);
            } else next_level.emplace_back(cur[lvl_size - 1]);
        }

        levels.push_back(std::move(next_level));
        lvl_size = levels.back().size();
    }
    return levels;
}

static MpzVector batchRemainders(const std::vector<MpzVector>& p_levels, const MpzVector& X) {
    const mpz_class& Z = p_levels.back()[0];
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

std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const U128Vector& X) {
    MpzVector X_mp;
    X_mp.reserve(X.size());
    for (u128 x : X) X_mp.push_back(u128_to_mpz(x));

    MpzVector rems = batchRemainders(p_levels, X_mp);

    u128 max_x = *std::max_element(X.begin(), X.end());
    unsigned bits = 128 - clz128(max_x);

    std::vector<size_t> smooth_cands;
    for (size_t i = 0; i < X.size(); ++i) {
        u128 x = X[i];
        u128 r = mpz_to_u128(rems[i]);
        // pow2mod's m_prime/r2_m are relative to x's odd part (not x
        // itself), computed via the montgomery interfaces. Each candidate
        // here is a distinct modulus, called once, so there's nothing to
        // reuse across iterations -- just compute and pass through.
        u128 m = x >> ctz128(x);
        u128 m_prime = mont::inverse(m);
        u128 r2_m = mont::r_squared_mod_n(m);
        u128 y = mont::pow2mod(r, bits, x, m_prime, r2_m);
        u128 g = mont::gcd(x, y);
        if (g == x) smooth_cands.push_back(i);
    }
    return smooth_cands;
}

namespace {
inline uint32_t find_pow(mpz_class d, const mpz_class& P) {
    uint32_t exp = 0;
    do {
        exp++;
        d /= P;
    } while (d % P == 0);
    return exp;
}
} // namespace

SparseList treeFactorize(const std::vector<MpzVector>& p_levels, u128 d_u128) {
    SparseList result;
    result.reserve(12);
    std::vector<std::pair<size_t, size_t>> stack;
    stack.reserve(p_levels.size());
    stack.emplace_back(p_levels.size() - 1, 0);

    mpz_class g;
    mpz_class d = u128_to_mpz(d_u128);
    size_t level, idx, right, left;
    uint32_t exp;
    while (!stack.empty()) {
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

namespace {

// Uniform random value in [0, L] via rejection sampling: draw a full 128-bit
// value, mask it down to L's bit range (mask is the caller-supplied all-ones
// bitmask covering [0, L]), and retry on the rare draw that still lands
// above L. Below 2^64, a single uniform_int_distribution draw suffices.
u128 drawT(u128 L, u128 mask) {
  thread_local std::random_device rd;
  thread_local std::mt19937_64 gen(rd());

  if (L <= std::numeric_limits<uint64_t>::max())
    return std::uniform_int_distribution<uint64_t>(0, (uint64_t)L)(gen);

  u128 num;

  do {
    uint64_t lo = std::uniform_int_distribution<uint64_t>(
        0, std::numeric_limits<uint64_t>::max())(gen);
    uint64_t hi = std::uniform_int_distribution<uint64_t>(
        0, std::numeric_limits<uint64_t>::max())(gen);
    num = (((u128)hi << 64) | (u128)lo) & mask;
  } while (num > L);

  return num;
}

} // namespace

void addRelations(
    u128 base,
    const ProblemParams& params,
    std::vector<SparseList>& M,
    U128Vector& X) {
  const size_t k = params.p_levels[0].size();
  size_t added = 0;

  while (added < k) {
    size_t batch_size = (size_t)std::ceil((double)(k - added) / params.smooth_density);

    U128Vector batch, t_batch;
    batch.reserve(batch_size);
    t_batch.reserve(batch_size);

    for (size_t iter = 0; iter < batch_size; iter++) {
      u128 t = drawT(params.p - 2, params.mask);
      u128 x = mont::powmod_odd(base, t, params.p, params.p_prime, params.r2);
      batch.push_back(x);
      t_batch.push_back(t);
    }

    auto smooth_x = smoothCandidates(params.p_levels, batch);

    for (size_t idx : smooth_x) {
      M.push_back(treeFactorize(params.p_levels, batch[idx]));
      X.push_back(t_batch[idx]);
    }

    added += smooth_x.size();
  }
}

} // namespace infra

// Disabled along with the LinBox/Givaro includes/typedefs above -- see the
// comment there. linSolveImpl/linSolve/dot_row/inv_mod/hensel_update/
// garner_step/rank_relation_gf2/crtSolve are the entire untrusted
// LinBox-based sparse solve + CRT/Hensel-lift pipeline; none of it is
// exercised by anything outside this disabled block.
#if 0



#endif // LinBox/Givaro disabled
