#include <gmpxx.h>
#include "smooth_algos.h"
#include "gauss_dream.h"
#include "infra.h"
// Unity-included, not a normal header include: LinBox's commentator/
// args-parser machinery is pulled in transitively by even lin_alg.h's base
// matrix/vector headers, and defines several symbols out-of-line without
// `inline` in this LinBox build -- so any second translation unit that also
// includes lin_alg.h and gets linked into the same binary hits ~60
// duplicate-symbol errors (see lin_alg.cpp's own top-of-file comment, and
// tests/internal/Makefile). Including the .cpp directly keeps this file and
// lin_alg's implementation in one translation unit and sidesteps that.
#include "lin_alg.cpp"
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

namespace {

inline uint32_t find_pow(mpz_class d, const mpz_class& P) {
    uint32_t exp = 0;
    do {
        exp++;
        d /= P;
    } while (d % P == 0);
    return exp;
}

int computeMaskBitlength(u128 p) {
  return 128 - clz128(p - 2);
}

u128 computeMask(int mask_bitlength) {
  return (mask_bitlength == 128) ? ~static_cast<u128>(0) : ((static_cast<u128>(1) << mask_bitlength) - 1);
}

// C = 1/sqrt(2): a fixed heuristic constant, not tuned or tabulated per-p.
u128 computeSmoothnessBound(u128 p) {
  const double C = 1.0 / std::sqrt(2.0);
  const double ln_p = salgo::mp_ln(p);
  const double ln_B = C * std::sqrt(ln_p * std::log(ln_p));
  return static_cast<u128>(std::ceil(std::exp(ln_B)));
}

std::vector<MpzVector> productTreeForBound(u128 B) {
  const std::vector<uint32_t> factor_base = gauss::sieveTo(static_cast<uint32_t>(B));
  return buildProductTree(MpzVector(factor_base.begin(), factor_base.end()));
}

// Uniform random value in [0, L] via rejection sampling: draw a full 128-bit
// value, mask it down to L's bit range (mask is the caller-supplied all-ones
// bitmask covering [0, L]), and retry on the rare draw that still lands
// above L. Below 2^64, a single uniform_int_distribution draw suffices.
u128 drawT(u128 L, u128 mask) {
  thread_local std::random_device rd;
  thread_local std::mt19937_64 gen(rd());

  if (L <= std::numeric_limits<uint64_t>::max())
    return std::uniform_int_distribution<uint64_t>(0, static_cast<uint64_t>(L))(gen);

  u128 num;

  do {
    const uint64_t lo = std::uniform_int_distribution<uint64_t>(
        0, std::numeric_limits<uint64_t>::max())(gen);
    const uint64_t hi = std::uniform_int_distribution<uint64_t>(
        0, std::numeric_limits<uint64_t>::max())(gen);
    num = ((static_cast<u128>(hi) << 64) | static_cast<u128>(lo)) & mask;
  } while (num > L);

  return num;
}

MpzVector batchRemainders(const std::vector<MpzVector>& p_levels, const MpzVector& X) {
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

} // namespace

std::vector<MpzVector> buildProductTree(MpzVector level) {
    size_t lvl_size = level.size();
    size_t height = 1 + (lvl_size > 1 ? (32 - __builtin_clz(static_cast<unsigned>(lvl_size) - 1)) : 0);
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

std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const U128Vector& X) {
    MpzVector X_mp;
    X_mp.reserve(X.size());
    for (u128 x : X) X_mp.push_back(u128_to_mpz(x));

    const MpzVector rems = batchRemainders(p_levels, X_mp);

    const u128 max_x = *std::max_element(X.begin(), X.end());
    const unsigned bits = 128 - clz128(max_x);

    std::vector<size_t> smooth_cands;
    for (size_t i = 0; i < X.size(); ++i) {
        const u128 x = X[i];
        const u128 r = mpz_to_u128(rems[i]);
        // pow2mod's m_prime/r2_m are relative to x's odd part (not x
        // itself), computed via the montgomery interfaces. Each candidate
        // here is a distinct modulus, called once, so there's nothing to
        // reuse across iterations -- just compute and pass through.
        const u128 m = x >> ctz128(x);
        const u128 m_prime = mont::inverse(m);
        const u128 r2_m = mont::r_squared_mod_n(m);
        const u128 y = mont::pow2mod(r, bits, x, m_prime, r2_m);
        const u128 g = mont::gcd(x, y);
        if (g == x) smooth_cands.push_back(i);
    }
    return smooth_cands;
}

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

ProblemParams::ProblemParams(u128 p) :
  p(p),
  p_prime(mont::inverse(p)),
  r2(mont::r_squared_mod_n(p)),
  mask_bitlength(computeMaskBitlength(p)),
  mask(computeMask(mask_bitlength)),
  B(computeSmoothnessBound(p)),
  // logDickman returns ln(rho(u)); addRelations divides batch sizes by
  // smooth_density directly, so this needs the actual probability, not
  // its log.
  smooth_density(std::exp(salgo::logDickman(salgo::mp_ln(p) / salgo::mp_ln(B)))),
  p_levels(productTreeForBound(B)),
  p_factorization(gauss::factorize(p - 1)) {}

constexpr int WIEDEMANN_THRESHOLD = 4000;

void addRelations(
    u128 base,
    const ProblemParams& params,
    RelationMatrix& M,
    U128Vector& X) {
  const size_t k = params.p_levels[0].size();
  size_t added = 0;

  while (added < k) {
    size_t batch_size = static_cast<size_t>(std::ceil(static_cast<double>(k - added) / params.smooth_density));

    U128Vector batch, t_batch;
    batch.reserve(batch_size);
    t_batch.reserve(batch_size);

    for (size_t iter = 0; iter < batch_size; iter++) {
      const u128 t = drawT(params.p - 2, params.mask);
      const u128 x = mont::powmod_odd(base, t, params.p, params.p_prime, params.r2);
      batch.push_back(x);
      t_batch.push_back(t);
    }

    const auto smooth_x = smoothCandidates(params.p_levels, batch);

    for (size_t idx : smooth_x) {
      M.push_back(treeFactorize(params.p_levels, batch[idx]));
      X.push_back(t_batch[idx]);
    }

    added += smooth_x.size();
  }
}

/**
  q^e via binary exponentiation. e is a factor-power exponent (bounded by
  log2 of p-1, so at most ~128 even in the most pathological u128 case, far
  fewer iterations in practice), and the result is one of p-1's own
  prime-power factors, so it never exceeds p-1 -- both comfortably exact in
  u128 with no risk of overflow.
*/
inline u128 binExp(u128 q, int e) {
    u128 base = q;
    u128 accum = 1;

    while (e) {
        if (e & 1) accum *= base;
        base *= base;
        e >>= 1;
    }

    return accum;
}

/**
  Modular inverse of P_prev mod n = q^e, where P_prev is the product of
  p-1's OTHER prime-power factors processed so far in crtSolve's CRT loop
  below -- guaranteed coprime to n, since p-1's distinct prime factors are
  never repeated across factors.

  q == 2 is handled separately: Euler's theorem needs phi(n), which for
  n = 2^e is a power of two itself and gcd(P_prev, n) = 1 is guaranteed the
  same way, but there's a cheaper route. mont::inverse(P_prev) computes
  P_prev^-1 mod 2^128 (valid since P_prev is odd whenever q == 2, being a
  product of the OTHER, necessarily-odd prime-power factors); a 2-adic
  inverse's low e bits are exactly its inverse mod 2^e (this is the same
  fact that makes Newton's-method doubling-precision inverse computation
  work at all), so masking down to e bits is exact, not an approximation.

  For odd q, Euler's theorem applies directly: phi(q^e) = q^(e-1)*(q-1), and
  P_prev^(phi-1) mod n is P_prev^-1 mod n whenever gcd(P_prev, n) = 1.
*/
inline u128 modInv(u128 P_prev, u128 q, int e, u128 n) {
    if (q != 2) {
        const u128 phi = binExp(q, e - 1) * (q - 1);   // Euler's totient of a prime power
        const u128 n_prime = mont::inverse(n);
        const u128 r2 = mont::r_squared_mod_n(n);

        return mont::powmod_odd(P_prev % n, phi - 1, n, n_prime, r2);
    } else {
        return mont::inverse(P_prev) & (((u128)1 << e) - 1);
    }
}

U128Vector crtSolve(
    const ProblemParams& params,
    const RelationMatrix& M,
    const U128Vector& X) {
    std::vector<std::pair<U128Vector, u128>> bases;
    bases.reserve(params.p_factorization.size());

    const size_t n_cols = params.p_levels[0].size();
    const lin_alg::Method method = (
        (n_cols == M.size() && n_cols >= WIEDEMANN_THRESHOLD) ?
        lin_alg::Method::Wiedemann :
        lin_alg::Method::SparseElimination);

    for (auto [q, e] : params.p_factorization) {
        U128Vector L;
        // henselLift can genuinely fail to reach q^e for a rank-deficient
        // base factor (see lin_alg.h's doc on henselLift for why); there is
        // no partial/best-effort result to salvage from that -- the honest
        // answer is "not solvable from this relation set", signaled by
        // returning empty, not by CRT-combining a corrupted L.
        if (lin_alg::henselLift(M, X, n_cols, q, e, L, method) != lin_alg::SolveStatus::OK)
            return {};
        bases.emplace_back(std::move(L), binExp(q, e));
    }

    U128Vector C;
    u128 P;

    std::tie(C, P) = bases[0];

    for (size_t iter = 1; iter < params.p_factorization.size(); ++iter) {
        auto& [L, n] = bases[iter];
        auto [q, e] = params.p_factorization[iter];

        const u128 P_inv = modInv(P, q, e, n);
        U128Vector t(L.size());

        for (size_t i = 0; i < L.size(); ++i) {
            // L[i] and (C[i] mod n) are each < n here, but n itself can be
            // close to p-1's own full size (a single large prime-power
            // factor), so their product with P_inv -- each also < n --
            // can exceed u128 (n^2 can be ~150 bits for a ~75-bit factor).
            // mont::mulmod_any reduces mod n via a full 256-bit product
            // instead of a native u128 multiply, so it doesn't overflow.
            // The subtraction has the same issue lin_alg::henselLift's
            // residual step does: nothing orders L[i] against C[i] mod n,
            // so an explicit wraparound replaces unsigned underflow. Unlike
            // henselLift's residual (where both operands share q^(k+1) as
            // their natural modulus), C[i]'s natural bound here is P, not
            // n, so it must be reduced mod n before the comparison even
            // makes sense.
            const u128 c_mod_n = C[i] % n;
            const u128 diff = (L[i] >= c_mod_n) ? (L[i] - c_mod_n) : (L[i] + n - c_mod_n);
            t[i] = mont::mulmod_any(diff, P_inv, n);
        }

        // P and n are distinct factors of p-1, so P*n is itself a divisor
        // of p-1 and fits comfortably in u128 -- unlike the product above,
        // this one is safe as a native multiply.
        for (size_t i = 0; i < L.size(); ++i) C[i] += P * t[i];

        P *= n;
    }

    return C;
}

} // namespace infra
