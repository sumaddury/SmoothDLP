#include <gmpxx.h>
#include "smooth_algos.h"
#include "gauss_dream.h"
#include <random>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>
#include <unordered_set>


std::map<int,double> prime_logs;

thread_local gmp_randstate_t rs;
thread_local bool rs_init = []{
    gmp_randinit_mt(rs);
    std::random_device rd;
    gmp_randseed_ui(rs, rd());
    return true;
}();

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

std::vector<mpz_class> batchRemainders(const vector<vector<mpz_class>>& p_levels,
                                        const vector<mpz_class>& X) {
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

std::vector<mpz_class> smoothCandidates(const vector<vector<mpz_class>>& p_levels, const vector<mpz_class>& X) {
    std::vector<mpz_class> rems = batchRemainders(p_levels, X);
    
    auto max_it = std::max_element(X.begin(), X.end());
    mpz_class max_x = *max_it;

    unsigned long bits = mpz_sizeinbase(max_x.get_mpz_t(), 2);
    mpz_class M;
    mpz_ui_pow_ui(M.get_mpz_t(), 2, bits);

    std::vector<mpz_class> smooth_cands;

    mpz_class y, g;
    for (size_t i = 0; i < X.size(); ++i) {
        const mpz_class& x_mpz = X[i];
        const mpz_class& r = rems[i];
        y = powm(r, M, x_mpz);
        g = gcd(x_mpz, y);

        if (g == x_mpz) smooth_cands.emplace_back(x_mpz);

    }
    return smooth_cands;
}

std::vector<std::pair<const mpz_class*, uint32_t>> treeFactorize(const vector<vector<mpz_class>>& p_levels, const mpz_class &d_mp) {
    std::vector<std::pair<const mpz_class*, uint32_t>> result;
    result.reserve(12);
    std::vector<pair<size_t,size_t>> stack;
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
        const mpz_class* Pptr = &p_levels[level][idx];
        g = gcd(d, *Pptr);
        if (g == 1) continue;
        if (level == 0) {
            exp = 0;
            do {
                exp++;
                d /= *Pptr;
            } while (d % *Pptr == 0);
            result.emplace_back(Pptr, exp);
        } else {
            left = 2 * idx;
            right = left + 1;
            stack.push_back(std::make_pair(level - 1, left));
            if (right < p_levels[level - 1].size()) stack.push_back(std::make_pair(level - 1, right));
        }
    }
    return result;
}