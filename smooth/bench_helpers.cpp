
#include <gmpxx.h>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>
#include "gauss_dream.h"
#include "smooth_algos.h"
#include "infra.h"
#include <chrono>
#include <iostream>
#include <cassert>
#include <random>

thread_local gmp_randstate_t rs;
thread_local bool rs_init = []{
    gmp_randinit_mt(rs);
    std::random_device rd;
    gmp_randseed_ui(rs, rd());
    return true;
}();

using SparseList = std::vector<std::pair<size_t, uint32_t>>;
using RelationMatrix = std::vector<SparseList>;
using MpzVector = std::vector<mpz_class>;
using FactorList = std::vector<std::pair<mpz_class, uint32_t>>;


mpz_class findGenerator(const mpz_class& p, FactorList factorization) {
    mpz_class a = 2;
    mpz_class b = p - 1;
    mpz_class range = b - a + 1;
    mpz_class g;
    mpz_class temp;
    mpz_class exp;

    while (true) {
        mpz_urandomm(g.get_mpz_t(), rs, range.get_mpz_t());
        g += a;
        bool valid = true;
        for (const auto& [q_mpz, e] : factorization) {
            exp = (p - 1) / q_mpz;
            mpz_powm (temp.get_mpz_t(), g.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
            if (temp == 1) {
                valid = false;
                break;
            }
        }

        if (valid) return g;
    }
}

inline uint32_t y_bound(const mpz_class& p, double C) {
    double L = mp_ln(p);
    double bound = C * std::exp(std::sqrt(0.5 * L * std::log(L)));
    return static_cast<uint32_t>(bound);
}

inline void verifyFactorization(SparseList factorization, MpzVector factorBase, mpz_class n) {
    mpz_class product = 1;
    mpz_class prime_power;
    for (auto const& [factor, exp] : factorization) {
        mpz_pow_ui(prime_power.get_mpz_t(), factorBase[factor].get_mpz_t(), static_cast<unsigned long>(exp));
        product *= prime_power;
    }
    assert((product == n) && "incorrect factorization found");
}

std::pair<RelationMatrix, MpzVector> naiveRelations(uint32_t k, const mpz_class& p, const mpz_class& g, 
    uint32_t y, uint32_t expectation, MpzVector factorBase, std::vector<MpzVector> levels) {
    RelationMatrix M;
    MpzVector X;
    M.reserve(k);
    X.reserve(k);

    mpz_class t_prime, b_hat;
    mpz_class range = p - 1;

    while (M.size() != k) {
        MpzVector candidates;
        candidates.reserve(expectation);
        for (size_t _ = 0; _ < expectation; ++_) {
            mpz_urandomm(t_prime.get_mpz_t(), rs, range.get_mpz_t());
            candidates.emplace_back(t_prime);
        }

        auto smooth = smoothCandidates(levels, candidates, p, g);
        for (const mpz_class& t : smooth) {
            if (M.size() >= k) break;
            mpz_powm(b_hat.get_mpz_t(), g.get_mpz_t(), t.get_mpz_t(), p.get_mpz_t());
            assert(isSmooth(b_hat, y) && "incorrect smoothness found");
            auto factorization = treeFactorize(levels, b_hat);
            verifyFactorization(factorization, factorBase, b_hat);
            X.emplace_back(t);
            M.emplace_back(factorization);
        }
    }
    assert((M.size() == k && X.size() == k) && "wrong sizes");
    return std::make_pair(M, X);
}

double verifySolve(RelationMatrix M, MpzVector X, const mpz_class& p, const mpz_class& g, uint32_t y, 
    MpzVector factorBase, FactorList factorization) {
    
    auto t0 = std::chrono::high_resolution_clock::now();
    auto L = crtSolve(M, X, factorization);
    auto t1 = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(t1 - t0).count();

    for (size_t row = 0; row < M.size(); ++row) {
        mpz_class sum = 0;
        for (auto [col, coeff] : M[row])
            sum += mpz_class(coeff) * L[col];
        assert(( (sum - X[row]) % (p - 1) == 0 ) && "algebraic check failed");
    }

    for (size_t i = 0; i < L.size(); ++i) {
        mpz_class lhs;
        mpz_powm(lhs.get_mpz_t(), g.get_mpz_t(),
                L[i].get_mpz_t(), p.get_mpz_t());
        assert((lhs == factorBase[i]) && "discrete log check failed");
    }
    return duration;
}

inline std::ostream& operator<<(std::ostream& os, const std::pair<mpz_class, uint32_t>& p) {
    return os << "(" << p.first << ", " << p.second << ")";
}

inline std::ostream& operator<<(std::ostream& os, const std::vector<std::pair<mpz_class, uint32_t>>& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i + 1 != vec.size()) os << ", ";
    }
    return os << "]";
}

double benchSolve(const mpz_class& p, double C, double S) {
    std::cout << "Benchmark p=" << p << " | C=" << C << ", S=" << S << std::endl;

    FactorList factorization = factorize(p - 1);
    std::cout << "Factorized order to " << factorization << std::endl;

    mpz_class g = findGenerator(p, factorization);
    std::cout << "Found generator g=" << g << std::endl;

    uint32_t y = y_bound(p, C);
    double u = mp_ln(p) / std::log(y);
    double rho = std::exp(logDickman(u));
    uint32_t expectation = static_cast<uint32_t>(1.0 / rho);

    MpzVector factorBase = sieveTo(y);
    uint32_t k = factorBase.size();
    std::cout << "Using heuristics: y=" << y << ", expectation=" << expectation << ", k=" << k << std::endl;

    std::vector<MpzVector> levels = buildProductTree(factorBase);
    std::cout << "Built factor base tree" << std::endl;

    auto [M, X] = naiveRelations(k, p, g, y, expectation, factorBase, levels);
    std::cout << "Found relations" << std::endl;

    double duration = verifySolve(M, X, p, g, y, factorBase, factorization);
    return duration;
}

