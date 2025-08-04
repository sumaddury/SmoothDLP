
#pragma once
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

std::pair<RelationMatrix, MpzVector> naiveRelations(uint32_t k, const mpz_class& p, const mpz_class& g, 
    uint32_t y, uint32_t expectation, std::vector<MpzVector> levels, double S) {
    RelationMatrix M;
    MpzVector X;
    M.reserve(static_cast<size_t>(S) * k);
    X.reserve(static_cast<size_t>(S) * k);

    mpz_class range = p - 1;
    std::cout << "Using S=" << S << '\n';

    while (X.size() < S * k) {
        MpzVector b_cands;
        MpzVector t_cands;
        mpz_class t_prime, b_hat;

        b_cands.reserve(expectation);
        t_cands.reserve(expectation);
        for (size_t _ = 0; _ < expectation; ++_) {
            mpz_urandomm(t_prime.get_mpz_t(), rs, range.get_mpz_t());
            mpz_powm(b_hat.get_mpz_t(), g.get_mpz_t(), t_prime.get_mpz_t(), p.get_mpz_t());
            b_cands.emplace_back(std::move(b_hat));
            t_cands.emplace_back(std::move(t_prime));
        }
        
        auto smooth_idx = smoothCandidates(levels, b_cands);

        for (size_t idx : smooth_idx) {
            t_prime = t_cands[idx];
            b_hat = b_cands[idx];
            auto factorization = treeFactorize(levels, b_hat);
            std::cout << "row #" << M.size() << '\n';
            X.emplace_back(std::move(t_prime));
            M.emplace_back(factorization);
        }

    }
    return std::make_pair(M, X);
}

double verifySolve(RelationMatrix M, MpzVector X, const mpz_class& p, const mpz_class& g, uint32_t y, 
    MpzVector factorBase, FactorList factorization, size_t k) {
    
    auto t0 = std::chrono::high_resolution_clock::now();
    auto L = crtSolve(M, X, factorization, k);
    std::cout << "Finished solve" << std::endl;
    auto t1 = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(t1 - t0).count();

    // for (size_t row = 0; row < M.size(); ++row) {
    //     mpz_class sum = 0;
    //     for (auto [col, coeff] : M[row])
    //         sum += mpz_class(coeff) * L[col];
    //     assert(( (sum - X[row]) % (p - 1) == 0 ) && "algebraic check failed");
    // }

    // std::cout << "Finished alg check" << std::endl;

    // for (size_t i = 0; i < L.size(); ++i) {
    //     mpz_class lhs;
    //     mpz_powm(lhs.get_mpz_t(), g.get_mpz_t(),
    //             L[i].get_mpz_t(), p.get_mpz_t());
    //     assert((lhs == factorBase[i]) && "discrete log check failed");
    // }
    return duration;
}

inline std::ostream& operator<<(std::ostream& os, const std::pair<mpz_class, uint32_t>& p) {
    return os << "(" << p.first << ", " << p.second << ")";
}

inline std::ostream& operator<<(std::ostream& os, const FactorList& vec) {
    os << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i + 1 != vec.size()) os << ", ";
    }
    return os << "]";
}

inline std::ostream& operator<<(std::ostream& os,
                                const std::pair<std::size_t, std::uint32_t>& p)
{
    return os << "(" << p.first << ", " << p.second << ")";
}

inline std::ostream& operator<<(std::ostream& os, const SparseList& row)
{
    os << "[";
    for (std::size_t i = 0; i < row.size(); ++i) {
        os << row[i];
        if (i + 1 != row.size()) os << ", ";
    }
    return os << "]";
}

inline std::ostream& operator<<(std::ostream& os, const RelationMatrix& M)
{
    os << "[\n";
    for (std::size_t r = 0; r < M.size(); ++r) {
        os << "  " << r << ": " << M[r];
        if (r + 1 != M.size()) os << ",";
        os << '\n';
    }
    return os << "]";
}

inline std::ostream& operator<<(std::ostream& os, const MpzVector& vec)
{
    os << "[";
    for (std::size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i + 1 != vec.size()) os << ", ";
    }
    return os << "]";
}

inline void printRelationMatrix(const RelationMatrix& M)
{
    for (std::size_t row = 0; row < M.size(); ++row) {
        std::cout << row << ':';
        for (const auto& [col, val] : M[row])
            std::cout << " (" << col << ", " << val << ')';
        std::cout << '\n';
    }
}

double benchSolve(const mpz_class& p, double C, double S) {
    std::cout << "Benchmark p=" << p << " | C=" << C << std::endl;

    FactorList factorization = factorize(p - 1);
    std::cout << "Factorized order to " << factorization << std::endl;

    mpz_class g = findGenerator(p, factorization);
    std::cout << "Found generator g=" << g << std::endl;

    uint32_t y = y_bound(p, C);
    double u = mp_ln(p) / std::log(y);
    double rho = std::exp(logDickman(u));
    uint32_t expectation = static_cast<uint32_t>(1.0 / rho);

    auto it = std::upper_bound(full_primes_array.begin(), full_primes_array.end(), y);
    MpzVector factorBase;
    factorBase.reserve(std::distance(full_primes_array.begin(), it));
    std::transform(full_primes_array.begin(), it, std::back_inserter(factorBase),
                [](uint32_t val) { return mpz_class(val); });
    uint32_t k = factorBase.size();
    std::cout << "Using heuristics: y=" << y << ", expectation=" << expectation << ", k=" << k << std::endl;

    std::vector<MpzVector> levels = buildProductTree(factorBase);
    std::cout << "Built factor base tree" << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    auto [M, X] = naiveRelations(k, p, g, y, expectation, levels, S);
    auto t1 = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(t1 - t0).count();
    std::cout << "Duration: " << duration << std::endl;
    std::cout << "M= " << M << std::endl;
    std::cout << "X= " << X << std::endl;

    // double duration = verifySolve(M, X, p, g, y, factorBase, factorization);
    return 1.0;
}

