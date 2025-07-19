#pragma once
#include <gmpxx.h>
#include <cstdint>
#include <map>
#include <vector>


std::vector<std::vector<mpz_class>> buildProductTree(std::vector<mpz_class> level);

std::vector<mpz_class> smoothCandidates(const std::vector<std::vector<mpz_class>>& p_levels, const std::vector<mpz_class>& X, const mpz_class& P, const mpz_class& base);

std::vector<std::pair<size_t, uint32_t>> treeFactorize(const std::vector<std::vector<mpz_class>>& p_levels, const mpz_class &d_mp);

std::vector<mpz_class> crtSolve(const std::vector<std::vector<std::pair<size_t,uint32_t>>>& M_rows, const std::vector<mpz_class>& X_col,
                                const std::vector<std::pair<mpz_class,uint32_t>>& factorList);


