#pragma once
#include <gmpxx.h>
#include <cstdint>
#include <map>
#include <vector>

std::vector<std::vector<mpz_class>> buildProductTree(std::vector<mpz_class> level);

std::vector<mpz_class> batchRemainders(const std::vector<std::vector<mpz_class>>& p_levels, const std::vector<mpz_class>& X);

std::vector<mpz_class> smoothCandidates(const std::vector<std::vector<mpz_class>>& p_levels, const std::vector<mpz_class>& X);

std::vector<std::pair<size_t, uint32_t>> treeFactorize(const std::vector<std::vector<mpz_class>>& p_levels, const mpz_class &d_mp);


