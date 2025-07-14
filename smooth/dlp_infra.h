#pragma once
#include <gmpxx.h>
#include <cstdint>
#include <map>
#include <vector>

std::vector<std::vector<mpz_class>> buildProductTree(std::vector<mpz_class> level);

std::vector<mpz_class> batchRemainders(const vector<vector<mpz_class>>& p_levels,
                                        const vector<mpz_class>& X);

std::vector<mpz_class> smoothCandidates(const vector<vector<mpz_class>>& p_levels, const vector<mpz_class>& X);

std::vector<std::pair<const mpz_class*, uint32_t>> treeFactorize(const vector<vector<mpz_class>>& p_levels, const mpz_class &d_mp)


