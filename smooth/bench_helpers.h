#pragma once
#include <gmpxx.h>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>
#include <cassert>

using RelationMatrix = std::vector<std::vector<std::pair<size_t, uint32_t>>>;
using MpzVector = std::vector<mpz_class>;
using FactorList = std::vector<std::pair<mpz_class, uint32_t>>;

mpz_class findGenerator(const mpz_class& p, FactorList factorization);

std::pair<RelationMatrix, MpzVector> naiveRelations(uint32_t k, const mpz_class& p, const mpz_class& g, 
    uint32_t y, uint32_t expectation, MpzVector factorBase, std::vector<MpzVector> levels);

double verifySolve(RelationMatrix M, MpzVector X, const mpz_class& p, const mpz_class& g, uint32_t y, 
    MpzVector factorBase, FactorList factorization);

double benchSolve(const mpz_class& p, double C, double S);