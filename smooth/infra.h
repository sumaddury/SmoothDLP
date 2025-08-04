#pragma once
#include <gmpxx.h>
#include <cstdint>
#include <map>
#include <vector>

using SparseList = std::vector<std::pair<size_t, uint32_t>>;
using RelationMatrix = std::vector<SparseList>;
using MpzVector = std::vector<mpz_class>;
using FactorList = std::vector<std::pair<mpz_class, uint32_t>>;

std::vector<MpzVector> buildProductTree(MpzVector level);

std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const MpzVector& X);

SparseList treeFactorize(const std::vector<MpzVector>& p_levels, const mpz_class &d_mp);

std::size_t rank_relation_gf2(const RelationMatrix& rows, std::size_t k);

MpzVector crtSolve(RelationMatrix& M_rows, const MpzVector& X_col, const FactorList& factorList, const std::size_t k);


