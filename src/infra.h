#pragma once
#include <cstdint>
#include <vector>
#include "types.h"

std::vector<MpzVector> buildProductTree(MpzVector level);

std::vector<size_t> smoothCandidates(const std::vector<MpzVector>& p_levels, const U128Vector& X);

SparseList treeFactorize(const std::vector<MpzVector>& p_levels, u128 d);
