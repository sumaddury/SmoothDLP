#pragma once
#include <cstdint>
#include "types.h"

extern const uint32_t Y_SMOOTHNESS_BOUND;

bool isSmooth(u128 x, uint32_t y);

u128 log_mul(u128 x, double log_rho);

double mp_ln(u128 x);

double logDickman(double u);

u128 psiApprox(u128 x, uint64_t y);
