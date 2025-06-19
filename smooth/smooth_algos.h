#pragma once
#include <utility>
#include <gmpxx.h>

extern const uint32_t Y_SMOOTHNESS_BOUND;

bool isSmooth(const mpz_class &x_mp, uint32_t y);