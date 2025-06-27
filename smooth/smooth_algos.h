#pragma once
#include <utility>
#include <gmpxx.h>

extern const uint32_t Y_SMOOTHNESS_BOUND;

bool isSmooth(const mpz_class &x_mp, uint32_t y);

mpz_class log_mul(const mpz_class& x_mp, double log_rho);

double mp_ln(const mpz_class& x_mp);

double logDickman(double u);

mpz_class psiApprox(const mpz_class &x_mp, uint64_t y);