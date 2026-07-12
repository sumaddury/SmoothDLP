#pragma once
#include <cstdint>
#include "types.h"

/**
 * Largest smoothness bound isSmooth will accept for y.
 */
extern const uint32_t Y_SMOOTHNESS_BOUND;

/**
 * True iff x is y-smooth, i.e. every prime factor of x is <= y. Trial-
 * divides x by primes up to y (via full_primes_array), short-circuiting
 * early (x <= y, or the remaining cofactor is prime/too large to be smooth)
 * rather than fully factoring x.
 */
bool isSmooth(u128 x, uint32_t y);

/**
 * Computes x * exp(log_rho) for a (possibly very large) integer x, via GMP
 * integer scaling rather than naive floating-point multiplication: splits
 * exp(log_rho) into a fixed-digit decimal mantissa times a power of ten, so
 * precision isn't lost to x overflowing a double's exact-integer range.
 */
u128 log_mul(u128 x, double log_rho);

/**
 * Natural log of a (possibly >64-bit) integer x: keeps the top 53
 * significant bits as a double mantissa and adds back shift * ln(2) for
 * whatever was shifted off, since x itself may not fit in a double exactly.
 */
double mp_ln(u128 x);

/**
 * ln(rho(u)), where rho is Dickman's function -- rho(u) approximates the
 * density of y-smooth numbers near x, for u = ln(x)/ln(y). Piecewise:
 * closed form for u<=2, cubic-spline interpolation over a precomputed table
 * (loaded from dickman_table.bin) for 2<u<20, an asymptotic saddle-point
 * formula (via boost::math::expint) for u>=20. Nonincreasing in u.
 */
double logDickman(double u);

/**
 * Estimates Psi(x, y), the count of y-smooth integers up to x, as
 * x * rho(ln(x)/ln(y)) (see logDickman), computed via log_mul to preserve
 * precision for large x. Never exceeds x; nondecreasing in y.
 */
u128 psiApprox(u128 x, uint64_t y);
