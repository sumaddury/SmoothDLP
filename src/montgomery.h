#pragma once
#include "types.h"

namespace mont {

/**
  Full 128x128 -> 256-bit product: (hi:lo) = a * b, exact (no reduction).
*/
void mul128(u128 a, u128 b, u128& hi, u128& lo);

/**
  Montgomery REDC multiply: returns (a_bar * b_bar * R^-1) mod n, where
  R = 2^128. Requires odd n, and a_bar/b_bar already reduced (< n) --
  callers pass genuine Montgomery-form values in and get one back out.
  n_prime must be mont::inverse(n) (the *positive* inverse, n*n_prime == 1
  mod R); mulmod negates it internally to get the REDC constant n' with
  n*n' == -1 mod R, so callers never pre-negate it themselves.
*/
u128 mulmod(u128 a_bar, u128 b_bar, u128 n, u128 n_prime);

/**
  Modular inverse of odd n modulo 2^128: returns x such that n*x == 1
  (mod 2^128), via Newton's iteration (quadratic convergence, doubling
  correct bits each step) starting from x=n.
*/
u128 inverse(u128 n);

/**
  Reduces a 256-bit value (hi:lo) modulo n, for any n > 0: returns
  (hi*2^128 + lo) mod n. A plain full reduction (schoolbook/Knuth long
  division), not Montgomery REDC -- no n_prime involved.
*/
u128 reduce256(u128 hi, u128 lo, u128 n);

/**
  Returns R^2 mod n, where R = 2^128 -- the constant needed to convert an
  ordinary integer into Montgomery form (via one mulmod call with this as
  the other operand).
*/
u128 r_squared_mod_n(u128 n);

/**
  Greatest common divisor of a and b (binary/Stein's algorithm).
*/
u128 gcd(u128 a, u128 b);

/**
  Returns base^(2^bits) mod n for odd n -- i.e. `bits` successive modular
  squarings (the exponent is always exactly a power of two here, so plain
  repeated squaring suffices; no square-and-multiply needed). Computed via
  Montgomery form: convert base in, square `bits` times via mulmod, convert
  back out. n_prime = mont::inverse(n), r2 = mont::r_squared_mod_n(n) --
  both precomputed by the caller (rather than recomputed on every call, an
  expense than can matter when calling repeatedly against the same n).
*/
u128 pow2mod_odd(u128 base, unsigned bits, u128 n, u128 n_prime, u128 r2);

/**
  Returns (a * b) mod n for any n > 0 (odd or even), via mul128 + reduce256
  -- the general-purpose fallback when Montgomery form doesn't apply.
*/
u128 mulmod_any(u128 a, u128 b, u128 n);

/**
  Returns base^(2^bits) mod n for ANY n >= 1 (odd or even) -- same
  quantity as pow2mod_odd, generalized. Splits n = 2^e * m (e = ctz(n), m
  odd), solves the 2^e part with a cheap bitmask (no division), solves the
  m part via pow2mod_odd, and recombines with CRT (Garner's formula). When
  n is odd this degenerates to exactly pow2mod_odd. m_prime/r2_m are
  relative to *m* (n's odd part), not to n itself -- i.e.
  m_prime = mont::inverse(n >> ctz(n)), r2_m = mont::r_squared_mod_n(n >> ctz(n)).
*/
u128 pow2mod(u128 base, unsigned bits, u128 n, u128 m_prime, u128 r2_m);

/**
  Returns base^e mod n for odd n and an arbitrary exponent e (unlike
  pow2mod_odd, e need not be a power of two) via Montgomery-form
  square-and-multiply. n_prime/r2 follow the same hoisting convention as
  pow2mod_odd. n == 1 is special-cased to return 0 directly.
*/
u128 powmod_odd(u128 base, u128 e, u128 n, u128 n_prime, u128 r2);

} // namespace mont
