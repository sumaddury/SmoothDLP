#pragma once
#include "types.h"

namespace mont {

inline void mul128(u128 a, u128 b, u128& hi, u128& lo) {
    hi = 0;
    lo = 0;
}

inline u128 mulmod(u128 a_bar, u128 b_bar, u128 n, u128 n_prime) {
    return 0;
}

inline u128 inverse(u128 n) {
    return 0;
}

inline u128 reduce256(u128 hi, u128 lo, u128 n) {
    return 0;
}

inline u128 r_squared_mod_n(u128 n) {
    return 0;
}

inline u128 gcd(u128 a, u128 b) {
    return 0;
}

inline u128 pow2mod_odd(u128 base, unsigned bits, u128 n) {
    return 0;
}

inline u128 mulmod_any(u128 a, u128 b, u128 n) {
    return 0;
}

inline u128 pow2mod_any(u128 base, unsigned bits, u128 n) {
    return 0;
}

inline u128 pow2mod(u128 base, unsigned bits, u128 n) {
    return 0;
}

} // namespace mont
