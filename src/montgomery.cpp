#include "montgomery.h"

namespace mont {

void mul128(u128 a, u128 b, u128& hi, u128& lo) {
	const u128 a0 = (uint64_t)a;
	const u128 a1 = (a >> 64);
	const u128 b0 = (uint64_t)b;
	const u128 b1 = (b >> 64);

	const u128 p00 = a0 * b0;
	const u128 p01 = a0 * b1;
	const u128 p10 = a1 * b0;
	const u128 p11 = a1 * b1;

	const u128 r0 = (uint64_t)p00;
	const u128 m1 = (p00 >> 64) + ((uint64_t)p01) + ((uint64_t)p10);
	const u128 r1 = (uint64_t)m1;
	const u128 m2 = (p01 >> 64) + (p10 >> 64) + ((uint64_t)p11) + (m1 >> 64);
	const u128 r2 = (uint64_t)m2;
	const u128 r3 = (p11 >> 64) + (m2 >> 64);
	
	lo = ((r1 << 64) | r0);
	hi = ((r3 << 64) | r2);   
}

u128 gcd(u128 a, u128 b) {
	if (a == 0) return b;
	if (b == 0) return a;

	const int shift = ctz128(a | b);
	a >>= shift;
	b >>= shift;

	a >>= ctz128(a);
	while (b != 0) {
		b >>= ctz128(b);
		if (a > b) std::swap(a, b);
		b -= a;
	}
	
	return (a << shift);
}

u128 mulmod_any(u128 a, u128 b, u128 n) {
	u128 u_lo, u_hi;
	mul128(a, b, u_hi, u_lo);

	return reduce256(u_hi, u_lo, n);
}

u128 inverse(u128 n) {
	u128 x = n;

	for (int step = 0; step < 7; step++) x = x * ((u128)2 - n * x);

	return x;
}

u128 r_squared_mod_n(u128 n) {
	u128 r_reduced = (((u128)(-1) % n) + 1) % n;
	
	return mulmod_any(r_reduced, r_reduced, n);
}

u128 reduce256(u128 hi, u128 lo, u128 n) {
	std::array<uint64_t, 2> v_limbs{
		(uint64_t)n,
		(uint64_t)(n >> 64),
	};

	if (v_limbs[1] == 0) {
        uint64_t n_lo = v_limbs[0];
        uint64_t limbs[4] = {
            (uint64_t)(hi >> 64), (uint64_t)hi,
            (uint64_t)(lo >> 64), (uint64_t)lo
        };
        uint64_t r = 0;
        for (int i = 0; i < 4; i++) {
            u128 num = ((u128)r << 64) | limbs[i];
            r = (uint64_t)(num % n_lo);
        }
        return (u128)r;
    }

	const int d =  __builtin_clzll(v_limbs[1]);

	if (d) {
		v_limbs[1] <<= d;
		v_limbs[1] |= (v_limbs[0] >> (64 - d));
		v_limbs[0] <<= d;
	}

	std::array<uint64_t, 5> u_limbs{
		(uint64_t)lo,
		(uint64_t)(lo >> 64),
		(uint64_t)(hi),
		(uint64_t)(hi >> 64),
	};

	if (d) {
		for (int i = 3; i >= 0; i--) {
			u_limbs[i + 1] |= (u_limbs[i] >> (64 - d));
			u_limbs[i] <<= d;
		}
	}	
	
	std::array<uint64_t, 3> w_limbs{
		u_limbs[2], 
		u_limbs[3], 
		u_limbs[4],
	};

	for (int step = 1; step >= -1; step--) {
		uint64_t q_hat;
		u128 r_hat;
		u128 num = (((u128)w_limbs[2] << 64) | w_limbs[1]);
		
		if (w_limbs[2] >= v_limbs[1]) q_hat = (uint64_t)(-1);
		else q_hat = num / v_limbs[1];
		
		r_hat = num - (u128)q_hat * v_limbs[1];

		while (r_hat < ((u128)1 << 64) && (u128)q_hat * v_limbs[0] > (((u128)r_hat << 64) | w_limbs[0])) {
			q_hat--;
			r_hat += v_limbs[1];
		}
		
		std::array<uint64_t, 3> qV_limbs;
		uint64_t carry = 0;

		for (int i = 0; i < 3; i++) {
			u128 prod = (u128)q_hat * (i < 2 ? v_limbs[i] : 0) + carry;
			qV_limbs[i] = (uint64_t)prod;
			carry = (uint64_t)(prod >> 64);
		}

		uint64_t borrow = 0;
		
		for (int i = 0; i < 3; i++) {
			__int128 diff = (__int128)w_limbs[i] - (__int128)qV_limbs[i] - borrow;
			
			if (diff < 0) {
				diff += (__int128)1 << 64;
				borrow = 1;
			} else borrow = 0;
			
			w_limbs[i] = (uint64_t)diff;
		}

		if (borrow) {
            u128 sum = (u128)w_limbs[0] + v_limbs[0];
            w_limbs[0] = (uint64_t)sum;
            uint64_t c2 = (uint64_t)(sum >> 64);
            w_limbs[1] = (uint64_t)((u128)w_limbs[1] + v_limbs[1] + c2);
        }

		w_limbs[2] = w_limbs[1];
		w_limbs[1] = w_limbs[0];
		if (step >= 0) w_limbs[0] = u_limbs[step];
		
	}

	u128 output = ((u128)w_limbs[2] << 64) | w_limbs[1];
	if (d) output >>= d;

	return output;
}

u128 mulmod(u128 a_bar, u128 b_bar, u128 n, u128 n_prime) {

}

} // namespace mont