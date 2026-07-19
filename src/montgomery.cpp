#include "montgomery.h"

namespace mont {

void mul128(u128 a, u128 b, u128& hi, u128& lo) {
	const u128 a0 = static_cast<uint64_t>(a);
	const u128 a1 = (a >> 64);
	const u128 b0 = static_cast<uint64_t>(b);
	const u128 b1 = (b >> 64);

	const u128 p00 = a0 * b0;
	const u128 p01 = a0 * b1;
	const u128 p10 = a1 * b0;
	const u128 p11 = a1 * b1;

	const u128 r0 = static_cast<uint64_t>(p00);
	const u128 m1 = (p00 >> 64) + static_cast<uint64_t>(p01) + static_cast<uint64_t>(p10);
	const u128 r1 = static_cast<uint64_t>(m1);
	const u128 m2 = (p01 >> 64) + (p10 >> 64) + static_cast<uint64_t>(p11) + (m1 >> 64);
	const u128 r2 = static_cast<uint64_t>(m2);
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

	for (int step = 0; step < 7; step++) x = x * (static_cast<u128>(2) - n * x);

	return x;
}

u128 r_squared_mod_n(u128 n) {
	const u128 r_reduced = ((static_cast<u128>(-1) % n) + 1) % n;
	
	return mulmod_any(r_reduced, r_reduced, n);
}

u128 reduce256(u128 hi, u128 lo, u128 n) {
	std::array<uint64_t, 2> v_limbs{
		static_cast<uint64_t>(n),
		static_cast<uint64_t>(n >> 64),
	};

	if (v_limbs[1] == 0) {
        const uint64_t n_lo = v_limbs[0];
        uint64_t limbs[4] = {
            static_cast<uint64_t>(hi >> 64), static_cast<uint64_t>(hi),
            static_cast<uint64_t>(lo >> 64), static_cast<uint64_t>(lo)
        };
        uint64_t r = 0;
        for (int i = 0; i < 4; i++) {
            u128 num = (static_cast<u128>(r) << 64) | limbs[i];
            r = static_cast<uint64_t>(num % n_lo);
        }
        return static_cast<u128>(r);
    }

	const int d =  __builtin_clzll(v_limbs[1]);

	if (d) {
		v_limbs[1] <<= d;
		v_limbs[1] |= (v_limbs[0] >> (64 - d));
		v_limbs[0] <<= d;
	}

	std::array<uint64_t, 5> u_limbs{
		static_cast<uint64_t>(lo),
		static_cast<uint64_t>(lo >> 64),
		static_cast<uint64_t>(hi),
		static_cast<uint64_t>(hi >> 64),
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
		const u128 num = ((static_cast<u128>(w_limbs[2]) << 64) | w_limbs[1]);

		if (w_limbs[2] >= v_limbs[1]) q_hat = static_cast<uint64_t>(-1);
		else q_hat = num / v_limbs[1];

		r_hat = num - static_cast<u128>(q_hat) * v_limbs[1];

		while (r_hat < (static_cast<u128>(1) << 64) && static_cast<u128>(q_hat) * v_limbs[0] > ((r_hat << 64) | w_limbs[0])) {
			q_hat--;
			r_hat += v_limbs[1];
		}

		std::array<uint64_t, 3> qV_limbs;
		uint64_t carry = 0;

		for (int i = 0; i < 3; i++) {
			u128 prod = static_cast<u128>(q_hat) * (i < 2 ? v_limbs[i] : 0) + carry;
			qV_limbs[i] = static_cast<uint64_t>(prod);
			carry = static_cast<uint64_t>(prod >> 64);
		}

		uint64_t borrow = 0;

		for (int i = 0; i < 3; i++) {
			__int128 diff = static_cast<__int128>(w_limbs[i]) - static_cast<__int128>(qV_limbs[i]) - borrow;

			if (diff < 0) {
				diff += static_cast<__int128>(1) << 64;
				borrow = 1;
			} else borrow = 0;

			w_limbs[i] = static_cast<uint64_t>(diff);
		}

		if (borrow) {
            const u128 sum = static_cast<u128>(w_limbs[0]) + v_limbs[0];
            w_limbs[0] = static_cast<uint64_t>(sum);
            const uint64_t c2 = static_cast<uint64_t>(sum >> 64);
            w_limbs[1] = static_cast<uint64_t>(static_cast<u128>(w_limbs[1]) + v_limbs[1] + c2);
        }

		w_limbs[2] = w_limbs[1];
		w_limbs[1] = w_limbs[0];
		if (step >= 0) w_limbs[0] = u_limbs[step];
		
	}

	u128 output = (static_cast<u128>(w_limbs[2]) << 64) | w_limbs[1];
	if (d) output >>= d;
	
	return output;
	 
}

u128 mulmod(u128 a_bar, u128 b_bar, u128 n, u128 n_prime) {
	u128 T_hi, T_lo;
	mul128(a_bar, b_bar, T_hi, T_lo);

	const u128 m = -( T_lo * n_prime );

	u128 mn_hi, mn_lo;
	mul128(m, n, mn_hi, mn_lo);

	const u128 sum_lo = T_lo + mn_lo;
	const bool carry_lo = (sum_lo < T_lo) ? 1 : 0;

	const u128 s = T_hi + mn_hi;
	bool carry_out = (s < T_hi) ? 1 : 0;
	const u128 t_wrapped = s + carry_lo;
	carry_out |= (t_wrapped < s) ? 1 : 0;

	return (carry_out || t_wrapped >= n) ? t_wrapped - n : t_wrapped;
}

u128 pow2mod_odd(u128 base, unsigned bits, u128 n, u128 n_prime, u128 r2) {

	const u128 b0 = base % n;
	u128 x_bar = mulmod(b0, r2, n, n_prime);

	for (;bits; bits--) x_bar = mulmod(x_bar, x_bar, n, n_prime);

	return mulmod(x_bar, 1, n, n_prime);
}

u128 pow2mod(u128 base, unsigned bits, u128 n, u128 m_prime, u128 r2_m) {
	if (n == 1) return 0;
	
	const int e = ctz128(n);

	const u128 mask = (static_cast<u128>(1) << e) - 1;
	u128 x = base & mask;

	for (int k = 0; k < bits; k++) x = (x * x) & mask;

	const u128 m = (n >> e);

	if (m == 1) return x;

	const u128 r_m = pow2mod_odd(base, bits, m, m_prime, r2_m);
	const u128 m_inv = m_prime & mask;

	const u128 t = ((x - r_m) * m_inv) & mask;
	return r_m + m * t;
	
}

u128 powmod_odd(u128 base, u128 e, u128 n, u128 n_prime, u128 r2) {
	if (n == 1) return 0;

	const u128 b0 = base % n;
	u128 x_base = mulmod(b0, r2, n, n_prime);
	u128 x_bar = mulmod(1, r2, n, n_prime);

	while (e) {
		if (e & 1) x_bar = mulmod(x_bar, x_base, n, n_prime);

		x_base = mulmod(x_base, x_base, n, n_prime);

		e >>= 1;
	}

	return mulmod(x_bar, 1, n, n_prime);
}
	



} // namespace mont
