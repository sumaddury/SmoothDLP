#include <gmpxx.h>
#include "gauss_dream.h"
#include "montgomery.h"
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <utility>
#include <map>
#include <deque>
#include <bitset>
#include <random>

namespace gauss {

std::vector<uint32_t> full_primes_array;
static constexpr uint32_t SIEVE_BOUND = 5'000'000;
static constexpr uint32_t MAX_SMALL_PRIME = 1'000'000;
static std::bitset<SIEVE_BOUND + 1> prime_flag;
static constexpr u128 MAX_NAIVE = (u128)1'000'000'000'000ULL;
static constexpr std::array<uint32_t, 32> M = {{
    1, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131
}};
static std::size_t CAP_END;

static constexpr u128 DETERMINISTIC_MR_BOUND = (((u128)0x2be69ULL) << 64) | (u128)0x51adc5b22410a5fdULL;
static constexpr std::array<uint64_t, 7> SMALL_WITNESSES = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
static constexpr int RANDOM_MR_ROUNDS = 20;

static thread_local std::mt19937_64 witness_rng{std::random_device{}()};

struct SmallPrimesLoader {
    SmallPrimesLoader() {
        full_primes_array = sieveTo(MAX_SMALL_PRIME);

        CAP_END = std::upper_bound(
            full_primes_array.begin(),
            full_primes_array.end(),
            MAX_SMALL_PRIME
        ) - full_primes_array.begin();

        prime_flag.reset();
        for (uint32_t q : full_primes_array) {
            prime_flag.set(q);
        }
    }
} _small_primes_loader;

std::vector<uint32_t> sieveTo(uint64_t n) {
    std::vector<bool> composite(n + 1, false);
    uint64_t L = (uint64_t)std::sqrt((double)n);
    while ((L + 1) * (L + 1) <= n) ++L;
    while (L > 0 && L * L > n) --L;

    for (uint64_t i = 2; i <= L; ++i) {
        if (composite[i]) continue;
        for (uint64_t j = i * i; j <= n; j += i) composite[j] = true;
    }

    std::vector<uint32_t> primes;
    if (n >= 2) {
        double approx = double(n) / std::log(double(n > 2 ? n : 3));
        primes.reserve(static_cast<size_t>(approx));
    }
    for (uint64_t i = 2; i <= n; ++i) {
        if (!composite[i]) primes.push_back((uint32_t)i);
    }
    return primes;
}

namespace {

u128 random_u128() {
    uint64_t hi = witness_rng(), lo = witness_rng();
    return (((u128)hi) << 64) | (u128)lo;
}

bool miller_rabin_pass(u128 n, u128 d, unsigned s, u128 a,
                       u128 n_prime, u128 r2, u128 one_bar, u128 nm1_bar) {
    a %= n;
    if (a == 0) return true;

    u128 base_bar = mont::mulmod(a, r2, n, n_prime);
    u128 x_bar = one_bar;
    u128 e = d;
    while (e > 0) {
        if (e & 1) x_bar = mont::mulmod(x_bar, base_bar, n, n_prime);
        base_bar = mont::mulmod(base_bar, base_bar, n, n_prime);
        e >>= 1;
    }

    if (x_bar == one_bar || x_bar == nm1_bar) return true;
    for (unsigned i = 1; i < s; ++i) {
        x_bar = mont::mulmod(x_bar, x_bar, n, n_prime);
        if (x_bar == nm1_bar) return true;
    }
    return false;
}

u128 isqrt(u128 n) {
    if (n == 0) return 0;

    u128 x = ((u128)1) << (((127 - clz128(n)) / 2) + 1);

    while (true) {
        u128 xNext = (x + (n / x)) / 2;
        if (xNext >= x) return x;
        x = xNext;
    }
}

bool cube_le(u128 r, u128 n) {
    if (r == 0) return true;
    return r <= ((n / r) / r);
}

u128 icbrt(u128 n) {
    if (n == 0) return 0;

    u128 x = ((u128)1) << (((127 - clz128(n)) / 3) + 1);

    while (true) {
        u128 xNext = ((2 * x + ((n / x) / x))) / 3;
        if (xNext >= x) break;
        x = xNext;
    }

    if (!cube_le(x, n)) x--;
    else if (cube_le(x + 1, n)) x++;

    return x;
}

} // namespace

bool isPrime(u128 n) {
    if (n <= MAX_SMALL_PRIME) return prime_flag.test((uint32_t)n);
    if (n % 2 == 0) return false;

    u128 nm1 = n - 1;
    u128 d = nm1;
    unsigned s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }

    u128 n_prime = mont::inverse(n);
    u128 r2 = mont::r_squared_mod_n(n);
    u128 one_bar = mont::mulmod(1, r2, n, n_prime);
    u128 nm1_bar = mont::mulmod(nm1, r2, n, n_prime);

    if (n < DETERMINISTIC_MR_BOUND) {
        for (uint64_t a : SMALL_WITNESSES) {
            if ((u128)a >= n) continue;
            if (!miller_rabin_pass(n, d, s, a, n_prime, r2, one_bar, nm1_bar)) return false;
        }
    } else {
        u128 range = n - 3;
        for (int i = 0; i < RANDOM_MR_ROUNDS; ++i) {
            u128 a = 2 + (random_u128() % range);
            if (!miller_rabin_pass(n, d, s, a, n_prime, r2, one_bar, nm1_bar)) return false;
        }
    }
    return true;
}

std::pair<FactorList, u128> factorize_naive(u128 n) {
    FactorList output;
    output.reserve(128 - clz128(n == 0 ? 1 : n));
    using idx_t = std::vector<uint32_t>::size_type;
    idx_t start = 0, end = 0;
    u128 L;

    if (n <= MAX_NAIVE) {
        L = isqrt(n);
        end = std::upper_bound(
                  full_primes_array.begin(),
                  full_primes_array.begin() + CAP_END,
                  (uint32_t)L) - full_primes_array.begin();
    } else {
        end = static_cast<idx_t>(CAP_END);
    }

    bool break_cond = true;
    while (!isPrime(n)) {
        break_cond = true;
        for (idx_t i = start; i < end; ++i) {
            uint32_t p = full_primes_array[i];
            if (n % p != 0) continue;
            uint32_t e = 0;
            while (n % p == 0) {
                n /= p;
                e += 1;
            }
            output.emplace_back(p, e);
            if (n == 1) return std::make_pair(output, 1);
            start = i + 1;
            if (n <= MAX_NAIVE) {
                L = isqrt(n);
                end = std::upper_bound(
                          full_primes_array.begin() + start,
                          full_primes_array.begin() + CAP_END,
                          (uint32_t)L) - full_primes_array.begin();
            } else {
                end = static_cast<idx_t>(CAP_END);
            }
            break_cond = false;
            break;
        }
        if (break_cond) return std::make_pair(output, n);
    }
    output.emplace_back(n, 1);
    return std::make_pair(output, 1);
}

u128 squfof(u128 n_u128) {
    mpz_class n_mp = u128_to_mpz(n_u128);
    mpz_class h, mn, r, rn, b, a, c, temp, L, q, t, v, w, u;
    for (uint32_t m : M) {
        mn = m * n_mp;
        mpz_sqrt(r.get_mpz_t(), mn.get_mpz_t());
        rn = r;
        b = r;
        a = 1;
        h = rn;
        c = mn - h * h;

        temp = 2 * r;
        mpz_sqrt(L.get_mpz_t(), temp.get_mpz_t());
        uint64_t Lbound = 4 * L.get_ui();
        for (uint64_t i = 2; i <= Lbound; ++i) {
            std::swap(a, c);
            q = (rn + b) / a;
            t = b;
            b = q * a - b;
            c += q * (t - b);

            if (!(i & 1)) {
                mpz_sqrt(r.get_mpz_t(), c.get_mpz_t());
                if (r * r == c) {
                    q = (rn - b) / r;
                    v = q * r + b;
                    w = (mn - v * v) / r;
                    u = r;
                    while (true) {
                        std::swap(w, u);
                        r = v;
                        q = (rn + v) / u;
                        v = q * u - v;
                        if (v == r) break;
                        w += q * (r - v);
                    }
                    mpz_gcd(h.get_mpz_t(), n_mp.get_mpz_t(), u.get_mpz_t());
                    if (h != 1) return mpz_to_u128(h);
                }
            }
        }
    }
    return 1;
}

FactorList factorize(u128 n) {
    std::map<u128, uint32_t> factors;
    auto trial = factorize_naive(n);
    auto& small = trial.first;
    u128 rem = trial.second;
    if (rem == 1) return small;
    for (auto& pr : small) factors[pr.first] += pr.second;

    std::deque<u128> stack;
    stack.push_back(rem);

    while (!stack.empty()) {
        u128 m = stack.back();
        stack.pop_back();

        if (isPrime(m)) {
            factors[m] += 1;
            continue;
        }

        u128 r = isqrt(m);
        if (r * r == m) {
            stack.push_back(r);
            stack.push_back(r);
            continue;
        }

        u128 cr = icbrt(m);
        if (cr * cr * cr == m) {
            stack.push_back(cr);
            stack.push_back(cr);
            stack.push_back(cr);
            continue;
        }

        u128 f = squfof(m);
        stack.push_back(f);
        stack.push_back(m / f);
    }

    return FactorList(factors.begin(), factors.end());
}

} // namespace gauss
