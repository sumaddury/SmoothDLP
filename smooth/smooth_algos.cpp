#include <gmpxx.h>
#include "gauss_dream.h"
#include "smooth_algos.h"
#include <vector>
#include <array>
#include <cmath>
#include <unordered_set>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <cstring>
#include <limits>
#include <algorithm>
#include <utility>


constexpr uint32_t Y_SMOOTHNESS_BOUND = 10'000'000;

bool isSmooth(const mpz_class &x_mp, uint32_t y) {
    if (y > Y_SMOOTHNESS_BOUND ) throw std::invalid_argument("isSmooth: y exceeds max bound");
    if (x_mp <= y) return true;
    if (isPrime(x_mp)) return false;

    mpz_class x = x_mp;
    uint32_t p;
    auto it = std::upper_bound(full_primes_array.begin(), full_primes_array.end(), y);
    size_t idx = it - full_primes_array.begin();
    for (size_t i = 0; i < idx; ++i) {
        p = full_primes_array[i];
        if (!mpz_divisible_ui_p(x.get_mpz_t(), p)) continue;
        do {
            mpz_tdiv_q_ui(x.get_mpz_t(), x.get_mpz_t(), p);
        } while (mpz_divisible_ui_p(x.get_mpz_t(), p));
        if (x <= y) return true;
        if (x < p * p) return (x == 1);
        if (isPrime(x)) return false;
    }
    return (x == 1);
}

