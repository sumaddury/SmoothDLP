#include <gmpxx.h>
#include "smooth_algos.h"
#include "gauss_dream.h"
#include <boost/math/special_functions/expint.hpp>
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

namespace salgo {

constexpr uint32_t Y_SMOOTHNESS_BOUND = 33'554'432;
constexpr int DIGITS = 9;
constexpr int NEWTON_TRIALS = 5;
const double inv_ln10 = 1.0 / std::log(10.0);
const double scale = std::pow(10.0, DIGITS - 1);
const double LN2 = std::log(2.0);

static std::vector<double> U_LIST;
static std::vector<double> LOG_RHO_LIST;

struct DickmanTableLoader {
    DickmanTableLoader() {
        const char* src = __FILE__;
        const char* slash = std::strrchr(src, '/');
        const std::string dir = slash ? std::string(src, slash - src) : ".";

        const std::string path = dir + "/dickman_table.bin";
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            throw std::runtime_error("cannot open " + path);
        }

        uint32_t n;
        in.read(reinterpret_cast<char*>(&n), sizeof(n));
        if (!in) throw std::runtime_error("error reading entry count");

        U_LIST.resize(n);
        LOG_RHO_LIST.resize(n);

        for (uint32_t i = 0; i < n; ++i) {
            double u, lr;
            in.read(reinterpret_cast<char*>(&u), sizeof(u));
            in.read(reinterpret_cast<char*>(&lr), sizeof(lr));
            U_LIST[i] = u;
            LOG_RHO_LIST[i] = lr;
        }
    }
};
static DickmanTableLoader _dickman_table_loader;

bool isSmooth(u128 x, uint32_t y) {
    if (y > Y_SMOOTHNESS_BOUND) throw std::invalid_argument("isSmooth: y exceeds max bound");
    if (x <= y) return true;
    if (gauss::isPrime(x)) return false;

    u128 rem = x;
    const auto it = std::upper_bound(gauss::full_primes_array.begin(), gauss::full_primes_array.end(), y);
    const size_t idx = it - gauss::full_primes_array.begin();
    for (size_t i = 0; i < idx; ++i) {
        const uint32_t p = gauss::full_primes_array[i];
        if (rem % p != 0) continue;
        do {
            rem /= p;
        } while (rem % p == 0);
        if (rem <= y) return true;
        if (rem < static_cast<u128>(p) * p) return rem == 1;
        if (gauss::isPrime(rem)) return false;
    }
    return rem == 1;
}

u128 log_mul(u128 x, double log_rho) {
    const int e10 = static_cast<int>(std::floor(log_rho * inv_ln10));
    const double rem = log_rho - e10 * std::log(10.0);
    const double mant = std::exp(rem);

    const uint32_t coef = static_cast<uint32_t>(mant * scale + 0.5);
    const int shift = e10 - (DIGITS - 1);

    mpz_class z = u128_to_mpz(x) * coef;
    if (shift != 0) {
        mpz_class pow10;
        mpz_ui_pow_ui(pow10.get_mpz_t(), 10, static_cast<unsigned>(std::abs(shift)));
        if (shift > 0) {
            z *= pow10;
        } else {
            z /= pow10;
        }
    }
    return mpz_to_u128(z);
}

double logDickman(double u) {
    if (0 <= u && u <= 1) return 0.0;
    if (1 <= u && u <= 2) return std::log(1.0 - std::log(u));
    if (2 < u && u < 20) {
        const auto it = std::upper_bound(U_LIST.begin(), U_LIST.end(), u);
        const size_t idx = it - U_LIST.begin() - 1;
        const double u0 = U_LIST[idx], u1 = U_LIST[idx + 1];
        const double y0 = LOG_RHO_LIST[idx], y1 = LOG_RHO_LIST[idx + 1];
        const double h = u1 - u0;

        double d0, d1;
        if (idx > 0) {
            const double up = U_LIST[idx - 1], yp = LOG_RHO_LIST[idx - 1];
            d0 = (y1 - yp) / (u1 - up);
        } else {
            d0 = (y1 - y0) / h;
        }
        if (idx + 2 < U_LIST.size()) {
            const double un = U_LIST[idx + 2], yn = LOG_RHO_LIST[idx + 2];
            d1 = (yn - y0) / (un - u0);
        } else {
            d1 = (y1 - y0) / h;
        }

        const double t = (u - u0) / h;
        const double t2 = t * t;
        const double t3 = t2 * t;
        const double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        const double h10 = t3 - 2 * t2 + t;
        const double h01 = -2 * t3 + 3 * t2;
        const double h11 = t3 - t2;

        return (h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1);
    }

    double xi = std::log(u) + std::log(std::log(u));
    double ex;
    for (int _ = 1; _ <= NEWTON_TRIALS; ++_) {
        ex = std::exp(xi);
        xi -= (ex - 1.0 - u * xi) / (ex - u);
    }
    const double Ei_val = boost::math::expint(xi);
    return (Ei_val - u * xi - std::log(xi) - 0.5 * std::log(2 * (boost::math::constants::pi<double>()) * u));
}

double mp_ln(u128 x) {
    if (x == 0) return -std::numeric_limits<double>::infinity();
    const int bit_len = 128 - clz128(x);
    const int shift = bit_len > 53 ? bit_len - 53 : 0;
    const u128 mantissa_int = shift ? (x >> shift) : x;
    const double mantissa = static_cast<double>(mantissa_int);
    return std::log(mantissa) + shift * LN2;
}

u128 psiApprox(u128 x, uint64_t y) {
    const double u = mp_ln(x) / std::log(static_cast<double>(y));
    const double log_rho = logDickman(u);
    return log_mul(x, log_rho);
}

} // namespace salgo
