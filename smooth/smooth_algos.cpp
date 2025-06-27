#pragma once
#include <string>
#include <gmpxx.h>
#include <boost/math/special_functions/expint.hpp>
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
constexpr int DIGITS = 9;
constexpr int NEWTON_TRIALS = 5;
const double inv_ln10 = 1.0 / std::log(10.0);
const double scale = std::pow(10.0, DIGITS - 1);
const double LN2 = std::log(2.0);

static std::vector<double> U_LIST;
static std::vector<double> LOG_RHO_LIST;

struct DickmanTableLoader {
    DickmanTableLoader() {
        const char* src     = __FILE__;
        const char* slash   = std::strrchr(src, '/');
        std::string dir     = slash ? std::string(src, slash - src) : ".";

        std::string path = dir + "/dickman_table.bin";
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
            in.read(reinterpret_cast<char*>(&u),  sizeof(u));
            in.read(reinterpret_cast<char*>(&lr), sizeof(lr));
            U_LIST[i]       = u;
            LOG_RHO_LIST[i] = lr;
        }
    }
}; static DickmanTableLoader _dickman_table_loader;

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

mpz_class log_mul(const mpz_class& x_mp, double log_rho) {
    int e10 = static_cast<int>(std::floor(log_rho * inv_ln10));
    double rem = log_rho - e10 * std::log(10.0);
    double mant = std::exp(rem);

    uint32_t M = static_cast<uint32_t>(mant * scale + 0.5);
    int shift = e10 - (DIGITS - 1);

    mpz_class z = x_mp * M;
    if (shift != 0) {
        mpz_class pow10;
        mpz_ui_pow_ui(pow10.get_mpz_t(), 10, static_cast<unsigned>(std::abs(shift)));
        if (shift > 0) {
            z *= pow10;
        } else {
            z /= pow10;
        }
    }
    return z;
}

double logDickman(double u) {
    if (0 <= u && u <= 1) return 0.0;
    if (1 <= u && u <= 2) return std::log(1.0 - std::log(u));
    if (2 < u && u < 20) {
        auto it = std::upper_bound(U_LIST.begin(), U_LIST.end(), u);
        size_t idx = it - U_LIST.begin() - 1;
        double u0 = U_LIST[idx], u1 = U_LIST[idx + 1];
        double y0 = LOG_RHO_LIST[idx], y1 = LOG_RHO_LIST[idx + 1];
        double h = u1 - u0;

        double d0;
        double d1;
        if (idx > 0) {
            double up = U_LIST[idx - 1], yp = LOG_RHO_LIST[idx - 1];
            d0 = (y1 - yp) / (u1 - up);
        } else {
            d0 = (y1 - y0) / h;
        }
        if (idx + 2 < U_LIST.size()) {
            double un = U_LIST[idx + 2], yn = LOG_RHO_LIST[idx + 2];
            d1 = (yn - y0) / (un - u0);
        } else {
            d1 = (y1 - y0) / h;
        }
        
        double t = (u - u0) / h;
        double t2 = t*t;
        double t3 = t2*t;
        double h00 = 2.0*t3 - 3.0*t2 + 1.0;
        double h10 = t3 - 2*t2 + t;
        double h01 = -2*t3 + 3*t2;
        double h11 = t3 - t2;

        return (h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1);
    }

    double xi = std::log(u) + std::log(std::log(u));
    double ex;
    for (int _ = 1; _ <= NEWTON_TRIALS; ++_) {
        ex = std::exp(xi);
        xi -= (ex - 1.0 - u*xi) / (ex - u);
    }
    double Ei_val = boost::math::expint(xi);
    return (Ei_val - u*xi - std::log(xi) - 0.5*std::log(2*(boost::math::constants::pi<double>())*u));
}

double mp_ln(const mpz_class& x_mp) {
    signed long int exp2;
    double m = mpz_get_d_2exp(&exp2, x_mp.get_mpz_t());
    return std::log(m) + exp2 * LN2;
}

mpz_class psiApprox(const mpz_class &x_mp, uint64_t y) {
    double u = mp_ln(x_mp) / std::log(static_cast<double>(y));
    double log_rho = logDickman(u);
    return log_mul(x_mp, log_rho);
}