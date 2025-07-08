#include <gmpxx.h>
#include "smooth_algos.h"
#include <random>
#include <cstdint>

thread_local gmp_randstate_t rs;
thread_local bool rs_init = []{
    gmp_randinit_mt(rs);
    std::random_device rd;
    gmp_randseed_ui(rs, rd());
    return true;
}();

mpz_class randSmooth(uint32_t  y, const mpz_class& B) {
    if (B <= 1) throw std::invalid_argument("B must be ≥ 2");

    (void)rs_init;

    mpz_class x;
    for (;;) {
        do {
            mpz_urandomm(x.get_mpz_t(), rs, B.get_mpz_t());
        } while (x == 0);

        if (isSmooth(x, y)) return x;
    }
}

int legendreSymbol(const mpz_class& a_mp, const mpz_class& p_mp) {
    mpz_class a = a_mp % p_mp;
    if (a == 0 || a == 1) return a.get_si();
    if (a == p_mp - 1) {
        return (p_mp % 4 == 1) ? 1 : -1;
    }
    if (a == 2) {
        mpz_class mod_temp = p_mp % 8;
        return (mod_temp == 1 || mod_temp == 7) ? 1 : -1;
    }
    mpz_class symbol = mpz_legendre(a.get_mpz_t(), p_mp.get_mpz_t());
    return symbol.get_si();
}

