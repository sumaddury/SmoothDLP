#include <gmpxx.h>
#include <smooth_algos.h>
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