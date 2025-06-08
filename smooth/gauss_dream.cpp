#include <gmpxx.h>
#include "gauss_dream.h"
#include <vector>
#include <array>
#include <cmath>
#include <unordered_set>
#include <fstream>
#include <cstdint>
#include <iostream>
#include <cstring>


static std::vector<uint32_t> small_primes_array;
static std::unordered_set<uint32_t> small_primes_set;
static constexpr uint32_t MAX_SMALL_PRIME = 1'000'000;

struct SmallPrimesLoader {
  SmallPrimesLoader() {
    const char* src   = __FILE__;
    const char* slash = std::strrchr(src, '/');
    std::string dir   = slash ? std::string(src, slash - src) : ".";

    std::ifstream in(dir + "/small_primes.bin", std::ios::binary);
    uint32_t p;
    while (in.read(reinterpret_cast<char*>(&p), sizeof(p))) {
      small_primes_array.push_back(p);
    }
    small_primes_set.insert(
      small_primes_array.begin(), small_primes_array.end());
  }
} _loader;

std::vector<int> sieveTo(int n) {
    std::vector<bool> primality(n+1, true);
    primality[0] = primality[1] = false;

    int L = static_cast<int>(std::sqrt(n));
    for (int i = 2; i <= L; ++i) {
        if (!primality[i]) {
            continue;
        }
        for (int j = i*i; j <= n; j+=i) {
            primality[j] = false;
        }
    }
    std::vector<int> primes;
    primes.reserve(n / std::log(n));
    for (int i = 2; i <= n; ++i) {
        if (primality[i]) primes.push_back(i);
    }
    return primes;
}

bool isPrime(const mpz_class &n_mp) {
    if (n_mp <= MAX_SMALL_PRIME) {
        uint32_t u = static_cast<uint32_t>(n_mp.get_ui());
        return small_primes_set.count(u) > 0;
    }
    if (mpz_perfect_power_p(n_mp.get_mpz_t()) != 0) {
        return false;
    }
    
    mpz_class nm1 = n_mp - 1;
    unsigned s = mpz_scan1(nm1.get_mpz_t(), 0);
    mpz_class t = nm1 >> s;

    auto mr_body = [&](uint64_t a) -> bool {
        mpz_class g;
        mpz_gcd_ui(g.get_mpz_t(), n_mp.get_mpz_t(), a);
        if (g > 1) return false;

        mpz_class a_mp(static_cast<unsigned long>(a)), x;
        mpz_powm(x.get_mpz_t(), a_mp.get_mpz_t(), t.get_mpz_t(), n_mp.get_mpz_t());
        if (x == 1 || x == nm1) {
            return true;
        }
        for (unsigned i = 1; i < s; ++i) {
            x *= x; x %= n_mp;
            if (x == nm1) {
                return true;
            }
        }
        return false;
    };
    
    if (mpz_sizeinbase(n_mp.get_mpz_t(), 2) <= 64) {
        constexpr std::array<uint64_t,7> small_wits = {2,325,9375,28178,450775,9780504,1795265022};
        unsigned long n_ul = n_mp.get_ui();
        for (uint64_t a : small_wits) {
            if (a >= n_ul) continue;
            if (!mr_body(a)) return false;
        }
    } else {
        double ln_n = std::log(mpz_get_d(n_mp.get_mpz_t()));
        uint64_t M = static_cast<uint64_t>(std::floor(2 * ln_n * ln_n)) + 1;
        for (uint64_t a = 2; a < M; ++a) {
            if (!mr_body(a)) return false;
        }
    }
    return true;
}
