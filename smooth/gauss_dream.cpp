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
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <utility>
#include <map>
#include <deque>

static std::vector<uint32_t> small_primes_array;
static std::unordered_set<uint32_t> small_primes_set;
std::vector<uint32_t> full_primes_array;
static constexpr uint32_t MAX_SMALL_PRIME = 1'000'000;
static const mpz_class MAX_NAIVE("1000000000000");
static constexpr std::array<uint32_t, 32> M = {{
    1, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131
}};

struct SmallPrimesLoader {
    SmallPrimesLoader() {
        const char* src   = __FILE__;
        const char* slash = std::strrchr(src, '/');
        std::string dir   = slash ? std::string(src, slash - src) : ".";

        std::ifstream in(dir + "/primes_deltas.bin", std::ios::binary);
        if (!in) throw std::runtime_error("cannot open primes_deltas.bin");
        in.seekg(0, std::ios::end);
        std::streamsize N = in.tellg();
        in.seekg(0, std::ios::beg);

        std::vector<uint8_t> deltas;
        deltas.resize(size_t(N));
        if (!in.read(reinterpret_cast<char*>(deltas.data()), N))
            throw std::runtime_error("error reading primes_deltas.bin");

        full_primes_array.reserve(deltas.size());
        uint32_t p = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            p = (i == 0 ? deltas[0] : p + deltas[i]);
            full_primes_array.push_back(p);
        }

        auto it = std::upper_bound(
            full_primes_array.begin(),
            full_primes_array.end(),
            MAX_SMALL_PRIME
        );
        small_primes_array.assign(full_primes_array.begin(), it);
        small_primes_set.insert(
            small_primes_array.begin(),
            small_primes_array.end()
        );
    }
} _small_primes_loader;

std::vector<mpz_class> sieveTo(const mpz_class &n_mp) {
    if (n_mp > std::numeric_limits<size_t>::max()) {
        throw std::overflow_error("sieve bound too large");
    }
    size_t n = n_mp.get_ui();
    std::vector<bool> primality(n+1, true);
    primality[0] = primality[1] = false;

    mpz_class L_mp;
    mpz_sqrt(L_mp.get_mpz_t(), n_mp.get_mpz_t());
    size_t L = static_cast<size_t>(L_mp.get_ui());
    for (size_t i = 2; i <= L; ++i) {
        if (!primality[i]) continue;
        for (size_t j = i*i; j <= n; j+=i) {
            primality[j] = false;
        }
    }
    std::vector<mpz_class> primes;
    double approx = double(n) / std::log(double(n));
    primes.reserve(static_cast<size_t>(approx));
    for (size_t i = 2; i <= n; ++i) {
        if (primality[i]) primes.emplace_back(mpz_class(i));
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

std::pair<std::vector<std::pair<mpz_class, uint32_t>>, mpz_class> factorize_naive(const mpz_class &n_mp) {
    std::vector<std::pair<mpz_class, uint32_t>> output;
    output.reserve(mpz_sizeinbase(n_mp.get_mpz_t(), 2));
    using idx_t = decltype(small_primes_array)::size_type;
    idx_t start = 0, end = 0;
    mpz_class L_mp; 
    size_t    L;
    if (n_mp <= MAX_NAIVE) {
        mpz_sqrt(L_mp.get_mpz_t(), n_mp.get_mpz_t());
        L = static_cast<size_t>(L_mp.get_ui());
        auto it = std::upper_bound(
            small_primes_array.begin(),
            small_primes_array.end(),
            L
        );
        end = it - small_primes_array.begin();
    } else end = static_cast<idx_t>(small_primes_array.size());

    mpz_class n = n_mp;
    bool BREAK_COND = true;
    while (!isPrime(n)) {
        BREAK_COND = true;
        for (idx_t i = start; i < end; ++i) {
            uint32_t p = small_primes_array[i];
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
                mpz_sqrt(L_mp.get_mpz_t(), n.get_mpz_t());
                L = static_cast<size_t>(L_mp.get_ui());
                auto it = std::upper_bound(
                    small_primes_array.begin()+ start,
                    small_primes_array.begin()+ end,
                    L);
                end = it - small_primes_array.begin();
            } else end = static_cast<idx_t>(small_primes_array.size());
            BREAK_COND = false;
            break;
        }
        if (BREAK_COND) return std::make_pair(output, n);
    }
    output.emplace_back(n, 1);
    return std::make_pair(output, 1);
}

mpz_class squfof(const mpz_class &n_mp) {
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
            std::swap(a,c);
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
                    if (h != 1) return h;
                }
            }

        }
    }
    return 1;
}

std::vector<std::pair<mpz_class, uint32_t>> factorize(const mpz_class &n_mp) {
    std::map<mpz_class, uint32_t> FACTORS;
    auto trial = factorize_naive(n_mp);
    auto &small = trial.first;
    mpz_class rem = trial.second;
    if (rem == 1) {
        return small;
    }
    for (auto &pr : small) {
        FACTORS[ pr.first ] += pr.second;
    }
    std::deque<mpz_class> stack;
    stack.push_back(rem);

    mpz_class m;
    mpz_class f;
    mpz_class cr;
    mpz_class crem;
    while (!stack.empty()) {
        m = stack.back();
        stack.pop_back();

        if (isPrime(m)) {
            FACTORS[ m ] += 1;
            continue;
        }
        if (mpz_perfect_square_p(m.get_mpz_t()) != 0) {
            mpz_class r;
            mpz_sqrt(r.get_mpz_t(), m.get_mpz_t());
            stack.push_back(r);
            stack.push_back(r);
            continue;
        }
        mpz_rootrem(cr.get_mpz_t(), crem.get_mpz_t(), m.get_mpz_t(), 3);
        if (crem == 0) {
            stack.push_back(cr);
            stack.push_back(cr);
            stack.push_back(cr);
            continue;
        }

        f = squfof(m);
        stack.push_back(f);
        stack.push_back(m / f);
    }

    std::vector<std::pair<mpz_class,uint32_t>> items{FACTORS.begin(), FACTORS.end()};
    return items;
}