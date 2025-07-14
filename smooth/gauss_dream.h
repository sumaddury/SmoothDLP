#pragma once
#include <vector>
#include <string>
#include <utility>
#include <gmpxx.h>

extern std::vector<uint32_t> full_primes_array;

std::vector<mpz_class> sieveTo(const mpz_class &n_mp);

bool isPrime(const mpz_class &n_mp);

std::pair<std::vector<std::pair<mpz_class, uint32_t>>, mpz_class> factorize_naive(const mpz_class &n_mp, uint64_t B = MAX_NAIVE, size_t K_END = CAP_END);

mpz_class squfof(const mpz_class &n_mp);

std::vector<std::pair<mpz_class, uint32_t>> factorize(const mpz_class &n_mp);