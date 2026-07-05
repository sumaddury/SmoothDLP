#pragma once
#include <vector>
#include <utility>
#include <cstdint>
#include "types.h"

extern std::vector<uint32_t> full_primes_array;

std::vector<uint32_t> sieveTo(uint64_t n);

bool isPrime(u128 n);

std::pair<std::vector<std::pair<u128, uint32_t>>, u128> factorize_naive(u128 n);

u128 squfof(u128 n);

std::vector<std::pair<u128, uint32_t>> factorize(u128 n);
