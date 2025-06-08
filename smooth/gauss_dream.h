#pragma once
#include <vector>
#include <string>
#include <gmpxx.h>

std::vector<int> sieveTo(int n);

bool isPrime(const mpz_class &n_mp);
