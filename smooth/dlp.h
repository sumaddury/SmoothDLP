#include <gmpxx.h>
#include <random>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>
#include <unordered_map>
#include <string>

std::unordered_map<string, uint32_t> optimize(const mpz_class& P_mp, int T);