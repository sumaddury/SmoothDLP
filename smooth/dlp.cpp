#include <gmpxx.h>
#include <random>
#include <cstdint>
#include <map>
#include <vector>
#include <utility>

std::unordered_map<string, uint32_t> optimize(const mpz_class& P_mp, int T) {

}

"""
Suppose a user gives a P and a T
ignore T for now

P gives a range of y values to choose from 

for each y:
cost(y) = {
    k(y) = y / ln(y),
    solve_cost(k) = O(k^2),
    expected_batch(k, P, y) = k / dickman(ln(P) / ln(y)),
    batch_cost(y, batch_size) = batch_size * interpolate(y, batch_size)
    loss_function(y, expected_batch) = (expected_batch / P) * interpolate(y, expected_batch / P)
}
"""
