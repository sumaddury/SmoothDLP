import pytest

import smooth
from conftest import random_bitlength, reference_is_prime_small, reference_miller_rabin


def test_sieve_to_known_small_case():
    assert smooth.sieve_to(30) == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]


@pytest.mark.parametrize("n", [10, 97, 1000, 50_000, 200_000])
def test_sieve_to_matches_independent_reference(n):
    is_composite = [False] * (n + 1)
    expected = []
    for i in range(2, n + 1):
        if is_composite[i]:
            continue
        expected.append(i)
        for j in range(i * i, n + 1, i):
            is_composite[j] = True
    assert smooth.sieve_to(n) == expected


@pytest.mark.parametrize("p", [2, 3, 5, 7, 997, 7919, 104729])
def test_is_prime_known_primes(p):
    assert smooth.is_prime(p)


@pytest.mark.parametrize("c", [0, 1, 4, 6, 100, 999, 1_000_000])
def test_is_prime_known_composites(c):
    assert not smooth.is_prime(c)


def test_is_prime_exhaustive_small_range_vs_reference():
    mismatches = [n for n in range(200_000) if smooth.is_prime(n) != reference_is_prime_small(n)]
    assert mismatches == []


def test_is_prime_known_large_mersenne_primes():
    assert smooth.is_prime(2**61 - 1)
    assert smooth.is_prime(2**127 - 1)


def test_is_prime_known_large_composite():
    assert not smooth.is_prime((2**61 - 1) * 3)


@pytest.mark.parametrize("bits", [40, 64, 96, 127])
def test_is_prime_randomized_vs_oracle(rng, bits):
    for _ in range(100):
        n = random_bitlength(rng, bits)
        assert smooth.is_prime(n) == reference_miller_rabin(n, rng)


def test_is_prime_rejects_out_of_range():
    with pytest.raises(OverflowError):
        smooth.is_prime(2**128)
    with pytest.raises(OverflowError):
        smooth.is_prime(-5)


def _reconstruct(factors):
    n = 1
    for p, e in factors:
        n *= p**e
    return n


def test_factorize_known_small_composite():
    factors = sorted(smooth.factorize(360))
    assert factors == [(2, 3), (3, 2), (5, 1)]


def test_factorize_known_large_semiprime():
    p, q = 1_099_511_628_211, 2_199_023_255_579
    assert smooth.is_prime(p) and smooth.is_prime(q)
    n = p * q
    factors = smooth.factorize(n)
    assert _reconstruct(factors) == n
    assert sorted(f for f, _ in factors) == sorted([p, q])
    assert all(e == 1 for _, e in factors)


def _next_probable_prime(n: int, rng) -> int:
    n |= 1
    while not reference_miller_rabin(n, rng):
        n += 2
    return n


def test_factorize_randomized_reconstruction_and_primality(rng):
    for _ in range(30):
        n = 1
        for _ in range(2 + rng.randrange(3)):
            candidate = random_bitlength(rng, 12 + rng.randrange(10))
            n *= _next_probable_prime(candidate, rng)
        if n >= 2**127:
            continue

        factors = smooth.factorize(n)
        assert _reconstruct(factors) == n
        for f, _ in factors:
            assert reference_miller_rabin(f, rng)


