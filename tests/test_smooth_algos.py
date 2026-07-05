import math

import pytest

import smooth
from conftest import random_bitlength


def test_is_smooth_known_cases():
    assert smooth.is_smooth(77, 11)
    assert not smooth.is_smooth(77, 7)
    assert smooth.is_smooth(1, 5)
    assert smooth.is_smooth(97, 97)
    assert not smooth.is_smooth(97, 50)

    product = 2 * 3 * 5 * 7 * 11 * 13
    assert smooth.is_smooth(product, 13)
    assert not smooth.is_smooth(product, 11)


_SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def test_is_smooth_randomized_constructed_products(rng):
    for _ in range(50):
        bound = rng.choice(_SMALL_PRIMES)
        at_or_below = [p for p in _SMALL_PRIMES if p <= bound]
        above = [p for p in _SMALL_PRIMES if p > bound]

        x = 1
        for _ in range(1 + rng.randrange(5)):
            x *= rng.choice(at_or_below)
        assert smooth.is_smooth(x, bound)

        if above:
            assert not smooth.is_smooth(x * rng.choice(above), bound)


def test_log_dickman_known_values():
    assert math.isclose(smooth.log_dickman(0.0), 0.0, abs_tol=1e-12)
    assert math.isclose(smooth.log_dickman(0.5), 0.0, abs_tol=1e-12)
    assert math.isclose(smooth.log_dickman(1.0), 0.0, abs_tol=1e-12)

    expected = math.log(1.0 - math.log(2.0))
    assert math.isclose(smooth.log_dickman(2.0), expected, abs_tol=1e-3)


def test_log_dickman_is_monotonically_nonincreasing(rng):
    for _ in range(200):
        u1, u2 = sorted(rng.uniform(0.0, 40.0) for _ in range(2))
        if u2 - u1 < 1e-6:
            continue
        assert smooth.log_dickman(u1) >= smooth.log_dickman(u2) - 1e-9


@pytest.mark.parametrize("k", [1, 10, 64, 100, 127])
def test_mp_ln_known_powers_of_two(k):
    assert math.isclose(smooth.mp_ln(2**k), k * math.log(2.0), abs_tol=1e-9)


def test_mp_ln_randomized_against_python_reference(rng):
    for _ in range(50):
        x = random_bitlength(rng, 20 + rng.randrange(100))
        bit_len = x.bit_length()
        shift = max(0, bit_len - 53)
        mantissa = x >> shift
        expected = math.log(mantissa) + shift * math.log(2.0)
        assert math.isclose(smooth.mp_ln(x), expected, rel_tol=1e-9)


def test_log_mul_identity_at_zero():
    for x in [1, 2, 97, 2**60, 2**100 + 3]:
        assert smooth.log_mul(x, 0.0) == x


def test_log_mul_randomized_against_reference(rng):
    for _ in range(50):
        x = random_bitlength(rng, 20 + rng.randrange(100))
        log_rho = rng.uniform(-20.0, 0.0)
        expected = x * math.exp(log_rho)
        got = smooth.log_mul(x, log_rho)
        if expected > 1.0:
            assert abs(got - expected) / expected < 0.01


def test_psi_approx_never_exceeds_x(rng):
    for _ in range(50):
        x = random_bitlength(rng, 40 + rng.randrange(80))
        y = 4 + rng.randrange(1_000_000)
        assert smooth.psi_approx(x, y) <= x


def test_psi_approx_nondecreasing_in_y(rng):
    for _ in range(30):
        x = random_bitlength(rng, 40 + rng.randrange(80))
        y1 = 4 + rng.randrange(100_000)
        y2 = y1 + 1 + rng.randrange(100_000)
        assert smooth.psi_approx(x, y1) <= smooth.psi_approx(x, y2)


def test_psi_approx_matches_x_times_rho_within_tolerance(rng):
    for _ in range(30):
        x = random_bitlength(rng, 40 + rng.randrange(80))
        y = 1000 + rng.randrange(1_000_000)

        u = smooth.mp_ln(x) / math.log(y)
        expected = x * math.exp(smooth.log_dickman(u))
        got = smooth.psi_approx(x, y)

        if expected > 1.0:
            assert abs(got - expected) / expected < 0.01
