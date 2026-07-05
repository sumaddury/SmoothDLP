import math

import pytest

import smooth
from conftest import reference_gf2_rank, reference_miller_rabin

_SMALL_FACTOR_BASE = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
_OUTSIDE_FACTOR_BASE = [53, 59, 61, 67, 71, 73, 79]


def test_build_product_tree_known_small_case():
    levels = smooth.build_product_tree([2, 3, 5, 7])
    assert levels[-1] == [210]


def test_build_product_tree_root_matches_independent_product(rng):
    for _ in range(20):
        values = [rng.getrandbits(10 + rng.randrange(20)) | 1 for _ in range(1 + rng.randrange(40))]
        expected = math.prod(values)
        levels = smooth.build_product_tree(values)
        assert levels[-1][0] == expected


def test_smooth_candidates_known_smooth_and_nonsmooth():
    levels = smooth.build_product_tree(_SMALL_FACTOR_BASE)
    X = [
        2 * 3 * 5 * 7 * 11,
        2 * 3 * 53,
        41 * 43 * 47,
        59 * 61,
    ]
    assert sorted(smooth.smooth_candidates(levels, X)) == [0, 2]


def test_smooth_candidates_randomized_constructed_batch(rng):
    levels = smooth.build_product_tree(_SMALL_FACTOR_BASE)

    X, expected_smooth = [], []
    for _ in range(60):
        make_smooth = rng.random() < 0.5
        x = 1
        for _ in range(1 + rng.randrange(4)):
            x *= rng.choice(_SMALL_FACTOR_BASE)
        if not make_smooth:
            x *= rng.choice(_OUTSIDE_FACTOR_BASE)
        X.append(x)
        expected_smooth.append(make_smooth)

    got_idx = set(smooth.smooth_candidates(levels, X))
    for i, expect in enumerate(expected_smooth):
        assert (i in got_idx) == expect


def test_tree_factorize_known_composite():
    levels = smooth.build_product_tree(_SMALL_FACTOR_BASE)
    d = 8 * 3 * 11
    result = smooth.tree_factorize(levels, d)

    reconstructed = 1
    for idx, exp in result:
        reconstructed *= _SMALL_FACTOR_BASE[idx] ** exp
    assert reconstructed == d


def test_tree_factorize_randomized_reconstruction(rng):
    levels = smooth.build_product_tree(_SMALL_FACTOR_BASE)
    for _ in range(40):
        d = 1
        for _ in range(1 + rng.randrange(5)):
            d *= rng.choice(_SMALL_FACTOR_BASE)

        result = smooth.tree_factorize(levels, d)
        reconstructed = 1
        for idx, exp in result:
            reconstructed *= _SMALL_FACTOR_BASE[idx] ** exp
        assert reconstructed == d
