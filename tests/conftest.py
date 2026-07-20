import random

import pytest

smooth = pytest.importorskip("smooth")

_REQUIRED_NEW_API = [
    "build_product_tree",
    "smooth_candidates",
    "tree_factorize",
    "log_mul",
    "DLP",
]


def pytest_configure():
    missing = [name for name in _REQUIRED_NEW_API if not hasattr(smooth, name)]
    if missing:
        pytest.exit(f"smooth._core is missing {missing}", returncode=1)


@pytest.fixture
def rng():
    return random.Random(20260705)


def reference_is_prime_small(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def reference_gf2_rank(rows: list[list[int]], k: int) -> int:
    matrix = []
    for r in rows:
        row = [False] * k
        for c in r:
            row[c] = True
        matrix.append(row)

    rank = 0
    for col in range(k):
        pivot = None
        for r in range(rank, len(matrix)):
            if matrix[r][col]:
                pivot = r
                break
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for r in range(len(matrix)):
            if r != rank and matrix[r][col]:
                matrix[r] = [a ^ b for a, b in zip(matrix[r], matrix[rank])]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def random_bitlength(rng: random.Random, bits: int) -> int:
    v = rng.getrandbits(bits - 1) | (1 << (bits - 1))
    return v | 1


def reference_miller_rabin(n: int, rng: random.Random, rounds: int = 40) -> bool:
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p

    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1

    def is_composite_witness(a: int) -> bool:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            return False
        for _ in range(r - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                return False
        return True

    if n < (1 << 64):
        witnesses = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    else:
        witnesses = [rng.randrange(2, n - 1) for _ in range(rounds)]

    return not any(is_composite_witness(a % n) for a in witnesses if a % n not in (0, 1))
