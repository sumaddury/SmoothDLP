import pytest

import smooth


def _is_primitive_root(candidate: int, p: int) -> bool:
    # smooth.factorize here is test setup (finding valid g/b inputs), not
    # the thing under test -- same role gauss::factorize plays in the C++
    # side's test_dlp.cpp equivalent of this helper.
    order = p - 1
    return all(pow(candidate, order // q, p) != 1 for q, _ in smooth.factorize(order))


def _find_primitive_root(p: int, skip: int) -> int:
    cand = 2
    while cand < p:
        if cand != skip and _is_primitive_root(cand, p):
            return cand
        cand += 1
    raise ValueError(f"no primitive root found below {p}")


def _brute_force_discrete_log(g: int, b: int, p: int) -> int:
    val = 1
    for t in range(p - 1):
        if val == b:
            return t
        val = (val * g) % p
    raise ValueError("g is not a generator of b")


@pytest.mark.parametrize("p", [1009, 8191, 524287])
def test_dlp_solve_matches_brute_force(p):
    g = _find_primitive_root(p, 0)
    b = _find_primitive_root(p, g)

    want = _brute_force_discrete_log(g, b, p)
    got = smooth.DLP(p).solve(g, b)

    assert got == want
    assert pow(g, got, p) == b


def test_dlp_solve_self_consistent_for_large_prime():
    # p ~ 2^31: brute force is no longer practical, but g**x == b (mod p) is
    # still a complete correctness certificate on its own, since g is a
    # primitive root and there's exactly one x in [0, p-2] satisfying it.
    p = 2147483647
    g = _find_primitive_root(p, 0)
    b = _find_primitive_root(p, g)

    got = smooth.DLP(p).solve(g, b)

    assert got != 0  # g, b both primitive roots => a genuine x is never 0
    assert pow(g, got, p) == b


def test_dlp_solve_repeated_calls_on_same_instance_are_stable():
    # Regression coverage for a real concurrency bug found during
    # development: DLP::solve runs its g-side/b-side collection+solve
    # pipelines concurrently, and two threads calling crtSolve at once
    # crashed the underlying LinBox/Givaro/NTL/FFLAS-FFPACK stack
    # nondeterministically (not always on the first call). Repeating solve()
    # on one shared instance is what would have caught that.
    p = 1009
    g = _find_primitive_root(p, 0)
    b = _find_primitive_root(p, g)

    solver = smooth.DLP(p)
    for _ in range(15):
        got = solver.solve(g, b)
        assert pow(g, got, p) == b


def test_dlp_solve_inverse_relationship_when_swapping_roles():
    # b = g**x implies g = b**y with x*y === 1 (mod p-1), true here because
    # b is required to be a primitive root too (a value is a primitive
    # root's image under exponentiation by x iff x is itself a unit mod
    # p-1). Checking both directions catches a wrong-but-internally-
    # consistent pair of answers that a one-directional check would miss.
    p = 8191
    g = _find_primitive_root(p, 0)
    b = _find_primitive_root(p, g)

    solver = smooth.DLP(p)
    x = solver.solve(g, b)
    y = solver.solve(b, g)

    assert pow(g, x, p) == b
    assert pow(b, y, p) == g
    assert (x * y) % (p - 1) == 1
