import math
import struct
from sympy import Ei
from math import log, ceil, exp, pi
from bisect import bisect_left, bisect_right

_ANCHOR = [
    (2.0, 1e-3),
    (10.0, 1e-2),
    (20.0, 1e-1)
]

def build_delta_list():
    deltas = []
    for (u0, d0), (u1, d1) in zip(_ANCHOR[:-1], _ANCHOR[1:]):
        span = u1 - u0
        log_d0, log_d1 = math.log(d0), math.log(d1)

        u = u0
        while True:
            t = (u - u0) / span
            d = math.exp(log_d0 + t * (log_d1 - log_d0))

            if u + d >= u1 - 1e-12:
                deltas.append(u1 - u)
                break

            deltas.append(d)
            u += d
    return deltas

def build_tables(deltas):
    from sage.all import dickman_rho

    u_vals = [2.0]
    for d in deltas:
        u_vals.append(u_vals[-1] + d)

    log_rho_vals = [math.log(dickman_rho(u)) for u in u_vals]
    return u_vals, log_rho_vals