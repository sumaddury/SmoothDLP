import matplotlib.pyplot as plt
import time
from statistics import mean
import numpy as np
from math import isqrt

def sieveTo(n):
    primality = [False,False]+[True for c in range(n-1)]
    for i in range(2, isqrt(n)+1):
        if not primality[i]:
            continue
        for j in range(i*i, n + 1, i):
            primality[j] = False
    return [p for p in range(len(primality)) if primality[p]]

def primality_complexity(func, N, data=None):
    if data:
        test_primes = data
    else:
        print(f"Beginning sieveTo({N})...")
        test_primes = sieveTo(N)
        print(f"Found {len(test_primes)} primes.")

    times = []
    print("Collecting times...")
    for i, p in enumerate(test_primes):
        print(str(i)+', ', end='')
        t1 = time.time()
        func(n=p)
        t2 = time.time()
        times.append(t2 - t1)

    print("\nFinished collecting times.")
    return times, test_primes

def graph_complexity(times, test_primes, dpts=1000, log=False):
    delta = len(times) // dpts
    T = [mean(times[i : i + delta]) for i in range(0, len(times), delta)]
    P = [mean(test_primes[i : i + delta]) for i in range(0, len(test_primes), delta)]
    print("Finished concat.")
    print("Graphing...")

    fig = plt.figure(figsize=(10, 6), dpi=200)
    plt.scatter(P, T, marker='.', s=1, color='k')
    if log:
        plt.xscale('log', base=2)
        plt.yscale('log', base=2)
    xmin, xmax = min(P), max(P)
    ymin, ymax = min(T), max(T)
    print(xmin, xmax, ymin, ymax)
    
    if log:
        log2_min = np.log2(float(xmin))
        log2_max = np.log2(float(xmax))
        ticks_log2 = np.linspace(log2_min, log2_max, 10)
        xticks = 2 ** ticks_log2
        xtick_labels = [f"$2^{{{int(round(np.log2(float(x))))}}}$" for x in xticks]
        plt.xticks(xticks, xtick_labels)
    else:
        xticks = np.linspace(xmin, xmax, 10)
        def round_highest(x):
            if x == 0:
                return 0
            mag = 10 ** int(np.floor(np.log10(abs(float(x)))))
            return int(round(x / mag) * mag)
        xtick_labels = [str(round_highest(x)) for x in xticks]
        plt.xticks(xticks, xtick_labels)
    if log:
       log2_ymin = np.log2(float(ymin))
       log2_ymax = np.log2(float(ymax))
       yticks = 2 ** np.linspace(log2_ymin, log2_ymax, 10)
       ytick_labels = [f"$2^{{{int(round(np.log2(float(y))))}}}$" for y in yticks]
       plt.yticks(yticks, ytick_labels)
    else:
       yticks = np.linspace(ymin, ymax, 10)
       ytick_labels = [format(y, '.3e') for y in yticks]
       plt.yticks(yticks, ytick_labels)
    plt.tick_params(axis='both', direction='in')
    plt.xlabel('n')
    plt.ylabel('Average Time (sec)')
    plt.tight_layout()
    plt.show()