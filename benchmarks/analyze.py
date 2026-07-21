"""Summary statistics and plots for the DLP.solve() benchmark suite results.

Requires matplotlib (benchmarks/requirements.txt) -- not a dependency of the
`smooth` package itself, only of this analysis step.
"""

from __future__ import annotations

import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmarks import common as bc

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Fixed-order categorical palette (Okabe-Ito, colorblind-safe) -- one color
# per distinct series identity, never reassigned/cycled.
COLOR_TREND = "#0072B2"       # blue    -- average-case mean/trend
COLOR_REFERENCE = "#999999"   # gray    -- external reference point (paper)
COLOR_THRESHOLD = "#D55E00"   # vermillion -- Wiedemann-threshold straddle points
COLOR_CONTROL = "#009E73"     # green   -- generic 60-bit control group
COLOR_SAFE = "#0072B2"        # blue    -- safe_prime group
COLOR_SMOOTH = "#E69F00"      # orange  -- smooth_pm1 group

PAPER_REFERENCE_BITS = 70
PAPER_REFERENCE_SECONDS = 5000.0


# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------

def summarize_by(rows: list[dict], group_key: str) -> dict:
    groups: dict = {}
    for r in rows:
        if r["elapsed_s"] is None:
            continue
        groups.setdefault(r[group_key], []).append(r["elapsed_s"])

    summary = {}
    for key, times in sorted(groups.items(), key=lambda kv: str(kv[0])):
        summary[key] = {
            "n": len(times),
            "mean": statistics.mean(times),
            "median": statistics.median(times),
            "stdev": statistics.stdev(times) if len(times) > 1 else 0.0,
            "min": min(times),
            "max": max(times),
        }
    return summary


def print_summary_table(title: str, summary: dict) -> None:
    print(f"\n{title}")
    print(f"{'key':>16} {'n':>4} {'mean':>10} {'median':>10} {'stdev':>10} {'min':>10} {'max':>10}")
    for key, s in summary.items():
        print(
            f"{str(key):>16} {s['n']:>4} {s['mean']:>10.3f} {s['median']:>10.3f} "
            f"{s['stdev']:>10.3f} {s['min']:>10.3f} {s['max']:>10.3f}"
        )


def variance_decomposition(average_rows: list[dict]) -> dict:
    """Splits Benchmark A's timing spread, per bit length, into a
    *between-prime* component (stdev of each prime's mean solve time) and a
    *within-prime* component (mean stdev across the R pairs solved on the
    same prime). Large between-prime relative to within-prime is evidence
    that a prime's own structure (not sampling luck) drives its solve time;
    the reverse means solve time is mostly noise around a bit-length trend.
    """
    by_bits_prime = defaultdict(list)
    for r in average_rows:
        if r["elapsed_s"] is None:
            continue
        by_bits_prime[(r["bits"], r["prime_idx"])].append(r["elapsed_s"])

    by_bits = defaultdict(list)
    for (bits, _prime_idx), times in by_bits_prime.items():
        by_bits[bits].append(times)

    result = {}
    for bits, prime_groups in sorted(by_bits.items()):
        prime_means = [statistics.mean(g) for g in prime_groups]
        between = statistics.stdev(prime_means) if len(prime_means) > 1 else 0.0
        within_stdevs = [statistics.stdev(g) for g in prime_groups if len(g) > 1]
        within = statistics.mean(within_stdevs) if within_stdevs else 0.0
        result[bits] = {
            "between_prime_stdev": between,
            "within_prime_stdev": within,
            "num_primes": len(prime_groups),
        }
    return result


def print_variance_decomposition(decomp: dict) -> None:
    print("\nBenchmark A variance decomposition (seconds)")
    print(f"{'bits':>6} {'#primes':>8} {'between-prime stdev':>22} {'within-prime stdev':>22}")
    for bits, d in decomp.items():
        print(f"{bits:>6} {d['num_primes']:>8} {d['between_prime_stdev']:>22.3f} {d['within_prime_stdev']:>22.3f}")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_average_case_scaling(average_rows: list[dict], threshold_rows: list[dict], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=150)

    by_bits = defaultdict(list)
    for r in average_rows:
        if r["elapsed_s"] is not None:
            by_bits[r["bits"]].append(r["elapsed_s"])

    jitter_rng = random.Random("plot-jitter")
    for bits, times in by_bits.items():
        xs = [bits + jitter_rng.uniform(-0.6, 0.6) for _ in times]
        ax.scatter(xs, times, color=COLOR_TREND, alpha=0.35, s=22, linewidths=0, zorder=2)

    mean_bits = sorted(by_bits.keys())
    mean_times = [statistics.mean(by_bits[b]) for b in mean_bits]
    ax.plot(mean_bits, mean_times, color=COLOR_TREND, marker="o", linewidth=2, zorder=3, label="mean solve() time")

    ax.scatter(
        [PAPER_REFERENCE_BITS], [PAPER_REFERENCE_SECONDS],
        color=COLOR_REFERENCE, marker="D", s=70, zorder=4,
        label=f"paper (2024) reported ~{PAPER_REFERENCE_SECONDS:.0f}s @ {PAPER_REFERENCE_BITS}-bit",
    )

    threshold_by_bits = defaultdict(list)
    for r in threshold_rows:
        if r["elapsed_s"] is not None:
            threshold_by_bits[r["bits"]].append(r["elapsed_s"])
    if threshold_by_bits:
        tx = list(threshold_by_bits.keys())
        ty = [statistics.mean(v) for v in threshold_by_bits.values()]
        ax.scatter(
            tx, ty, color=COLOR_THRESHOLD, marker="*", s=220, zorder=5,
            label="Wiedemann-threshold straddle (single trial, 76b/84b)",
        )

    ax.set_yscale("log")
    ax.set_xlabel("prime bit length")
    ax.set_ylabel("DLP.solve() wall time (s, log scale)")
    ax.set_title("DLP.solve() wall time vs. prime bit length")
    ax.grid(True, which="both", alpha=0.25, color="#cccccc")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.legend(frameon=False, loc="upper left")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_pm1_shape_comparison(worst_rows: list[dict], average_rows: list[dict], out_path: Path) -> None:
    control_times = [r["elapsed_s"] for r in average_rows if r["bits"] == 60 and r["elapsed_s"] is not None]
    safe_times = [r["elapsed_s"] for r in worst_rows if r["label"] == "safe_prime" and r["elapsed_s"] is not None]
    smooth_times = [r["elapsed_s"] for r in worst_rows if r["label"] == "smooth_pm1" and r["elapsed_s"] is not None]

    groups = [
        ("generic\n(Benchmark A control)", control_times, COLOR_CONTROL),
        ("safe prime\n(p-1 = 2q)", safe_times, COLOR_SAFE),
        ("smooth p-1\n(many small factors)", smooth_times, COLOR_SMOOTH),
    ]

    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=150)
    jitter_rng = random.Random("plot-jitter-pm1")

    for i, (label, times, color) in enumerate(groups):
        if not times:
            continue
        xs = [i + jitter_rng.uniform(-0.12, 0.12) for _ in times]
        ax.scatter(xs, times, color=color, alpha=0.6, s=40, zorder=2)
        mean = statistics.mean(times)
        ax.hlines(mean, i - 0.25, i + 0.25, color=color, linewidth=2.5, zorder=3)

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylabel("DLP.solve() wall time (s)")
    ax.set_title("Solve time vs. p-1 factorization shape (60-bit primes)")
    ax.grid(True, axis="y", alpha=0.25, color="#cccccc")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def run_analysis(average_rows: list[dict], worst_rows: list[dict], results_dir: Path) -> None:
    threshold_rows = [r for r in worst_rows if r["label"].startswith("threshold_")]

    print_summary_table("Benchmark A -- average case, by bit length", summarize_by(average_rows, "bits"))
    print_variance_decomposition(variance_decomposition(average_rows))
    print_summary_table("Benchmark B -- p-1 shape / threshold, by label", summarize_by(worst_rows, "label"))

    timed_out = [r for r in (average_rows + worst_rows) if r["timed_out"]]
    if timed_out:
        print(f"\n{len(timed_out)} trial(s) hit their timeout cap (right-censored, excluded from stats above):")
        for r in timed_out:
            print(f"  bench={r['bench']} bits={r['bits']} p={r['p']}")

    results_dir.mkdir(parents=True, exist_ok=True)
    plot_average_case_scaling(average_rows, threshold_rows, results_dir / "average_case_scaling.png")
    plot_pm1_shape_comparison(worst_rows, average_rows, results_dir / "pm1_shape_comparison.png")
    print(f"\nwrote plots to {results_dir}")
