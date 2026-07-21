"""Runs the full DLP.solve() benchmark suite end to end: Benchmark A
(average-case bit-length sweep) + Benchmark B (worst-case structural stress
tests), writes raw CSVs, prints summary statistics, and saves plots.

Usage (from the repo root, with the smoothnumbers-dev conda env active):

    pip install -r benchmarks/requirements.txt   # once, for matplotlib
    python -m benchmarks.run_all

See benchmarks/README.md for the full methodology writeup.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmarks import analyze, bench_average_case, bench_worst_case, common as bc


def main() -> None:
    t_start = time.time()

    print("=" * 70)
    print("Benchmark A: average case (bit-length sweep)")
    print("=" * 70)
    average_rows = bench_average_case.run()
    bc.write_csv(average_rows, bc.RESULTS_DIR / "average_case.csv")

    print("\n" + "=" * 70)
    print("Benchmark B: worst case (p-1 shape + Wiedemann-threshold straddle)")
    print("=" * 70)
    worst_rows = bench_worst_case.run()
    bc.write_csv(worst_rows, bc.RESULTS_DIR / "worst_case.csv")

    bc.write_csv(average_rows + worst_rows, bc.RESULTS_DIR / "all_results.csv")

    print("\n" + "=" * 70)
    print("Analysis")
    print("=" * 70)
    analyze.run_analysis(average_rows, worst_rows, bc.RESULTS_DIR)

    print(f"\nTotal suite wall time: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    main()
