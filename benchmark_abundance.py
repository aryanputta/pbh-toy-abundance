"""Benchmark the toy PBH abundance code against a direct reference formula.

This is an internal numerical benchmark, not a comparison to GRChombo or the
Einstein Toolkit. External numerical-relativity outputs can be added as CSV
files later and compared with the same plotting style.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from toy_abundance import abundance_curve


def direct_gaussian_tail(amplitude: float, mass_solar: np.ndarray, delta_c: float = 0.45) -> np.ndarray:
    """Reference estimate using the unsmoothed peak variance at every mass."""
    variance = np.full_like(mass_solar, (4.0 / 9.0) ** 2 * amplitude)
    erfc_vec = np.vectorize(math.erfc)
    beta = 0.5 * erfc_vec(delta_c / np.sqrt(2.0 * variance))
    return 2.7e8 * beta * mass_solar**-0.5


def run_benchmark(amplitude: float = 0.015, width: float = 0.5) -> dict[str, np.ndarray]:
    curve = abundance_curve(amplitude=amplitude, width=width, n_points=120)
    reference = direct_gaussian_tail(amplitude, curve["mass_solar"])
    ratio = np.divide(curve["f_pbh"], reference, out=np.zeros_like(curve["f_pbh"]), where=reference > 0)
    return {
        "mass_solar": curve["mass_solar"],
        "toy_f_pbh": curve["f_pbh"],
        "reference_f_pbh": reference,
        "toy_to_reference_ratio": ratio,
    }


def write_csv(data: dict[str, np.ndarray], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(data.keys())
        for row in zip(*data.values()):
            writer.writerow([f"{value:.8e}" for value in row])


def make_plot(data: dict[str, np.ndarray], path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 6.2), sharex=True)
    axes[0].loglog(data["mass_solar"], np.maximum(data["toy_f_pbh"], 1e-300), label="smoothed toy model")
    axes[0].loglog(
        data["mass_solar"],
        np.maximum(data["reference_f_pbh"], 1e-300),
        "--",
        label="direct Gaussian-tail reference",
    )
    axes[0].set_ylabel("$f_{PBH}$")
    axes[0].set_title("Toy abundance benchmark comparison")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend()
    axes[1].semilogx(data["mass_solar"], data["toy_to_reference_ratio"])
    axes[1].axhline(1.0, color="black", linewidth=1.0, linestyle="--")
    axes[1].set_xlabel("PBH mass [$M_\\odot$]")
    axes[1].set_ylabel("toy / reference")
    axes[1].grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    results_dir = Path("results")
    data = run_benchmark()
    write_csv(data, results_dir / "abundance_benchmark.csv")
    make_plot(data, results_dir / "abundance_benchmark_compare.png")
    print(f"benchmark_points={len(data['mass_solar'])}")
    print(f"finite_ratios={np.isfinite(data['toy_to_reference_ratio']).all()}")


if __name__ == "__main__":
    main()

