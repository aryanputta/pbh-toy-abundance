"""Toy PBH abundance estimates from a lognormal curvature spectrum.

The formulas here are deliberately lightweight. They are useful for checking
parameter dependence and code structure, not for precision PBH constraints.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

M_SUN_G = 1.98847e33
K_REF_MPC = 2.9e5
M_REF_SOLAR = 30.0


def curvature_power(k: np.ndarray, amplitude: float, k0: float, width: float) -> np.ndarray:
    """Lognormal primordial curvature power spectrum."""
    if amplitude <= 0:
        raise ValueError("amplitude must be positive")
    if k0 <= 0:
        raise ValueError("k0 must be positive")
    if width <= 0:
        raise ValueError("width must be positive")
    return amplitude * np.exp(-0.5 * (np.log(k / k0) / width) ** 2)


def horizon_mass_solar(k: np.ndarray) -> np.ndarray:
    """Approximate horizon mass in solar masses for radiation-era re-entry."""
    return M_REF_SOLAR * (k / K_REF_MPC) ** -2


def density_variance(
    k_eval: np.ndarray,
    amplitude: float,
    k0: float,
    width: float,
    smoothing_width: float = 0.7,
) -> np.ndarray:
    """Estimate smoothed density variance around each k using a log-Gaussian window."""
    log_grid = np.linspace(np.log(k_eval.min()) - 5 * width, np.log(k_eval.max()) + 5 * width, 900)
    k_grid = np.exp(log_grid)
    power = curvature_power(k_grid, amplitude, k0, width)
    out = []
    for k in k_eval:
        window = np.exp(-0.5 * (np.log(k_grid / k) / smoothing_width) ** 2)
        norm = np.trapezoid(window, x=log_grid)
        smoothed_power = np.trapezoid(power * window, x=log_grid) / norm
        out.append((4.0 / 9.0) ** 2 * smoothed_power)
    return np.asarray(out)


def collapse_fraction(variance: np.ndarray, delta_c: float = 0.45) -> np.ndarray:
    """Gaussian-tail collapse fraction beta."""
    if delta_c <= 0:
        raise ValueError("delta_c must be positive")
    variance = np.asarray(variance)
    if np.any(variance <= 0):
        raise ValueError("variance values must be positive")
    erfc_vec = np.vectorize(math.erfc)
    return 0.5 * erfc_vec(delta_c / np.sqrt(2.0 * variance))


def present_fraction(beta: np.ndarray, mass_solar: np.ndarray, gamma: float = 0.2) -> np.ndarray:
    """Approximate present-day PBH dark matter fraction f_PBH."""
    return 2.7e8 * beta * math.sqrt(gamma / 0.2) * mass_solar ** -0.5


def abundance_curve(
    amplitude: float = 0.02,
    width: float = 0.5,
    k0: float = K_REF_MPC,
    n_points: int = 160,
) -> dict[str, np.ndarray]:
    """Compute mass, beta, and f_PBH arrays for one spectrum."""
    k = np.geomspace(k0 / 50.0, k0 * 50.0, n_points)
    mass = horizon_mass_solar(k)
    variance = density_variance(k, amplitude, k0, width)
    beta = collapse_fraction(variance)
    f_pbh = present_fraction(beta, mass)
    order = np.argsort(mass)
    return {
        "k_mpc": k[order],
        "mass_solar": mass[order],
        "variance": variance[order],
        "beta": beta[order],
        "f_pbh": f_pbh[order],
    }


def write_csv(curve: dict[str, np.ndarray], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(curve.keys())
        for row in zip(*curve.values()):
            writer.writerow([f"{value:.8e}" for value in row])


def make_plots(results_dir: Path) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for amp in [0.006, 0.01, 0.015, 0.02]:
        curve = abundance_curve(amplitude=amp, width=0.5)
        ax.loglog(curve["mass_solar"], np.maximum(curve["f_pbh"], 1e-300), label=f"A={amp:g}")
    ax.set_xlabel("PBH mass [$M_\\odot$]")
    ax.set_ylabel("Toy present fraction $f_{PBH}$")
    ax.set_title("Abundance sensitivity to spectrum amplitude")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(results_dir / "abundance_vs_amplitude.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for width in [0.25, 0.4, 0.7, 1.0]:
        curve = abundance_curve(amplitude=0.015, width=width)
        ax.loglog(curve["mass_solar"], np.maximum(curve["f_pbh"], 1e-300), label=f"width={width:g}")
    ax.set_xlabel("PBH mass [$M_\\odot$]")
    ax.set_ylabel("Toy present fraction $f_{PBH}$")
    ax.set_title("Abundance sensitivity to spectrum width")
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(results_dir / "abundance_vs_width.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--amplitude", type=float, default=0.015)
    parser.add_argument("--width", type=float, default=0.5)
    parser.add_argument("--k0", type=float, default=K_REF_MPC)
    parser.add_argument("--make-plots", action="store_true")
    args = parser.parse_args()

    results_dir = Path("results")
    curve = abundance_curve(args.amplitude, args.width, args.k0)
    write_csv(curve, results_dir / "toy_abundance_curve.csv")
    if args.make_plots:
        make_plots(results_dir)

    peak = int(np.argmax(curve["f_pbh"]))
    print(f"peak_mass_solar={curve['mass_solar'][peak]:.6e}")
    print(f"peak_f_pbh={curve['f_pbh'][peak]:.6e}")


if __name__ == "__main__":
    main()
