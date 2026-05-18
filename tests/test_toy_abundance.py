import numpy as np

from benchmark_abundance import run_benchmark
from toy_abundance import abundance_curve, collapse_fraction


def test_collapse_fraction_increases_with_variance():
    beta = collapse_fraction(np.array([1e-4, 1e-3, 1e-2]))
    assert np.all(np.diff(beta) > 0)


def test_abundance_curve_is_finite_and_positive():
    curve = abundance_curve(amplitude=0.015, width=0.5, n_points=32)
    assert set(curve) == {"k_mpc", "mass_solar", "variance", "beta", "f_pbh"}
    assert np.all(np.isfinite(curve["f_pbh"]))
    assert np.all(curve["mass_solar"] > 0)


def test_larger_amplitude_gives_larger_peak_abundance():
    low = abundance_curve(amplitude=0.008, width=0.5, n_points=40)
    high = abundance_curve(amplitude=0.015, width=0.5, n_points=40)
    assert high["f_pbh"].max() > low["f_pbh"].max()


def test_benchmark_ratios_are_finite():
    benchmark = run_benchmark()
    assert np.all(np.isfinite(benchmark["toy_to_reference_ratio"]))
    assert benchmark["toy_f_pbh"].max() > 0
