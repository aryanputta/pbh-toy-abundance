# Toy PBH Abundance Calculator

This repo is a compact, reproducible toy model for primordial black hole (PBH)
abundance estimates from a lognormal primordial curvature power spectrum.

The goal is computational clarity: encode the physics assumptions, expose the
parameters, generate plots, and make the limitations obvious.

## Physics Model

Input spectrum:

```text
P_R(k) = A exp[-log(k / k0)^2 / (2 width^2)]
```

The code maps wavenumber to horizon mass using the common radiation-era scaling
`M ~ k^-2`, estimates a smoothed density variance with a Gaussian window, and
uses a Press-Schechter-like Gaussian tail:

```text
beta(M) = 0.5 erfc(delta_c / sqrt(2 sigma^2))
```

This is a toy approximation. It is not a substitute for full transfer
functions, non-Gaussian statistics, critical collapse, or numerical relativity.

## Run

```bash
python toy_abundance.py --make-plots
pytest
```

Generated files are written to `results/`.

## What To Look For

- Larger spectrum amplitude should increase PBH abundance.
- Wider spectra spread abundance across a broader mass range.
- Outputs are useful for code demonstration and method discussion, not for
  publication-level exclusion claims.

