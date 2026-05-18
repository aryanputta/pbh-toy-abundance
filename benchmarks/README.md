# External Benchmark Hook

This folder is for future numerical-relativity comparison data.

Expected CSV format for an external PBH collapse or abundance reference:

```csv
mass_solar,f_pbh,source
1.0e-12,1.0e-3,grchombo-example
```

The current benchmark in `benchmark_abundance.py` compares against a direct
analytic Gaussian-tail reference. It does not claim validation against
GRChombo, Einstein Toolkit, or any full numerical-relativity code.

