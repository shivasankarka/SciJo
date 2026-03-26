# Changelog

All notable changes to this project will be documented in this file.

## v0.2.0
### Added
- `simpson` integration method in the `integrate` module. (commit: 34644a6)
- `bisect` (bisection) method for `root_scalar`. (commit: 62d1bd6)
- `secant` method for derivative-free root finding. (commit: fbe74ab)
- `RootResults` result struct for root-finding routines. (commit: d2fafb9)
- `root_scalar` entry-point with Newton's method. (commit: 51c3dfe)
- `optimize.__init__` now re-exports all public symbols (`root_scalar`, `newton`, `bisect`, `secant`).
- `jacobian` for vector-valued functions. (commit: b03b1e7)
- `romb` (Romberg) integration method in `integrate.fixed_sample`. (commit: 3d33047)
- `quad` QNG Gauss-Kronrod integrator and tests. (commits: a2275df, 0e35c17)
- Cooley–Tukey FFT/IFFT implementation and tests. (commits: bf025e0, f7a58ba)
- Test runner script and CI workflow. (commits: 6063209, ec07599, 14d29a0)
- Benchmark scripts for derivatives/performance comparisons. (commits: 8718e0d, 5c32a4e, 9bdb3bc)
- Top-level docs in `docs/`. (commit: e062ee7)
- Developer guide (`docs/developer_guide.md`) covering file headers, docstrings, naming, testing, and contribution conventions.

### Changed
- `root_scalar` accepts additional optional arguments for method selection and tolerances. (commit: 52dc2a2)
- Bumped package version to v0.2. (commit: 2600da2)
- `trapezoid` moved under `integrate.fixed_sample` (file rename). (commits: e9f49b7, 368b685)
- Cleaned up and reorganized constants, differentiate, integrate, interpolate, fft, and optimize modules. (commits: e7b8819, ef79f2b, b7f8a52, 179087d, 68dc6be, 20b48af, 870b23d, dfe0c4c, 6c50822, 80d33f5, e203451)
- Updated Mojo version and pixi dependency configs. (commits: 2910b9f, 2603c3c, 6c978e1, bf2eee5, dc0f7e4)
- Standardized license header block across all files in `differentiate`, `fft`, `integrate`, `interpolate`, and `optimize` submodules.
- Added module-level docstrings to all `__init__.mojo` and implementation files.
- Rewrote all function and struct docstrings to follow the Mojo docstring style guide (consistent `Parameters:`, `Args:`, `Returns:`, `Raises:` sections).
- Added per-field docstrings to all result structs (`DiffResult`, `IntegralResult`, `RootResults`).
- Standardized `Args:` label everywhere (replaced `Arguments:` in `integrate.fixed_sample`).
- Standardized parameter descriptions (e.g., `"The floating-point data type."` across all modules).
- Renamed misleading `central_diff` variable to `diff_estimate` in forward/backward derivative functions.

### Fixed
- Fixed `Scalar[Self.dtype]` typo in `generate_backward_finite_difference_table` return docstring (standalone function, not a method).
- Fixed `trapezoid` overload to avoid array copies. (commit: b5ec86c)
- Formatting/typo fixes in `quad` and docstrings. (commits: 5972a82, 579083f)

### Removed

### Security
- No security-related changes in this release.

## [0.1] - initial release
- Initial public release (baseline for v0.2)
