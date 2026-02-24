# Changelog

All notable changes to this project will be documented in this file.

The format is based on "Keep a Changelog" and follows Semantic Versioning.

## [Unreleased] - v0.2
### Added
- `simpson` integration method in the `integrate` module. (commit: 34644a6)
- `bisect` (bisection) method for `root_scalar`. (commit: 62d1bd6)
- `secant` method for derivative-free root finding. (commit: fbe74ab)
- `RootResults` result struct for root-finding routines. (commit: d2fafb9)
- `root_scalar` entry-point with Newton's method. (commit: 51c3dfe)
- `optimize.__init__` now re-exports all public symbols (`root_scalar`, `newton`, `bisect`, `secant`).
- Developer guide (`docs/developer_guide.md`) covering file headers, docstrings, naming, testing, and contribution conventions.

### Changed
- `root_scalar` accepts additional optional arguments for method selection and tolerances. (commit: 52dc2a2)
- Bumped package version to v0.2. (commit: 2600da2)
- Standardized license header block across all files in `differentiate`, `fft`, `integrate`, `interpolate`, and `optimize` submodules.
- Added module-level docstrings to all `__init__.mojo` and implementation files.
- Rewrote all function and struct docstrings to follow the Mojo docstring style guide (consistent `Parameters:`, `Args:`, `Returns:`, `Raises:` sections).
- Added per-field docstrings to all result structs (`DiffResult`, `IntegralResult`, `RootResults`).
- Standardized `Args:` label everywhere (replaced `Arguments:` in `integrate.fixed_sample`).
- Standardized parameter descriptions (e.g., `"The floating-point data type."` across all modules).
- Renamed misleading `central_diff` variable to `diff_estimate` in forward/backward derivative functions.

### Fixed
- Fixed `Scalar[Self.dtype]` typo in `generate_backward_finite_difference_table` return docstring (standalone function, not a method).

### Removed

### Security
- No security-related changes in this release.

## [0.1] - initial release
- Initial public release (baseline for v0.2)
