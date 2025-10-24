# Changelog

All notable changes to this project will be documented in this file.

The format is based on "Keep a Changelog" and follows Semantic Versioning.

## [Unreleased] - v0.2
### Added
- Add `simpson` integration method to the `integration` module. (commit: 34644a6)
- Add `bisect` (bisection) method to `root_scalar`, with improved docstrings and usage examples. (commit: 62d1bd6)
- Add `secant` method to `root_scalar` for derivative-free root finding. (commit: fbe74ab)
- Introduce `RootResults` (a result object for root-finding routines) to standardize return values and make results easier to inspect. (commit: d2fafb9)
- Add a basic `root_scalar` optimizer implementing the Newton method to provide a simple optimization entry-point. (commit: 51c3dfe)

### Changed
- `root_scalar` now accepts additional optional arguments to configure methods and tolerances; this makes the API more flexible while keeping sensible defaults. (commit: 52dc2a2)
- Bumped package version reference to v0.2. (commit: 2600da2)
- Improved docstrings and examples for root-finding routines (notably `bisect`) to clarify behavior and return types.

### Fixed
- Documentation and docstring fixes associated with the new root-finding methods (included alongside method additions).

### Deprecated
- No deprecations in this release.

### Removed
- No removals in this release.

### Security
- No security-related changes in this release.

## [0.1] - initial release
- Initial public release (baseline for v0.2)
