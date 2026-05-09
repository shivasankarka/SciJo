# ===----------------------------------------------------------------------=== #
# SciJo: Interpolate module for Mojo
# Distributed under the Apache 2.0 License.
# ===----------------------------------------------------------------------=== #
"""Interpolate Module (`scijo.interpolate`)
===========================================
Provides tools for interpolating data. It includes linear interpolation
functions and callable interpolator objects for both single-point and
array-based evaluation.

Available Functions
-------------------
- `interp1d`            — Create a callable linear interpolator or interpolate directly.
- `LinearInterpolator`  — A reusable callable interpolation object.

Examples
--------
    ```mojo
    from scijo.interpolate import interp1d

    var x = nm.arange[f64](0.0, 1.0, 0.5)
    var y = x * x
    var interp = interp1d(x, y, bounds_error=False, fill_value=0.0)
    var yq = interp(Scalar[f64](0.25))
    ```
"""

from .interpolate import interp1d, LinearInterpolator
