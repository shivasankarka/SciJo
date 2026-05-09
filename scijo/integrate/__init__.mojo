# ===----------------------------------------------------------------------=== #
# SciJo: Integrate module for Mojo
# Distributed under the Apache 2.0 License.
# ===----------------------------------------------------------------------=== #
"""Integrate Module (`scijo.integrate`)
=======================================

Provides tools for numerical integration and quadrature. It includes adaptive
and non-adaptive methods for computing definite integrals, as well as
fixed-sample integration rules for discrete data.

Available Functions
-------------------
- `quad`        — General-purpose adaptive quadrature (Gauss-Kronrod).
- `trapezoid`   — Composite trapezoidal rule for discrete data.
- `simpson`     — Simpson's rule for discrete data.
- `romb`        — Romberg integration with Richardson extrapolation.

Examples
--------
    ```mojo
    from scijo.integrate import quad, trapezoid

    def integrand[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
        return x * x

    var result = quad[f64, integrand](0.0, 1.0, None)
    ```
"""

from .quadrature import quad
from .fixed_sample import trapezoid, simpson, romb
