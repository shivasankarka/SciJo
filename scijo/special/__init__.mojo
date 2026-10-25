# ===----------------------------------------------------------------------=== #
# SciJo: Special functions module for Mojo
# Distributed under the Apache 2.0 License.
# ===----------------------------------------------------------------------=== #
"""Special Functions (`scijo.special`)
======================================
Gamma, error, beta, digamma, regularized incomplete gamma/beta, and Bessel
functions.

Every function here has two overloads: a concrete one taking and returning
a plain `Scalar[dtype]`, and a generic one over `scijo.types.NumericComparable`.
The concrete overload is a thin wrapper around the generic one via `Primal`;
the generic algorithm is the only implementation, so a power user who wants an exact
derivative just calls the same function on a `Dual` instead of a `Primal`.

Examples
--------
    ```mojo
    from scijo.special import gamma
    from scijo.prelude import *

    var y = gamma(Scalar[f64](2.5))  # 1.3293403881791370...
    ```
"""

# ===----------------------------------------------------------------------=== #
# SciJo
# ===----------------------------------------------------------------------=== #
from .bessel import (
    bessel_i0,
    bessel_i1,
    bessel_j0,
    bessel_j1,
    bessel_k0,
    bessel_k1,
    bessel_y0,
    bessel_y1,
)
from .beta import (
    beta,
    lnbeta,
)
from .beta_inc import beta_inc
from .erf import (
    erf,
    erfc,
)
from .gamma import (
    gamma,
    lngamma,
)
from .gamma_inc import (
    gamma_inc,
    gamma_incc,
)
from .psi import digamma
