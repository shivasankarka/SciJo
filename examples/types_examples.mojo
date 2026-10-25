"""
A tour of `scijo.types`: a function written once, generically over
`Numeric`, works unchanged on a plain `Primal`, on a `Complex` value, and
on a `Dual` for an exact derivative.

```bash
mojo run -I . -I tests/ examples/types_examples.mojo
```
"""

# ===----------------------------------------------------------------------=== #
# SciJo
# ===----------------------------------------------------------------------=== #
from scijo.types import Complex, Numeric, Primal, variable, variable2
from scijo.prelude import *


def main() raises:
    plain_primal_call()
    exact_derivative_via_dual()
    second_derivative_via_nested_dual()
    same_kernel_on_complex()


# ===----------------------------------------------------------------------=== #
# The one kernel every example below reuses, written once against `Numeric`.
# ===----------------------------------------------------------------------=== #


def f[T: Numeric](x: T) -> T:
    """F(x) = x^3 - 2x^2 + x - 5, f'(x) = 3x^2 - 4x + 1, f''(x) = 6x - 4."""
    return x * x * x - 2.0 * x * x + x - 5.0


# ===----------------------------------------------------------------------=== #
# Level 1: a plain `Primal`, no derivative involved.
# ===----------------------------------------------------------------------=== #


def plain_primal_call() raises:
    print("=" * 80)
    print("LEVEL 1: f on a plain Primal")
    print("=" * 80)

    var result = f(Primal[f64](3.0))
    print("f(3.0)           =", result.v, " (expect 7.0)")


# ===----------------------------------------------------------------------=== #
# Level 2: the same `f`, called through `Dual` for a free exact derivative -
# no hand-written `f'` anywhere.
# ===----------------------------------------------------------------------=== #


def exact_derivative_via_dual() raises:
    print()
    print("=" * 80)
    print("LEVEL 2: the same f, called through Dual for a free derivative")
    print("=" * 80)

    var x0 = 3.0
    var result = f(variable(x0))
    print("f(3.0)           =", result.value.v, " (expect 7.0)")
    print("f'(3.0)          =", result.deriv.v, " (expect 3*9 - 4*3 + 1 = 16)")


# ===----------------------------------------------------------------------=== #
# Level 2, nested: `variable2` seeds a `Dual[Dual[Primal]]`, so the same
# chain rule applied twice gets the second derivative for free too.
# ===----------------------------------------------------------------------=== #


def second_derivative_via_nested_dual() raises:
    print()
    print("=" * 80)
    print("Nested Dual: f''(x) from the same f, no new code for the order")
    print("=" * 80)

    var x0 = 3.0
    var result = f(variable2(x0))
    print("f(3.0)           =", result.value.value.v, " (expect 7.0)")
    print("f'(3.0)          =", result.deriv.value.v, " (expect 16.0)")
    print("f''(3.0)         =", result.deriv.deriv.v, " (expect 6*3 - 4 = 14)")


# ===----------------------------------------------------------------------=== #
# Level 2, complex: the same `f`, called on a `Complex` value - still no
# rewrite, since `Complex` conforms to `Numeric` too.
# ===----------------------------------------------------------------------=== #


def same_kernel_on_complex() raises:
    print()
    print("=" * 80)
    print("Complex: the same f, called on a Complex value")
    print("=" * 80)

    var z = Complex[f64](1.0, 2.0)
    var result = f(z)
    print("f(1+2i)          =", result.re, "+", result.im, "i")
