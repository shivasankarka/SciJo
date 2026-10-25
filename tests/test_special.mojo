# ===----------------------------------------------------------------------=== #
# Stdlib
# ===----------------------------------------------------------------------=== #
from std.math import exp, pi, sqrt
from std.testing import assert_almost_equal, TestSuite

# ===----------------------------------------------------------------------=== #
# SciJo
# ===----------------------------------------------------------------------=== #
from scijo import f64
from scijo.special import (
    beta,
    bessel_i0,
    bessel_i1,
    bessel_j0,
    bessel_j1,
    bessel_k0,
    bessel_k1,
    bessel_y0,
    bessel_y1,
    digamma,
    erf,
    erfc,
    gamma,
    lnbeta,
    lngamma,
)
from scijo.types import NumericComparable, variable


# ===----------------------------------------------------------------------=== #
# Tier 1 (`scijo.special`, plain `Scalar`): values against the scipy oracle
# (`pixi run -e test python3 -c "from scipy.special import ..."`, recorded
# in each assertion's comment) and against identities that need no oracle.
# ===----------------------------------------------------------------------=== #


def test_gamma_matches_oracle_and_factorial() raises:
    assert_almost_equal(gamma(Scalar[f64](5.0)), 24.0, msg="gamma(5) == 4!")
    assert_almost_equal(
        gamma(Scalar[f64](2.5)), 1.329340388179137, msg="gamma(2.5)"
    )
    assert_almost_equal(gamma(Scalar[f64](0.5)), sqrt(pi), msg="gamma(0.5)")


def test_lngamma_matches_oracle_and_stays_finite() raises:
    assert_almost_equal(
        lngamma(Scalar[f64](171.0)),
        706.5730622457875,
        msg="lngamma(171), where gamma(171) already overflows f64",
    )


def test_digamma_matches_oracle() raises:
    assert_almost_equal(
        digamma(Scalar[f64](1.0)), -0.5772156649015329, msg="digamma(1)"
    )
    assert_almost_equal(
        digamma(Scalar[f64](2.5)), 0.7031566406452432, msg="digamma(2.5)"
    )


def test_erf_matches_oracle_and_odd_symmetry() raises:
    assert_almost_equal(
        erf(Scalar[f64](1.0)), 0.8427007929497148, msg="erf(1)"
    )
    assert_almost_equal(
        erf(Scalar[f64](-1.0)), -0.8427007929497148, msg="erf is odd"
    )
    assert_almost_equal(
        erf(Scalar[f64](1.0)) + erfc(Scalar[f64](1.0)),
        1.0,
        msg="erf(x) + erfc(x) == 1",
    )
    assert_almost_equal(
        erfc(Scalar[f64](5.0)),
        1.5374597944280353e-12,
        atol=1e-20,
        msg="erfc(5), where 1 - erf(5) would lose all its digits",
    )


def test_beta_matches_oracle_and_gamma_relation() raises:
    assert_almost_equal(
        beta(Scalar[f64](2.0), Scalar[f64](3.0)),
        0.08333333333333333,
        msg="beta(2, 3)",
    )
    assert_almost_equal(
        lnbeta(Scalar[f64](2.5), Scalar[f64](4.5)),
        -3.8408317706947397,
        msg="lnbeta(2.5, 4.5)",
    )
    assert_almost_equal(
        beta(Scalar[f64](2.5), Scalar[f64](4.5)),
        beta(Scalar[f64](4.5), Scalar[f64](2.5)),
        msg="beta(a, b) == beta(b, a)",
    )


def test_bessel_matches_oracle() raises:
    assert_almost_equal(bessel_j0(Scalar[f64](2.0)), 0.22389077914123562)
    assert_almost_equal(bessel_j1(Scalar[f64](2.0)), 0.5767248077568734)
    assert_almost_equal(bessel_y0(Scalar[f64](2.0)), 0.5103756726497451)
    assert_almost_equal(bessel_y1(Scalar[f64](2.0)), -0.10703243154093758)
    assert_almost_equal(bessel_i0(Scalar[f64](2.0)), 2.279585302336067)
    assert_almost_equal(bessel_i1(Scalar[f64](2.0)), 1.5906368546373295)
    assert_almost_equal(bessel_k0(Scalar[f64](2.0)), 0.1138938727495334)
    assert_almost_equal(bessel_k1(Scalar[f64](2.0)), 0.13986588181652246)


def test_bessel_wronskian_identity() raises:
    """`J0(x)*Y1(x) - J1(x)*Y0(x) == -2/(pi*x)`, true for every `x`, and
    independent of the series/asymptotic split each side was computed with
    (checked at a series-branch and an asymptotic-branch point)."""
    for x in [2.0, 30.0]:
        var lhs = bessel_j0(Scalar[f64](x)) * bessel_y1(Scalar[f64](x)) - (
            bessel_j1(Scalar[f64](x)) * bessel_y0(Scalar[f64](x))
        )
        var rhs = -2.0 / (pi * x)
        assert_almost_equal(
            lhs, rhs, msg="Bessel Wronskian at x=" + String(x)
        )


def test_bessel_asymptotic_branch_matches_oracle() raises:
    """Past each function's series/asymptotic threshold (`x >= 15` for
    `J`/`Y`/`I`, `x >= 8` for `K`), it's the asymptotic expansion being
    checked here, not the series."""
    assert_almost_equal(
        bessel_j0(Scalar[f64](20.0)), 0.16702466434058313, msg="J0(20)"
    )
    assert_almost_equal(
        bessel_k0(Scalar[f64](10.0)),
        1.778006231616765e-05,
        atol=1e-12,
        msg="K0(10)",
    )


# ===----------------------------------------------------------------------=== #
# Tier 2 (the generic `scijo.special` overloads, run through `Dual`): every
# derivative here falls out of `Dual`'s own chain-rule arithmetic - there is
# no hand-written derivative formula anywhere in `scijo.special` to call.
# The expected values below are still the closed-form derivative, used only
# to check the *result* is exact, the same "no finite-difference tolerance
# needed" pattern as `tests/test_types.mojo`.
# ===----------------------------------------------------------------------=== #


def erf_fn[T: NumericComparable](x: T) -> T:
    return erf(x)


def gamma_fn[T: NumericComparable](x: T) -> T:
    return gamma(x)


def digamma_fn[T: NumericComparable](x: T) -> T:
    return digamma(x)


def bessel_j0_fn[T: NumericComparable](x: T) -> T:
    return bessel_j0(x)


def bessel_j1_fn[T: NumericComparable](x: T) -> T:
    return bessel_j1(x)


def bessel_y0_fn[T: NumericComparable](x: T) -> T:
    return bessel_y0(x)


def bessel_i1_fn[T: NumericComparable](x: T) -> T:
    return bessel_i1(x)


def bessel_k1_fn[T: NumericComparable](x: T) -> T:
    return bessel_k1(x)


def test_dual_erf_derivative() raises:
    var x0 = 0.8
    var x = variable(x0)
    var result = erf_fn(x)
    assert_almost_equal(result.value.v, erf(Scalar[f64](x0)), msg="erf(0.8)")
    var expected_deriv = (2.0 / sqrt(pi)) * exp(-(x0 * x0))
    assert_almost_equal(
        result.deriv.v,
        expected_deriv,
        msg="d/dx erf(x) = (2/sqrt(pi)) * exp(-x^2)",
    )


def test_dual_gamma_derivative() raises:
    """`d/dx Gamma(x) = Gamma(x) * psi(x)`, previously impossible to give
    `Dual` without an unbounded `digamma -> trigamma -> ...` dependency
    chain - no longer an issue, since there is no hand-written derivative
    formula to need one."""
    var x0 = 2.5
    var x = variable(x0)
    var result = gamma_fn(x)
    assert_almost_equal(result.value.v, gamma(Scalar[f64](x0)), msg="gamma(2.5)")
    var expected_deriv = gamma(Scalar[f64](x0)) * digamma(Scalar[f64](x0))
    assert_almost_equal(
        result.deriv.v, expected_deriv, msg="d/dx Gamma(x) = Gamma(x)*psi(x)"
    )


def test_dual_digamma_derivative() raises:
    """`d/dx psi(x) = trigamma(x)`, approximated here via a central
    finite difference since `scijo.special` has no `trigamma` oracle to
    compare against directly."""
    var x0 = 2.5
    var h = 1e-5
    var x = variable(x0)
    var result = digamma_fn(x)
    assert_almost_equal(result.value.v, digamma(Scalar[f64](x0)), msg="digamma(2.5)")
    var expected_deriv = (
        digamma(Scalar[f64](x0 + h)) - digamma(Scalar[f64](x0 - h))
    ) / (2.0 * h)
    assert_almost_equal(
        result.deriv.v, expected_deriv, atol=1e-6, msg="d/dx psi(x) = trigamma(x)"
    )


def test_dual_bessel_j0_derivative() raises:
    """`d/dx J0(x) = -J1(x)`."""
    var x0 = 3.0
    var x = variable(x0)
    var result = bessel_j0_fn(x)
    assert_almost_equal(result.value.v, bessel_j0(Scalar[f64](x0)))
    assert_almost_equal(result.deriv.v, -bessel_j1(Scalar[f64](x0)))


def test_dual_bessel_j1_derivative() raises:
    """`d/dx J1(x) = J0(x) - J1(x)/x`."""
    var x0 = 3.0
    var x = variable(x0)
    var result = bessel_j1_fn(x)
    assert_almost_equal(result.value.v, bessel_j1(Scalar[f64](x0)))
    var expected = bessel_j0(Scalar[f64](x0)) - bessel_j1(
        Scalar[f64](x0)
    ) / x0
    assert_almost_equal(result.deriv.v, expected)


def test_dual_bessel_y0_derivative() raises:
    """`d/dx Y0(x) = -Y1(x)`."""
    var x0 = 3.0
    var x = variable(x0)
    var result = bessel_y0_fn(x)
    assert_almost_equal(result.value.v, bessel_y0(Scalar[f64](x0)))
    assert_almost_equal(result.deriv.v, -bessel_y1(Scalar[f64](x0)))


def test_dual_bessel_i1_derivative() raises:
    """`d/dx I1(x) = I0(x) - I1(x)/x`."""
    var x0 = 2.0
    var x = variable(x0)
    var result = bessel_i1_fn(x)
    assert_almost_equal(result.value.v, bessel_i1(Scalar[f64](x0)))
    var expected = bessel_i0(Scalar[f64](x0)) - bessel_i1(
        Scalar[f64](x0)
    ) / x0
    assert_almost_equal(result.deriv.v, expected)


def test_dual_bessel_k1_derivative() raises:
    """`d/dx K1(x) = -K0(x) - K1(x)/x`."""
    var x0 = 2.0
    var x = variable(x0)
    var result = bessel_k1_fn(x)
    assert_almost_equal(result.value.v, bessel_k1(Scalar[f64](x0)))
    var expected = -bessel_k0(Scalar[f64](x0)) - bessel_k1(
        Scalar[f64](x0)
    ) / x0
    assert_almost_equal(result.deriv.v, expected)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
