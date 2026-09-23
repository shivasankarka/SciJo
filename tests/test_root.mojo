"""
Tests for multi-dimensional root-finding (scijo.optimize.root).

To run: `mojo test tests/test_root.mojo -I .` from the project root directory.
"""

# ===----------------------------------------------------------------------=== #
# Stdlib
# ===----------------------------------------------------------------------=== #
from std.python import (
    Python,
    PythonObject,
)
from std.testing import (
    assert_almost_equal,
    assert_equal,
    assert_true,
    TestSuite,
)

# ===----------------------------------------------------------------------=== #
# External
# ===----------------------------------------------------------------------=== #
import numojo as nm

# ===----------------------------------------------------------------------=== #
# SciJo
# ===----------------------------------------------------------------------=== #
import scijo as sj
from scijo.linalg import matmul
from scijo.optimize import root


def test_root_linear_system() raises:
    """Solving a linear system: f(x) = Ax - b, root at x = A^-1 b."""

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        var A = nm.fromstring[dtype]("[[3, 1], [1, 2]]")
        var b = nm.fromstring[dtype]("[9, 8]")
        return matmul(A, x) - b

    var x0 = nm.fromstring[sj.f64]("[0, 0]")
    var result = root[sj.f64, f](x0, atol=1e-10, rtol=1e-12)
    assert_equal(result.success, True)
    assert_almost_equal(result.x.item(0), 2.0, atol=1e-8)
    assert_almost_equal(result.x.item(1), 3.0, atol=1e-8)
    assert_equal(result.method, "newton")


def test_root_nonlinear_system() raises:
    """Classic 2-eqn nonlinear system: x^2 + y^2 = 4, x*y = 1."""

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        var out = nm.zeros[dtype](nm.Shape(2))
        out.store(0, x.item(0) * x.item(0) + x.item(1) * x.item(1) - 4.0)
        out.store(1, x.item(0) * x.item(1) - 1.0)
        return out^

    var x0 = nm.fromstring[sj.f64]("[1.5, 1.0]")
    var result = root[sj.f64, f](x0, atol=1e-10, rtol=1e-12)
    assert_equal(result.success, True)
    # Check the residual, rather than a specific root (there are 4 by symmetry).
    assert_almost_equal(result.fun.item(0), 0.0, atol=1e-6)
    assert_almost_equal(result.fun.item(1), 0.0, atol=1e-6)


def test_root_with_analytic_jacobian() raises:
    """Same linear system, but with an analytic Jacobian supplied."""

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        var A = nm.fromstring[dtype]("[[3, 1], [1, 2]]")
        var b = nm.fromstring[dtype]("[9, 8]")
        return matmul(A, x) - b

    @parameter
    def jac[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        return nm.fromstring[dtype]("[[3, 1], [1, 2]]")

    var x0 = nm.fromstring[sj.f64]("[0, 0]")
    var result = root[sj.f64, f, jac](x0, atol=1e-10, rtol=1e-12)
    assert_equal(result.success, True)
    assert_almost_equal(result.x.item(0), 2.0, atol=1e-8)
    assert_almost_equal(result.x.item(1), 3.0, atol=1e-8)
    # An analytic Jacobian converges in a single Newton step on a linear system.
    assert_equal(result.nit, 1)


def test_root_broyden_linear_system() raises:
    """Broyden's method on a linear system: same problem and expected root
    as `test_root_linear_system`, without ever forming a Jacobian.
    """

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        var A = nm.fromstring[dtype]("[[3, 1], [1, 2]]")
        var b = nm.fromstring[dtype]("[9, 8]")
        return matmul(A, x) - b

    var x0 = nm.fromstring[sj.f64]("[0, 0]")
    var result = root[sj.f64, f, method="broyden1"](x0, atol=1e-10, rtol=1e-12)
    assert_equal(result.success, True)
    assert_almost_equal(result.x.item(0), 2.0, atol=1e-6)
    assert_almost_equal(result.x.item(1), 3.0, atol=1e-6)
    assert_equal(result.method, "broyden1")


def test_root_broyden_nonlinear_system() raises:
    """Broyden's method on the same nonlinear system as
    `test_root_nonlinear_system`: check the residual, not a specific root
    (there are 4 by symmetry).
    """

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        var out = nm.zeros[dtype](nm.Shape(2))
        out.store(0, x.item(0) * x.item(0) + x.item(1) * x.item(1) - 4.0)
        out.store(1, x.item(0) * x.item(1) - 1.0)
        return out^

    var x0 = nm.fromstring[sj.f64]("[1.5, 1.0]")
    var result = root[sj.f64, f, method="broyden1"](
        x0, atol=1e-10, rtol=1e-12, maxiter=200
    )
    assert_equal(result.success, True)
    assert_almost_equal(result.fun.item(0), 0.0, atol=1e-6)
    assert_almost_equal(result.fun.item(1), 0.0, atol=1e-6)


def test_root_broyden_matches_scipy() raises:
    """Broyden's method converges to the same root as
    `scipy.optimize.root(method="broyden1")`.

    Uses the linear system from `test_root_linear_system`/
    `test_root_broyden_linear_system` (a single, globally unique root)
    rather than a nonlinear system, since a nonlinear system can have
    several valid roots and scipy's `broyden1` - which damps its steps and
    scales its initial Jacobian estimate differently from the plain
    textbook update here - is not guaranteed to walk toward the same one
    from a shared starting point even when both solvers converge correctly.
    """
    var opt: PythonObject
    try:
        opt = Python.import_module("scipy.optimize")
    except:
        print("SciPy not available, skipping root broyden1 comparison test")
        return

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        var A = nm.fromstring[dtype]("[[3, 1], [1, 2]]")
        var b = nm.fromstring[dtype]("[9, 8]")
        return matmul(A, x) - b

    var x0 = nm.fromstring[sj.f64]("[0, 0]")
    var result = root[sj.f64, f, method="broyden1"](x0, atol=1e-10, rtol=1e-12)
    assert_equal(result.success, True)

    # scipy.optimize.root needs a Python-callable objective; a Mojo closure
    # can't be passed directly as one, so build the (numpy-only) equivalent
    # from a Python source string instead.
    var py_f = Python.evaluate(
        "lambda x: [3*x[0] + x[1] - 9, x[0] + 2*x[1] - 8]"
    )
    var py_x0 = Python.list(0.0, 0.0)
    var py_result = opt.root(py_f, py_x0, method="broyden1")
    assert_true(
        Bool(py=py_result.success), "scipy broyden1 should itself converge"
    )
    var py_x = py_result.x

    assert_almost_equal(
        result.x.item(0),
        Float64(py=py_x[0]),
        atol=1e-4,
        msg="Broyden root x0 should match scipy",
    )
    assert_almost_equal(
        result.x.item(1),
        Float64(py=py_x[1]),
        atol=1e-4,
        msg="Broyden root x1 should match scipy",
    )


def test_root_broyden_shape_mismatch_raises() raises:
    """F(x) with a different length than x should raise for Broyden too."""

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        return nm.zeros[dtype](nm.Shape(3))

    var x0 = nm.fromstring[sj.f64]("[0, 0]")
    var raised = False
    try:
        _ = root[sj.f64, f, method="broyden1"](x0)
    except:
        raised = True
    assert_equal(raised, True)


def test_root_shape_mismatch_raises() raises:
    """F(x) with a different length than x should raise, not misbehave."""

    @parameter
    def f[
        dtype: DType
    ](
        x: nm.NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> nm.NDArray[dtype]:
        return nm.zeros[dtype](nm.Shape(3))

    var x0 = nm.fromstring[sj.f64]("[0, 0]")
    var raised = False
    try:
        _ = root[sj.f64, f](x0)
    except:
        raised = True
    assert_equal(raised, True)


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
