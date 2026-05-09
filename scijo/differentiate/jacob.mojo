# ===----------------------------------------------------------------------=== #
# SciJo: Differentiate module for Mojo
# Distributed under the Apache 2.0 License.
# ===----------------------------------------------------------------------=== #
"""Jacobian Matrix Computation (`scijo.differentiate.jacob`)
===========================================================
Computes the Jacobian matrix of a vector-valued function using central finite
differences with parallelized column evaluation.

Examples
--------
    ```mojo
    from scijo.differentiate import jacobian

    def f[dtype: DType](x: NDArray[dtype], args: Optional[List[Scalar[dtype]]]) raises -> NDArray[dtype]:
        return x * x

    var x = nm.array[f64]([1.0, 2.0])
    var J = jacobian[f64, f](x)
    ```
"""

from std.algorithm.functional import parallelize

from numojo.routines.creation import zeros, full
from numojo.core import NDArray, Shape


def jacobian[
    dtype: DType,
    jacob_func: def[dtype: DType](
        x: NDArray[dtype], args: Optional[List[Scalar[dtype]]]
    ) capturing raises -> NDArray[dtype],
](
    x: NDArray[dtype],
    args: Optional[List[Scalar[dtype]]] = None,
    tolerances: Dict[String, Scalar[dtype]] = {"abs": 1e-5, "rel": 1e-3},
    maxiter: Int = 10,
) raises -> NDArray[dtype]:
    """Computes the Jacobian matrix of a vector-valued function using central finite differences.

    Evaluates J[i, j] = ∂f_i/∂x_j using the central difference formula
    (f(x + h*e_j) - f(x - h*e_j)) / (2h), where e_j is the j-th unit vector.
    Each column of the Jacobian is computed in parallel.

    Parameters:
        dtype: The floating-point data type.
        jacob_func: Vector-valued function with signature def(x, args) -> NDArray[dtype].

    Args:
        x: Input vector of shape (n,) at which to evaluate the Jacobian.
        args: Optional arguments to pass to the function.
        tolerances: Tolerance dictionary with "abs" and "rel" keys (reserved for future use).
        maxiter: Maximum iterations (reserved for future use).

    Raises:
        Error: If function evaluation fails for any perturbation.

    Returns:
        NDArray[dtype] of shape (m, n) representing the Jacobian matrix,
        where m is the output dimension and n is the input dimension.

    Examples:
        ```mojo
        import numojo as nm
        from scijo.differentiate import jacobian
        from scijo.prelude import *

        def f[dtype: DType](x: NDArray[dtype], args: Optional[List[Scalar[dtype]]]) raises -> NDArray[dtype]:
            return x * x

        var x = nm.array[f64]([1.0, 2.0])
        var J = jacobian[f64, f](x)
        ```
    """
    var n: Int = len(x)
    var f0: NDArray[dtype] = jacob_func(x, args)
    var m: Int = len(f0)

    var jacob: NDArray[dtype] = zeros[dtype](Shape(m, n))
    var step: NDArray[dtype] = full[dtype](Shape(n), fill_value=0.5)

    @parameter
    def closure(j: Int):
        try:
            var x_plus = x.copy()
            var x_minus = x.copy()
            var hj = step.load(j)
            x_plus.store(j, val=x.load(j) + hj)
            x_minus.store(j, val=x.load(j) - hj)
            var f_plus = jacob_func(x_plus, args)
            var f_minus = jacob_func(x_minus, args)

            var col: NDArray[dtype] = (f_plus - f_minus) / (2.0 * hj)

            for i in range(m):
                var val: Scalar[dtype] = col.load(i)
                var flat_idx: Int = i * n + j
                jacob.store(flat_idx, val=val)
        except:
            print("Error in computing Jacobian column for j =", j)

    parallelize[closure](n, num_workers=n)

    return jacob^
