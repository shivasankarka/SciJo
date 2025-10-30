from numojo.core.ndarray import NDArray, NDArrayShape
from numojo.core.error import *
import numojo as nm


fn trapezoid[
    dtype: DType
](
    y: NDArray[dtype],
    dx: Scalar[dtype] = 1.0,
    axis: Int = -1,
) raises -> Scalar[
    dtype
]:
    """
    Integrate along the given axis using the composite trapezoidal rule.
    Integrates y(x) along each 1d slice on the given axis, computing ∫ y(x) dx. Uses evenly spaced points with spacing dx.

    Parameters:
        dtype: The data type of the input arrays and the output scalar.
               Must be a floating-point type.

    Arguments:
        y: Input array to integrate. Must be 1-D for this implementation.
        dx: The spacing between sample points. The default is 1.0.
        axis: The axis along which to integrate. Currently only supports 1-D arrays,
              so this parameter is ignored.

    Returns:
        Scalar[dtype]: Definite integral of y as approximated by the trapezoidal rule. Returns 0.0 for arrays with fewer than 2 elements.

    Raises:
        Error(ShapeError): If y is not 1-D.
        Error(ValueError): If y is empty.

    Examples:
        >>> var y = nm.fromstring[DType.float64]("[1, 2, 3]")
        >>> var result = trapezoid[DType.float64](y)  # Returns 4.0
        >>> var result = trapezoid[DType.float64](y, dx=2.0)  # Returns 8.0.
    """
    constrained[
        dtype.is_floating_point(),
        (
            "trapezoid[dtype: DType](y, dx, axis): dtype must be a"
            " floating-point type."
        ),
    ]()

    if y.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected y to be 1-D, received ndim={}."
                ).format(y.ndim),
                suggestion="Pass a 1-D NDArray for y (e.g. shape (N,)).",
                location="trapezoid(y, dx=1.0)",
            )
        )

    if y.size == 0:
        raise Error(
            ValueError(
                message="Cannot integrate over an empty array.",
                suggestion="Provide a non-empty array for y.",
                location="trapezoid(y, dx=1.0)",
            )
        )

    if y.size == 1:
        return Scalar[dtype](0.0)

    var integral: Scalar[dtype] = 0.0
    for i in range(y.size - 1):
        var y_i = y.item(i)
        var y_i1 = y.item(i + 1)
        integral += (y_i + y_i1) * dx * 0.5

    return integral


fn trapezoid[
    dtype: DType
](
    y: NDArray[dtype],
    x: NDArray[dtype],
    axis: Int = -1,
) raises -> Scalar[
    dtype
]:
    """
    Integrate along the given axis using the composite trapezoidal rule. Integrates y(x) along each 1d slice on the given axis, computing ∫ y(x) dx. When x is specified, this integrates along the parametric curve.

    Parameters:
        dtype: The data type of the input arrays and the output scalar.
               Must be a floating-point type.

    Arguments:
        y: Input array to integrate. Must be 1-D for this implementation.
        x: Array of sample points corresponding to the y values.
        axis: The axis along which to integrate. Currently only supports 1-D arrays.

    Returns:
        Scalar[dtype]: Definite integral of y as approximated by the trapezoidal rule. Returns 0.0 for arrays with fewer than 2 elements.

    Raises:
        Error(ShapeError): If y or x are not 1-D, or if their sizes differ.
        Error(ValueError): If y is empty.

    Examples:
        >>> var y = nm.fromstring[DType.float64]("[1, 2, 3]")
        >>> var x = nm.fromstring[DType.float64]("[4, 6, 8]")
        >>> var result = trapezoid[DType.float64](y, x)  # Returns 8.0.
    """
    constrained[
        dtype.is_floating_point(),
        (
            "trapezoid[dtype: DType](y, x, axis): dtype must be a"
            " floating-point type."
        ),
    ]()

    if y.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected y to be 1-D, received ndim={}."
                ).format(y.ndim),
                suggestion="Pass a 1-D NDArray for y (e.g. shape (N,)).",
                location="trapezoid(y, x)",
            )
        )

    if y.size == 0:
        raise Error(
            ValueError(
                message="Cannot integrate over an empty array.",
                suggestion="Provide a non-empty array for y.",
                location="trapezoid(y, x)",
            )
        )

    if y.size == 1:
        return Scalar[dtype](0.0)

    if x.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected x to be 1-D, received ndim={}."
                ).format(x.ndim),
                suggestion="Provide a 1-D NDArray for x.",
                location="trapezoid(y, x)",
            )
        )

    if y.size != x.size:
        raise Error(
            ShapeError(
                message=(
                    String("Size mismatch: y.size={} != x.size={}.").format(
                        y.size, x.size
                    )
                ),
                suggestion="Ensure x and y have identical lengths.",
                location="trapezoid(y, x)",
            )
        )

    var integral: Scalar[dtype] = 0.0
    for i in range(y.size - 1):
        var y_i = y.item(i)
        var y_i1 = y.item(i + 1)
        var x_i = x.item(i)
        var x_i1 = x.item(i + 1)
        var dx_segment = x_i1 - x_i
        integral += (y_i + y_i1) * dx_segment * 0.5

    return integral

fn simpson[
    dtype: DType
](
    y: NDArray[dtype],
    dx: Scalar[dtype] = 1.0,
    axis: Int = -1,
) raises -> Scalar[dtype]:
    """
    Integrate along the given axis using Simpson's rule. Integrates y(x) along each
    1d slice on the given axis. Uses evenly spaced points with spacing dx.

    Parameters:
        dtype: The data type of the input arrays and the output scalar.
               Must be a floating-point type.

    Args:
        y: Input array to integrate. Must be 1-D for this implementation.
        dx: The spacing between sample points. The default is 1.0.
        axis: The axis along which to integrate. Currently only supports 1-D arrays,
              so this parameter is ignored.

    Returns:
        Scalar[dtype]: Definite integral of y as approximated by Simpson's rule.
    """
    if y.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected y to be 1-D, received ndim={}."
                ).format(y.ndim),
                suggestion="Pass a 1-D NDArray for y (e.g. shape (N,)). Only 1-D arrays are supported currently.",
                location="simpson(y, dx=1.0)",
            )
        )
    var integral: Scalar[dtype] = 0.0
    alias multiplier: Scalar[dtype] = 1.0 / 6.0
    for i in range(0, y.size - 1, 2):
        integral +=  multiplier * dx * 2 * (
            y.item(i) + 4.0 * y.item(i + 1) + y.item(i + 2)
        )

    return integral

fn simpson[
    dtype: DType
](
    y: NDArray[dtype],
    x: NDArray[dtype],
    axis: Int = -1,
) raises -> Scalar[dtype]:
    """
    Integrate along the given axis using Simpson's rule. Integrates y(x) along each
    1d slice on the given axis. When x is specified, this integrates along the
    parametric curve.

    Parameters:
        dtype: The data type of the input arrays and the output scalar.
               Must be a floating-point type.

    Args:
        y: Input array to integrate. Must be 1-D for this implementation.
        x: Array of sample points corresponding to the y values.
        axis: The axis along which to integrate. Currently only supports 1-D arrays.

    Returns:
        Scalar[dtype]: Definite integral of y as approximated by Simpson's rule.
    """
    if y.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected y to be 1-D, received ndim={}."
                ).format(y.ndim),
                suggestion="Pass a 1-D NDArray for y (e.g. shape (N,)). Only 1-D arrays are supported currently.",
                location="simpson(y, x)",
            )
        )

    if x.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected x to be 1-D, received ndim={}."
                ).format(x.ndim),
                suggestion="Provide a 1-D NDArray for x.",
                location="simpson(y, x)",
            )
        )

    if y.size != x.size:
        raise Error(
            ShapeError(
                message=(
                    String("Size mismatch: y.size={} != x.size={}.").format(
                        y.size, x.size
                    )
                ),
                suggestion="Ensure x and y have identical lengths.",
                location="simpson(y, x)",
            )
        )

    var integral: Scalar[dtype] = 0.0
    alias multiplier: Scalar[dtype] = 1.0 / 6.0
    for i in range(1, y.size - 1, 2):
        var dx_segment = x.item(i + 1) - x.item(i - 1)
        integral += multiplier * dx_segment * (
            y.item(i - 1) + 4.0 * y.item(i) + y.item(i + 1)
            )

    if y.size % 2 == 0:
        var y_n1 = y.item(y.size - 2)
        var y_n = y.item(y.size - 1)
        var x_n1 = x.item(x.size - 2)
        var x_n = x.item(x.size - 1)
        var dx_last = x_n - x_n1
        integral += (y_n1 + y_n) * dx_last * 0.5

    return integral

# TODO: fix the loop implementation.
fn romb[dtype: DType](y: NDArray[dtype], dx: Scalar[dtype] = 1.0, axis: Int = -1) raises -> Scalar[dtype]:
    """
    Integrate along the given axis using Romberg integration. Integrates y(x) along each 1d slice on the given axis, computing ∫ y(x) dx. Uses evenly spaced points with spacing dx.

    Parameters:
        dtype: The data type of the input arrays and the output scalar.
               Must be a floating-point type.

    Args:
        y: Input array to integrate. Must be 1-D for this implementation.
        dx: The spacing between sample points. The default is 1.0.
        axis: The axis along which to integrate. Currently only supports 1-D arrays,
              so this parameter is ignored.

    Returns:
        Scalar[dtype]: Definite integral of y as approximated by Romberg integration.
    """
    var maxiter: Int = 10
    if y.ndim != 1:
        raise Error(
            ShapeError(
                message=String(
                    "Expected y to be 1-D, received ndim={}."
                ).format(y.ndim),
                suggestion="Pass a 1-D NDArray for y (e.g. shape (N,)). Only 1-D arrays are supported currently.",
                location="romb(y, dx=1.0)",
            )
        )

    var step: Scalar[dtype] = dx
    # var Rone: List[Scalar[dtype]] = List[Scalar[dtype]](capacity= maxiter)
    # var Rtwo: List[Scalar[dtype]] = List[Scalar[dtype]](capacity= maxiter)
    var Rone: NDArray[dtype] = nm.zeros[dtype](NDArrayShape(maxiter))
    var Rtwo: NDArray[dtype] = nm.zeros[dtype](NDArrayShape(maxiter))

    var R1: UnsafePointer[Scalar[dtype]] = Rone.unsafe_ptr()
    var R2: UnsafePointer[Scalar[dtype]] = Rtwo.unsafe_ptr()

    R1[0] = 0.5 * dx * (y.item(0) + y.item(y.size - 1))

    for i in range(1, maxiter):
        step /= 2.0
        var c: Scalar[dtype] = 0
        var ep: Int = 2 * (i - 1)
        for j in range(1, ep + 1):
            c += y.item(Int(2 * j - 1))
        R2[0] = step * c + 0.5 * R1[0]

        for j in range(1, i + 1):
            var const: Scalar[dtype] = Scalar[dtype](4.0)**j
            R2[j] = (const  * R2[j - 1] - R1[j - 1]) / (const - 1.0)

        if i > 1 and abs(R1[i-1] - R2[i]) < Scalar[dtype](1e-6):
            return R2[i]

        var temp: UnsafePointer[Scalar[dtype]] = R1
        R1 = R2
        R2 = temp

    print("Initial R1: ", Rone)
    print("Initial R2: ", Rtwo)

    return Rone.item(maxiter - 1)
