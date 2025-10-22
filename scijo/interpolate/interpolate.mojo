"""
===----------------------------------------------------------------------===
Interpolate Module - Implements interpolation functions
Last updated: 2025-10-14
===----------------------------------------------------------------------===
"""

from numojo import zeros

from .utility import _binary_search, _validate_interpolation_input

# ! In a lot of places, I am using pure pointer method to retrieve values from NDArray, not so safe. Also causes some problems especially with .item() method of NDArray.


# TODO: Add extrapolation and fill_value handling to LinearInterpolator
# ! Array needs to be sorted for binary search to work
struct LinearInterpolator[dtype: DType = DType.float64](Copyable, Movable):
    """
    A callable linear interpolation object similar to scipy.interpolate.interp1d.

    This struct stores the interpolation data (x, y) and provides a callable interface
    that can interpolate single values or arrays of values efficiently using binary search.

    The interpolator performs linear interpolation between adjacent data points.
    For points outside the data range, behavior is controlled by bounds_error and fill_value.

    Example:
        ```mojo
        import scijo as sj
        from scijo.interpolate import LinearInterpolator
        import numojo as nm

        fn main() raises:
            var x = nm.arange[sj.f64](0, 10, 1)  # [0, 1, 2, ..., 9]
            var y = x * x  # quadratic function: [0, 1, 4, 9, ..., 81]
            var f = LinearInterpolator(x, y)

            # Interpolate single value
            var result = f(3.5)  # returns 12.25
            print("f(3.5) =", result)

            # Interpolate array of values
            var xi = nm.arange[sj.f64](0.5, 9.5, 0.5)
            var yi = f(xi)  # returns array of interpolated values
            print("Interpolated values:", yi)
        ```
    """

    var x: NDArray[dtype]
    var y: NDArray[dtype]
    var bounds_error: Bool
    var fill_value: Optional[Scalar[dtype]]

    fn __init__(
        out self,
        x: NDArray[dtype],
        y: NDArray[dtype],
        bounds_error: Bool = True,
        fill_value: Optional[Scalar[dtype]] = None,
    ) raises:
        """
        Initialize the linear interpolator.

        Args:
            x: The x-coordinates of the data points, must be strictly increasing.
            y: The y-coordinates of the data points, same length as x.
            bounds_error: If True, raise error when interpolating outside bounds.
                          If False, use fill_value or extrapolate linearly.
            fill_value: Value to use for points outside the data range when
                       bounds_error is False. If None, extrapolate linearly.

        Raises:
            Error: If x and y have different lengths, have fewer than 2 points,
                   or x is not strictly increasing.
        """
        _validate_interpolation_input(x, y)

        self.x = x.copy()
        self.y = y.copy()
        self.bounds_error = bounds_error
        self.fill_value = fill_value

    fn __call__(self, xi: Scalar[dtype]) raises -> Scalar[dtype]:
        """
        Interpolate a single value.

        Args:
            xi: The point at which to interpolate.

        Returns:
            The interpolated value at xi.

        Raises:
            Error: If bounds_error is True and xi is outside data range.
        """

        var x_min = self.x._buf.ptr[0]
        var x_max = self.x._buf.ptr[self.x.size - 1]

        if xi < x_min or xi > x_max:
            if self.bounds_error:
                raise Error(
                    "Interpolation point "
                    + String(xi)
                    + " is outside data range ["
                    + String(x_min)
                    + ", "
                    + String(x_max)
                    + "]"
                )
            elif self.fill_value:
                return self.fill_value.value()

        if xi <= x_min:
            if self.fill_value:
                return self.fill_value.value()
            return self.y._buf.ptr[0]
        elif xi >= x_max:
            if self.fill_value:
                return self.fill_value.value()
            return self.y._buf.ptr[self.y.size - 1]

        var j: Int = _binary_search(self.x, xi)

        var x0: Scalar[dtype] = self.x._buf.ptr[j - 1]
        var x1: Scalar[dtype] = self.x._buf.ptr[j]
        var y0: Scalar[dtype] = self.y._buf.ptr[j - 1]
        var y1: Scalar[dtype] = self.y._buf.ptr[j]

        var slope: Scalar[dtype] = (y1 - y0) / (x1 - x0)
        return y0 + slope * (xi - x0)

    fn __call__(self, xi: NDArray[dtype]) raises -> NDArray[dtype]:
        """
        Interpolate an array of values.

        Args:
            xi: Array of points at which to interpolate.

        Returns:
            Array of interpolated values with the same shape as xi.

        Raises:
            Error: If bounds_error is True and any point in xi is outside data range.
        """
        var result: NDArray[dtype] = zeros[dtype](xi.shape)
        var x_min: Scalar[dtype] = self.x._buf.ptr[0]
        var x_max: Scalar[dtype] = self.x._buf.ptr[self.x.size - 1]

        for i in range(xi.size):
            var x_val: Scalar[dtype] = xi._buf.ptr[i]

            if x_val < x_min or x_val > x_max:
                if self.bounds_error:
                    raise Error(
                        "Interpolation point "
                        + String(x_val)
                        + " is outside data range ["
                        + String(x_min)
                        + ", "
                        + String(x_max)
                        + "]"
                    )
                elif self.fill_value:
                    result._buf.ptr[i] = self.fill_value.value()
                    continue

            if x_val <= x_min:
                if self.fill_value and x_val < x_min:
                    result._buf.ptr[i] = self.fill_value.value()
                else:
                    result._buf.ptr[i] = self.y._buf.ptr[0]
                continue
            elif x_val >= x_max:
                if self.fill_value and x_val > x_max:
                    result._buf.ptr[i] = self.fill_value.value()
                else:
                    result._buf.ptr[i] = self.y._buf.ptr[self.y.size - 1]
                continue

            var j: Int = _binary_search(self.x, x_val)

            var x0: Scalar[dtype] = self.x._buf.ptr[j - 1]
            var x1: Scalar[dtype] = self.x._buf.ptr[j]
            var y0: Scalar[dtype] = self.y._buf.ptr[j - 1]
            var y1: Scalar[dtype] = self.y._buf.ptr[j]

            var slope: Scalar[dtype] = (y1 - y0) / (x1 - x0)
            result._buf.ptr[i] = y0 + slope * (x_val - x0)

        return result^


# TODO: Add more interpolation methods like 'quadratic', 'cubic'.
# TODO: Add both interpolate and extrapolate fill methods.
fn interp1d[
    dtype: DType = DType.float64
](
    x: NDArray[dtype],
    y: NDArray[dtype],
    bounds_error: Bool = True,
    fill_value: Optional[Scalar[dtype]] = None,
) raises -> LinearInterpolator[dtype]:
    """
    Interpolates the values of y at new points using linear interpolation.

    Parameters:
        dtype: The element type (default: DType.float64).

    Args:
        x: The x-coordinates of the data points, must be strictly increasing.
        y: The y-coordinates of the data points, same length as x.
        bounds_error: If True, raise error when interpolating outside bounds.
                      If False, use fill_value or extrapolate linearly.
        fill_value: Value to use for points outside the data range when
                   bounds_error is False. If None, extrapolate linearly.

    Returns:
        A callable LinearInterpolator object.

    Raises:
        Error: If x and y have different lengths, have fewer than 2 points,
               or x is not strictly increasing.

    Example:
        ```mojo
        import scijo as sj
        from scijo.interpolate import interp1d
        import numojo as nm

        fn main() raises:
            # Create data points
            var x = nm.arange[sj.f64](0, 10, 1)  # [0, 1, 2, ..., 9]
            var y = x * x  # [0, 1, 4, 9, ..., 81]

            # Create interpolation function
            var f = interp1d(x, y)

            # Interpolate single value
            var result = f(3.5)  # Should be approximately 12.25
            print("f(3.5) =", result)

            # Interpolate array of values
            var xi = nm.arange[sj.f64](0.5, 9.5, 0.5)
            var yi = f(xi)
            print("Interpolated values:", yi)

            # Create interpolator with bounds handling
            var f_fill = interp1d(x, y, bounds_error=False, fill_value=Scalar[sj.f64](-999.0))
            var out_of_bounds = f_fill(15.0)  # Returns -999.0
            print("Out of bounds value:", out_of_bounds)
        ```
    """
    return LinearInterpolator[dtype](x, y, bounds_error, fill_value)


fn interp1d[
    dtype: DType = DType.float64,
    type: String = "linear",
    fill_method: String = "interpolate",
](
    xi: NDArray[dtype],
    x: NDArray[dtype],
    y: NDArray[dtype],
) raises -> NDArray[
    dtype
]:
    """
    Interpolate the values of y at the points xi using specified method. Similar to numpy interp.

    This is a functional interface that directly returns interpolated values
    without creating a reusable interpolator object.

    Parameters:
        dtype: The element type (default: DType.float64).
        type: The interpolation method ("linear" is currently supported).
        fill_method: How to handle out-of-bounds values:
                    - "interpolate": Clamp to boundary values
                    - "extrapolate": Linearly extrapolate beyond boundaries

    Args:
        xi: Array of points at which to interpolate.
        x: Array of x-coordinates of data points, must be strictly increasing.
        y: Array of y-coordinates of data points, same length as x.

    Returns:
        The interpolated values of y at the points xi as an NDArray of dtype.

    Raises:
        Error: If x and y have different lengths, have fewer than 2 points,
               x is not strictly increasing, or invalid method/fill_method specified.

    Example:
        ```mojo
        import scijo as sj
        from scijo.interpolate import interp1d
        import numojo as nm

        fn main() raises:
            var x = nm.arange[sj.f64](0, 5, 1)    # [0, 1, 2, 3, 4]
            var y = x * x                         # [0, 1, 4, 9, 16]
            var xi = nm.arange[sj.f64](0.5, 4.0, 0.5)  # [0.5, 1.5, 2.5, 3.5]

            # Interpolate with clamping at boundaries
            var yi = interp1d(xi, x, y, type="linear", fill_method="interpolate")

            # Interpolate with extrapolation beyond boundaries
            var ye = interp1d(xi, x, y, type="linear", fill_method="extrapolate")
        ```
    """
    _validate_interpolation_input(x, y)

    @parameter
    if type == "linear" and fill_method == "extrapolate":
        return _interp1d_linear_extrapolate(xi, x, y)
    elif type == "linear" and fill_method == "interpolate":
        return _interp1d_linear_interpolate(xi, x, y)
    else:
        raise Error(
            String(
                "Invalid interpolation method: {} with fill_method: {}."
                " Supported: type='linear' with fill_method='interpolate' or"
                " 'extrapolate'"
            ).format(type, fill_method)
        )


fn _interp1d_linear_interpolate[
    dtype: DType
](xi: NDArray[dtype], x: NDArray[dtype], y: NDArray[dtype]) raises -> NDArray[
    dtype
]:
    """
    Linear interpolation with boundary clamping.

    For points outside the data range, returns the boundary values (y[0] or y[-1]).

    Parameters:
        dtype: The element type.

    Args:
        xi: Array of interpolation points.
        x: Array of x-coordinates (must be sorted).
        y: Array of y-coordinates.

    Returns:
        Array of interpolated values.
    """
    var result: NDArray[dtype] = NDArray[dtype](xi.shape)
    var x_min: Scalar[dtype] = x._buf.ptr[0]
    var x_max: Scalar[dtype] = x._buf.ptr[x.size - 1]

    for i in range(xi.size):
        var xi_val: Scalar[dtype] = xi._buf.ptr[i]

        if xi_val <= x_min:
            result.itemset(i, y._buf.ptr[0])
        elif xi_val >= x_max:
            result.itemset(i, y._buf.ptr[y.size - 1])
        else:
            var j: Int = _binary_search(x, xi_val)

            var x0: Scalar[dtype] = x._buf.ptr[j - 1]
            var x1: Scalar[dtype] = x._buf.ptr[j]
            var y0: Scalar[dtype] = y._buf.ptr[j - 1]
            var y1: Scalar[dtype] = y._buf.ptr[j]
            var t: Scalar[dtype] = (xi_val - x0) / (x1 - x0)
            result._buf.ptr[i] = y0 + t * (y1 - y0)

    return result^


fn _interp1d_linear_extrapolate[
    dtype: DType
](xi: NDArray[dtype], x: NDArray[dtype], y: NDArray[dtype]) raises -> NDArray[
    dtype
]:
    """
    Linear interpolation with linear extrapolation beyond boundaries.

    For points outside the data range, extrapolates using the slope of the
    nearest boundary segment.

    Parameters:
        dtype: The element type.

    Args:
        xi: Array of interpolation points.
        x: Array of x-coordinates (must be sorted).
        y: Array of y-coordinates.

    Returns:
        Array of interpolated/extrapolated values.
    """
    var result: NDArray[dtype] = NDArray[dtype](xi.shape)
    var x_min: Scalar[dtype] = x._buf.ptr[0]
    var x_max: Scalar[dtype] = x._buf.ptr[x.size - 1]

    for i in range(xi.size):
        var xi_val: Scalar[dtype] = xi._buf.ptr[i]

        if xi_val < x_min:
            var slope = (y._buf.ptr[1] - y._buf.ptr[0]) / (
                x._buf.ptr[1] - x._buf.ptr[0]
            )
            result.itemset(i, y._buf.ptr[0] + slope * (xi_val - x._buf.ptr[0]))
        elif xi_val > x_max:
            var slope = (y._buf.ptr[y.size - 1] - y._buf.ptr[y.size - 2]) / (
                x._buf.ptr[x.size - 1] - x._buf.ptr[x.size - 2]
            )
            result.itemset(
                i,
                y._buf.ptr[y.size - 1]
                + slope * (xi_val - x._buf.ptr[x.size - 1]),
            )
        else:
            if xi_val == x_min:
                result.itemset(i, y._buf.ptr[0])
            elif xi_val == x_max:
                result.itemset(i, y._buf.ptr[y.size - 1])
            else:
                var j: Int = _binary_search(x, xi_val)

                var x0: Scalar[dtype] = x._buf.ptr[j - 1]
                var x1: Scalar[dtype] = x._buf.ptr[j]
                var y0: Scalar[dtype] = y._buf.ptr[j - 1]
                var y1: Scalar[dtype] = y._buf.ptr[j]
                var t: Scalar[dtype] = (xi_val - x0) / (x1 - x0)
                result._buf.ptr[i] = y0 + t * (y1 - y0)

    return result^


#  Higher-order interpolation methods

# fn _interp1d_quadratic_interpolate[dtype: DType](
#     xi: NDArray[dtype], x: NDArray[dtype], y: NDArray[dtype]
# ) raises -> NDArray[dtype]:
#     """Quadratic interpolation with boundary clamping."""
#     # Implementation would use 3-point Lagrange interpolation
#     pass

# fn _interp1d_cubic_interpolate[dtype: DType](
#     xi: NDArray[dtype], x: NDArray[dtype], y: NDArray[dtype]
# ) raises -> NDArray[dtype]:
#     """Cubic interpolation with boundary clamping."""
#     # Implementation would use 4-point Lagrange interpolation or cubic splines
#     pass
