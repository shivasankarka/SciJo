fn _binary_search[
    dtype: DType
](x: NDArray[dtype], value: Scalar[dtype]) raises -> Int:
    """
    Binary search to find the interval containing the interpolation point.

    Args:
        x: Sorted array of x-coordinates.
        value: The value to search for.

    Returns:
        Index j such that x[j-1] <= value < x[j], or appropriate boundary index.
    """
    var left: Int = 0
    var right: Int = x.size - 1

    while right - left > 1:
        var mid = (left + right) // 2
        if x._buf.ptr[mid] <= value:
            left = mid
        else:
            right = mid

    return right


fn _validate_interpolation_input[
    dtype: DType
](x: NDArray[dtype], y: NDArray[dtype]) raises:
    """
    Validate input arrays for interpolation.

    Args:
        x: Array of x-coordinates.
        y: Array of y-coordinates.

    Raises:
        Error: If input arrays are invalid.
    """
    if x.size != y.size:
        raise Error("x and y arrays must have the same length")

    if x.size < 2:
        raise Error("x and y arrays must have at least 2 points")

    for i in range(1, x.size):
        if x.item(i) <= x.item(i - 1):
            raise Error("x array must be strictly increasing")
