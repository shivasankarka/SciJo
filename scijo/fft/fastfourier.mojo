from numojo.core.complex import ComplexNDArray, ComplexSIMD, ComplexDType
from numojo.core.ndarray import NDArray
from numojo.core.ndshape import NDArrayShape
from numojo.routines.constants import Constants
from numojo.core.item import Item

from math import sin, cos
from algorithm.functional import vectorize
from sys import simd_width_of


# Cooley-Tukey FFT algorithm
fn fft[
    dtype: ComplexDType = ComplexDType.float64
](arr: ComplexNDArray[dtype]) raises -> ComplexNDArray[dtype]:
    """Computes the Fast Fourier Transform using the Cooley-Tukey algorithm.

    The FFT decomposes the DFT computation by recursively breaking down the
    transform of size N into two transforms of size N/2, making it significantly
    faster than the naive O(N²) DFT computation for large arrays.

    Parameters:
        dtype: The data type of the complex elements (ComplexDType).

    Args:
        arr: Input complex array to transform. Must be 1-dimensional with length
             that is a power of 2.

    Returns:
        ComplexNDArray containing the Fast Fourier Transform of the input array.
        The output has the same shape and dtype as the input.

    Raises:
        Error: If the input array is not 1-dimensional.
        Error: If the array length is not a power of 2.

    Example:
        ```mojo
        import scijo as sj
        import numojo as nm
        var input = nm.routines.creation.ones[nm.cf32](nm.Shape(8))
        var transformed = sj.fft.fft[nm.cf32](input)
        ```

    Note:
        This implementation currently only supports 1D arrays with power-of-2 lengths.
        For general-purpose FFT computations, consider using more sophisticated
        algorithms that handle arbitrary sizes.
    """
    if arr.ndim != 1:
        raise Error("FFT currently only supports 1D arrays")

    var n: Int = arr.shape[0]
    if n <= 1:
        return ComplexNDArray[dtype](re=arr._re.copy(), im=arr._im.copy())

    if (n & (n - 1)) != 0:
        raise Error(
            "FFT currently only supports arrays with length that is a power"
            " of 2"
        )

    var half_size = n // 2
    var even_indices = ComplexNDArray[dtype](NDArrayShape(half_size))
    var odd_indices = ComplexNDArray[dtype](NDArrayShape(half_size))

    for i in range(half_size):
        even_indices[Item(i)] = arr[Item(2 * i)]
        odd_indices[Item(i)] = arr[Item(2 * i + 1)]

    var even_fft = fft[dtype](even_indices)
    var odd_fft = fft[dtype](odd_indices)

    var result = ComplexNDArray[dtype](arr.shape)

    for k in range(half_size):
        var angle = (
            -2.0
            * Constants.pi
            * Scalar[dtype._dtype](k)
            / Scalar[dtype._dtype](n)
        )
        var twiddle = ComplexSIMD[dtype](
            cos(angle).cast[dtype._dtype](), sin(angle).cast[dtype._dtype]()
        )

        var twiddle_odd = twiddle * odd_fft[Item(k)]

        result[Item(k)] = even_fft[Item(k)] + twiddle_odd
        result[Item(k + half_size)] = even_fft[Item(k)] - twiddle_odd

    return result^


fn _ifft_unnormalized[
    dtype: ComplexDType = ComplexDType.float64
](arr: ComplexNDArray[dtype]) raises -> ComplexNDArray[dtype]:
    """Internal unnormalized inverse FFT helper function."""
    if arr.ndim != 1:
        raise Error("FFT currently only supports 1D arrays")

    var n: Int = arr.shape[0]
    if n <= 1:
        return ComplexNDArray[dtype](re=arr._re.copy(), im=arr._im.copy())

    if (n & (n - 1)) != 0:
        raise Error(
            "FFT currently only supports arrays with length that is a power"
            " of 2"
        )

    var half_size = n // 2
    var even_indices = ComplexNDArray[dtype](NDArrayShape(half_size))
    var odd_indices = ComplexNDArray[dtype](NDArrayShape(half_size))

    for i in range(half_size):
        even_indices[Item(i)] = arr[Item(2 * i)]
        odd_indices[Item(i)] = arr[Item(2 * i + 1)]

    var even_ifft = _ifft_unnormalized[dtype](even_indices)
    var odd_ifft = _ifft_unnormalized[dtype](odd_indices)

    var result = ComplexNDArray[dtype](arr.shape)

    for k in range(half_size):
        var angle = (
            2.0
            * Constants.pi
            * Scalar[dtype._dtype](k)
            / Scalar[dtype._dtype](n)
        )
        var twiddle = ComplexSIMD[dtype](
            cos(angle).cast[dtype._dtype](), sin(angle).cast[dtype._dtype]()
        )

        var twiddle_odd = twiddle * odd_ifft[Item(k)]

        result[Item(k)] = even_ifft[Item(k)] + twiddle_odd
        result[Item(k + half_size)] = even_ifft[Item(k)] - twiddle_odd

    return result^


fn ifft[
    dtype: ComplexDType = ComplexDType.float64
](arr: ComplexNDArray[dtype]) raises -> ComplexNDArray[dtype]:
    """Computes the Inverse Fast Fourier Transform using the Cooley-Tukey algorithm.

    Parameters:
        dtype: The data type of the complex elements (ComplexDType).

    Args:
        arr: Input complex array to transform. Must be 1-dimensional with length
             that is a power of 2.

    Returns:
        ComplexNDArray containing the Inverse Fast Fourier Transform of the input array.
        The output has the same shape and dtype as the input.

    Raises:
        Error: If the input array is not 1-dimensional.
        Error: If the array length is not a power of 2.

    Example:
        ```mojo
        from numojo.prelude import *
        from scijo.fft import ifft
        var input = nm.zeros[nm.cf32](nm.Shape(8))
        var transformed = ifft(input)
        ```

    Note:
        This implementation currently only supports 1D arrays with power-of-2 lengths.
    """
    var n: Int = arr.shape[0]
    var result = _ifft_unnormalized[dtype](arr)

    var inv_n = CScalar[dtype](1.0, 1.0) / CScalar[dtype](n, n)
    # result = result * inv_n # ! There's some mistakes in NuMojo mul overload.
    for i in range(n):
        result.store[width=1](i, result.load[width=1](i) * inv_n)
    # TODO: Fix problem with vectorize skipping some elements
    # alias simd_width = simd_width_of[dtype._dtype]()
    # @parameter
    # fn closure[width: Int](i: Int):
    #     try:
    #         result.store[width=width](i, result.load[width=width](i) * inv_n)
    #     except:
    #         pass
    # vectorize[closure, simd_width](n)
    return result^
