# ===----------------------------------------------------------------------=== #
# SciJo: FFT module for Mojo
# Distributed under the Apache 2.0 License.
# ===----------------------------------------------------------------------=== #
"""FFT Module (`scijo.fft`)
===========================
Provides Fast Fourier Transform operations for complex and real arrays.
It includes forward and inverse FFT using the Cooley-Tukey algorithm.

Available Functions
-------------------
- `fft`   — Compute the forward FFT (complex input).
- `ifft`  — Compute the inverse FFT (complex input).
- `rfft`  — Compute the FFT of a real array, returning N//2+1 bins.
- `irfft` — Compute the inverse FFT returning a real array.

Examples
--------
    ```mojo
    from scijo.fft import fft, ifft, rfft, irfft

    var arr = nm.linspace[cf32](CScalar[cf32](0, 0), CScalar[cf32](10, 10), num=10)
    var fft_arr = fft(arr)
    var time = ifft(fft_arr)
    ```
"""

from .fastfourier import fft, ifft, rfft, irfft
