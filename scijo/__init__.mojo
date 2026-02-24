# ===----------------------------------------------------------------------=== #
# Scijo: A Scientific Computing Library for Mojo
# Distributed under the Apache 2.0 License with LLVM Exceptions.
# See LICENSE and the LLVM License for more information.
# https://github.com/Mojo-Numerics-and-Algorithms-group/NuMojo/blob/main/LICENSE
# https://llvm.org/LICENSE.txt
#  ===----------------------------------------------------------------------=== #
"""
SciJo Top-Level Package (`numojo`)
==================================

Welcome to SciJo, a scientific computing library built for the Mojo programming language.

This top-level package exposes the core components of SciJo, including array types, error handling, and type definitions,
as well as a suite of modules for advanced numerical tasks.

Available Modules
-----------------
- `constants`: Common mathematical and physical constants used throughout the library.
- `differentiate`: Tools for numerical differentiation and gradient computation.
- `integrate`: Numerical integration routines for single and multi-dimensional problems.
- `fft`: Fast Fourier Transform algorithms for signal processing and spectral analysis.
- `interpolate`: Interpolation methods for estimating values between data points.
- `optimize`: Optimization algorithms for solving minimization and maximization problems.

Explore the documentation for each module to get started with scientific computing in Mojo using SciJo.
"""

from numojo.core.error import NumojoError
from numojo.core import NDArray, Item, NDArrayShape, NDArrayStrides, Shape

from numojo.core import (
    i8,
    i64,
    i128,
    i256,
    int,
    u8,
    u16,
    u32,
    u64,
    u128,
    u256,
    uint,
    bf16,
    f16,
    f32,
    f64,
    boolean,
)

from numojo.core.type_aliases import CScalar, ComplexScalar
from numojo.core import ComplexDType

from numojo.core import (
    ComplexSIMD,
    ci8,
    ci64,
    ci128,
    ci256,
    cint,
    cu8,
    cu16,
    cu32,
    cu64,
    cu128,
    cu256,
    cuint,
    cbf16,
    cf16,
    cf32,
    cf64,
    cboolean,
    cinvalid,
)
