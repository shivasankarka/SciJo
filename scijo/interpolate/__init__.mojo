# ===----------------------------------------------------------------------=== #
# Scijo: Interpolate
# Distributed under the Apache 2.0 License with LLVM Exceptions.
# See LICENSE and the LLVM License for more information.
# https://github.com/Mojo-Numerics-and-Algorithms-group/NuMojo/blob/main/LICENSE
# https://llvm.org/LICENSE.txt
#  ===----------------------------------------------------------------------=== #
"""Interpolate Module (scijo.interpolate)

The `interpolate` module provides tools for interpolating data. It includes
linear interpolation functions and callable interpolator objects for both
single-point and array-based evaluation.
"""

from .interpolate import interp1d, LinearInterpolator
