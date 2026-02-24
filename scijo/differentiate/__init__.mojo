# ===----------------------------------------------------------------------=== #
# Scijo: Differentiate
# Distributed under the Apache 2.0 License with LLVM Exceptions.
# See LICENSE and the LLVM License for more information.
# https://github.com/Mojo-Numerics-and-Algorithms-group/NuMojo/blob/main/LICENSE
# https://llvm.org/LICENSE.txt
#  ===----------------------------------------------------------------------=== #
"""Differentiate Module (scijo.differentiate)

The `differentiate` module provides tools for numerical differentiation and gradient computation.
It includes functions for calculating derivatives, Jacobians, and other related operations essential for scientific computing tasks such as optimization, sensitivity analysis, and solving differential equations.
"""
from .derivative import derivative
from .jacobian import jacobian
