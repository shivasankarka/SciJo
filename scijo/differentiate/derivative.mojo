"""
SciJo - Scientific Computing Library for Mojo
==============================================

Differentiate Module - Numerical Differentiation
------------------------------------------------

This module implements numerical differentiation using finite difference methods.
It provides functions to compute first-order derivatives of scalar functions using
central, forward, and backward finite difference schemes with adaptive step sizing
and Richardson extrapolation for improved accuracy.

The implementation is designed to mimic SciPy's derivative function behavior while
providing additional flexibility in choosing the finite difference method and
controlling convergence parameters taking advantage of Mojo's metaprogramming techniques.

Author: Shivasankar K.A
Version: 0.1.0
Date: July 2025

Mathematical Background:
- Central differences: f'(x) ≈ Σ(c_i * f(x + i*h)) / h, symmetric around x
- Forward differences: f'(x) ≈ Σ(c_i * f(x + i*h)) / h, using x and forward points
- Backward differences: f'(x) ≈ Σ(c_i * f(x - i*h)) / h, using x and backward points

References:
- Fornberg, B. (1988). Generation of Finite Difference Formulas on Arbitrarily
  Spaced Grids. Mathematics of Computation, 51(184), 699-706.
- Scipy.derivative documentation
- Wikipedia: Finite difference coefficient
"""

from .numojo.prelude import *
from .utility import (
    Result,
    generate_central_finite_difference_table,
    generate_forward_finite_difference_table,
    generate_backward_finite_difference_table,
)

fn _derivative_central_difference[
    dtype: DType,
    func: fn (x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[
        dtype
    ],
](
    x0: Scalar[dtype],
    args: Optional[List[Scalar[dtype]]],
    first_order_coefficients: Dict[Int, List[Scalar[dtype]]],
    tolerance: Dict[String, Scalar[dtype]] = {"atol": 1e-6, "rtol": 1e-6},
    order: Int = 8,
    initial_step: Scalar[dtype] = 0.5,
    step_factor: Scalar[dtype] = 2.0,
    max_iter: Int = 10,
) raises -> Result[dtype]:
    """Computes first derivative using central finite difference method.

    This function implements the central finite difference method for computing
    first-order derivatives. It uses symmetric stencils around the evaluation
    point and employs Richardson extrapolation with adaptive step size reduction
    to achieve high accuracy.

    The central difference method provides the highest accuracy for smooth functions
    when function evaluations are available on both sides of the evaluation point.
    The algorithm iteratively refines the step size until convergence is achieved
    within the specified tolerance.

    Parameters:
        dtype: The floating-point data type (e.g., DType.float32, DType.float64).
        func: Function to differentiate with signature fn(x, args) -> Scalar[dtype].

    Args:
        x0: Point at which to evaluate the derivative.
        args: Optional arguments to pass to the function.
        first_order_coefficients: Precomputed finite difference coefficient table.
        tolerance: Convergence tolerances with "atol" and "rtol" keys.
        order: Accuracy order (2, 4, 6, or 8). Higher orders use more function evaluations.
        initial_step: Initial step size for finite differences.
        step_factor: Factor by which to reduce step size in each iteration.
        max_iter: Maximum number of Richardson extrapolation iterations.

    Returns:
        Result[dtype] containing the computed derivative, convergence information,
        and diagnostic data including number of iterations and function evaluations.

    Raises:
        Error: If the specified order is not supported (must be 2, 4, 6, or 8).

    Note:
        Central differences require function evaluations on both sides of x0.
        For boundary points or when this is not possible, use forward or
        backward differences instead.
    """
    var central_diff: Scalar[dtype] = 0.0
    var prev_diff: Scalar[dtype] = 0.0
    var atol: Scalar[dtype] = tolerance["atol"]
    var rtol: Scalar[dtype] = tolerance["rtol"]
    
    # Validate tolerance parameters
    if atol < 0:
        raise Error(
            "SciJo Derivative (Central): Invalid absolute tolerance.\n"
            "  Expected: atol ≥ 0\n"
            "  Got: atol = " + String(atol) + "\n"
            "  Note: Absolute tolerance must be non-negative for convergence testing."
        )
    if rtol < 0:
        raise Error(
            "SciJo Derivative (Central): Invalid relative tolerance.\n"
            "  Expected: rtol ≥ 0\n"
            "  Got: rtol = " + String(rtol) + "\n"
            "  Note: Relative tolerance must be non-negative for convergence testing."
        )
    
    # Validate step parameters
    if initial_step <= 0:
        raise Error(
            "SciJo Derivative (Central): Invalid initial step size.\n"
            "  Expected: initial_step > 0\n"
            "  Got: initial_step = " + String(initial_step) + "\n"
            "  Note: Step size must be positive for finite difference computation."
        )
    if step_factor <= 1:
        raise Error(
            "SciJo Derivative (Central): Invalid step reduction factor.\n"
            "  Expected: step_factor > 1\n"
            "  Got: step_factor = " + String(step_factor) + "\n"
            "  Note: Step factor must be > 1 for Richardson extrapolation convergence."
        )
    
    # Validate iteration parameters
    if max_iter <= 0:
        raise Error(
            "SciJo Derivative (Central): Invalid maximum iterations.\n"
            "  Expected: max_iter > 0\n"
            "  Got: max_iter = " + String(max_iter) + "\n"
            "  Note: At least one iteration is required for derivative computation."
        )
    
    var coefficients: List[Scalar[dtype]]
    if order in (2, 4, 6, 8):
        coefficients = first_order_coefficients[order]
    else:
        raise Error(
            "SciJo Derivative (Central): Invalid accuracy order specified.\n"
            "  Expected: order ∈ {2, 4, 6, 8}\n"
            "  Got: order = " + String(order) + "\n"
            "  Note: Higher orders provide better accuracy but require more function evaluations.\n"
            "  Available orders map to truncation errors: 2→O(h²), 4→O(h⁴), 6→O(h⁶), 8→O(h⁸)"
        )
    var step: Scalar[dtype] = initial_step

    for i in range(max_iter):
        central_diff = 0.0  # Reset for each iteration
        var j: Int = 0
        for ref coeff in coefficients:
            central_diff += coeff * func(x0 + step * (j - len(coefficients) // 2), args)
            j += 1
        central_diff /= step
        if i > 0:
            var diff_change = abs(central_diff - prev_diff)
            var tolerance_threshold = atol + rtol * abs(central_diff)
            if diff_change < tolerance_threshold:
                return Result[dtype](
                    success=True,
                    df=central_diff,
                    error=diff_change,
                    nit=i + 1,
                    nfev=(i + 1) * len(coefficients),
                )

        prev_diff = central_diff
        step /= step_factor

    return Result[dtype](
        success=False,
        df=central_diff,
        error=0.0,
        nit=max_iter,
        nfev=max_iter * len(coefficients),
    )


fn _derivative_forward_difference[
    dtype: DType,
    func: fn (x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[
        dtype
    ],
](
    x0: Scalar[dtype],
    args: Optional[List[Scalar[dtype]]],
    first_order_coefficients: Dict[Int, List[Scalar[dtype]]],
    tolerance: Dict[String, Scalar[dtype]] = {"atol": 1e-6, "rtol": 1e-6},
    order: Int = 8,
    initial_step: Scalar[dtype] = 0.5,
    step_factor: Scalar[dtype] = 2.0,
    max_iter: Int = 10,
) raises -> Result[dtype]:
    """Computes first derivative using forward finite difference method.

    This function implements the forward finite difference method for computing
    first-order derivatives. It uses function evaluations only at the evaluation
    point and points in the forward direction, making it suitable for boundary
    conditions or when backward evaluations are not available.

    The forward difference method is essential at the left boundary of domains
    or when the function is undefined for x < x0. While generally less accurate
    than central differences, higher-order forward differences can achieve
    comparable accuracy with more function evaluations.

    Parameters:
        dtype: The floating-point data type (e.g., DType.float32, DType.float64).
        func: Function to differentiate with signature fn(x, args) -> Scalar[dtype].

    Args:
        x0: Point at which to evaluate the derivative.
        args: Optional arguments to pass to the function.
        first_order_coefficients: Precomputed finite difference coefficient table.
        tolerance: Convergence tolerances with "atol" and "rtol" keys.
        order: Accuracy order (1, 2, 3, 4, 5, or 6). Higher orders use more function evaluations.
        initial_step: Initial step size for finite differences.
        step_factor: Factor by which to reduce step size in each iteration.
        max_iter: Maximum number of Richardson extrapolation iterations.

    Returns:
        Result[dtype] containing the computed derivative, convergence information,
        and diagnostic data including number of iterations and function evaluations.

    Raises:
        Error: If the specified order is not supported (must be 1, 2, 3, 4, 5, or 6).

    Note:
        Forward differences only require function evaluations at x0 and x0+h, x0+2h, etc.
        This makes them ideal for boundary value problems or when the function
        domain has constraints.
    """
    var central_diff: Scalar[dtype] = 0.0
    var prev_diff: Scalar[dtype] = 0.0
    var atol: Scalar[dtype] = tolerance["atol"]
    var rtol: Scalar[dtype] = tolerance["rtol"]
    
    # Validate tolerance parameters
    if atol < 0:
        raise Error(
            "SciJo Derivative (Forward): Invalid absolute tolerance.\n"
            "  Expected: atol ≥ 0\n"
            "  Got: atol = " + String(atol) + "\n"
            "  Note: Absolute tolerance must be non-negative for convergence testing."
        )
    if rtol < 0:
        raise Error(
            "SciJo Derivative (Forward): Invalid relative tolerance.\n"
            "  Expected: rtol ≥ 0\n"
            "  Got: rtol = " + String(rtol) + "\n"
            "  Note: Relative tolerance must be non-negative for convergence testing."
        )
    
    # Validate step parameters
    if initial_step <= 0:
        raise Error(
            "SciJo Derivative (Forward): Invalid initial step size.\n"
            "  Expected: initial_step > 0\n"
            "  Got: initial_step = " + String(initial_step) + "\n"
            "  Note: Step size must be positive for finite difference computation."
        )
    if step_factor <= 1:
        raise Error(
            "SciJo Derivative (Forward): Invalid step reduction factor.\n"
            "  Expected: step_factor > 1\n"
            "  Got: step_factor = " + String(step_factor) + "\n"
            "  Note: Step factor must be > 1 for Richardson extrapolation convergence."
        )
    
    # Validate iteration parameters
    if max_iter <= 0:
        raise Error(
            "SciJo Derivative (Forward): Invalid maximum iterations.\n"
            "  Expected: max_iter > 0\n"
            "  Got: max_iter = " + String(max_iter) + "\n"
            "  Note: At least one iteration is required for derivative computation."
        )
    
    if order in (1, 2, 3, 4, 5, 6):
        coefficients = first_order_coefficients[order]
    else:
        raise Error(
            "SciJo Derivative (Forward): Invalid accuracy order specified.\n"
            "  Expected: order ∈ {1, 2, 3, 4, 5, 6}\n"
            "  Got: order = " + String(order) + "\n"
            "  Note: Forward differences are limited to order 6 due to numerical stability.\n"
            "  Available orders map to truncation errors: 1→O(h), 2→O(h²), ..., 6→O(h⁶)"
        )
    var step: Scalar[dtype] = initial_step

    for i in range(max_iter):
        central_diff = 0.0  # Reset for each iteration
        var j: Int = 0
        for ref coeff in coefficients:
            central_diff += coeff * func(x0 + step * (j), args)
            j += 1
        central_diff /= step
        if i > 0:
            var diff_change = abs(central_diff - prev_diff)
            var tolerance_threshold = atol + rtol * abs(central_diff)
            if diff_change < tolerance_threshold:
                return Result[dtype](
                    success=True,
                    df=central_diff,
                    error=diff_change,
                    nit=i + 1,
                    nfev=(i + 1) * len(coefficients),
                )

        prev_diff = central_diff
        step /= step_factor

    return Result[dtype](
        success=False,
        df=central_diff,
        error=0.0,
        nit=max_iter,
        nfev=max_iter * len(coefficients),
    )


fn _derivative_backward_difference[
    dtype: DType,
    func: fn (x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[
        dtype
    ],
](
    x0: Scalar[dtype],
    args: Optional[List[Scalar[dtype]]],
    first_order_coefficients: Dict[Int, List[Scalar[dtype]]],
    tolerance: Dict[String, Scalar[dtype]] = {"atol": 1e-6, "rtol": 1e-6},
    order: Int = 8,
    initial_step: Scalar[dtype] = 0.5,
    step_factor: Scalar[dtype] = 2.0,
    max_iter: Int = 10,
) raises -> Result[dtype]:
    """Computes first derivative using backward finite difference method.

    This function implements the backward finite difference method for computing
    first-order derivatives. It uses function evaluations only at the evaluation
    point and points in the backward direction, making it suitable for boundary
    conditions or when forward evaluations are not available.

    The backward difference method is essential at the right boundary of domains
    or when the function is undefined for x > x0. The coefficients are derived
    from forward differences by reversing the stencil and adjusting signs
    appropriately for first derivatives.

    Parameters:
        dtype: The floating-point data type (e.g., DType.float32, DType.float64).
        func: Function to differentiate with signature fn(x, args) -> Scalar[dtype].

    Args:
        x0: Point at which to evaluate the derivative.
        args: Optional arguments to pass to the function.
        first_order_coefficients: Precomputed finite difference coefficient table.
        tolerance: Convergence tolerances with "atol" and "rtol" keys.
        order: Accuracy order (1, 2, 3, 4, 5, or 6). Higher orders use more function evaluations.
        initial_step: Initial step size for finite differences.
        step_factor: Factor by which to reduce step size in each iteration.
        max_iter: Maximum number of Richardson extrapolation iterations.

    Returns:
        Result[dtype] containing the computed derivative, convergence information,
        and diagnostic data including number of iterations and function evaluations.

    Raises:
        Error: If the specified order is not supported (must be 1, 2, 3, 4, 5, or 6).

    Note:
        Backward differences only require function evaluations at x0 and x0-h, x0-2h, etc.
        This makes them ideal for right boundary conditions or when the function
        domain has forward constraints.
    """
    var central_diff: Scalar[dtype] = 0.0
    var prev_diff: Scalar[dtype] = 0.0
    var atol: Scalar[dtype] = tolerance["atol"]
    var rtol: Scalar[dtype] = tolerance["rtol"]
    
    # Validate tolerance parameters
    if atol < 0:
        raise Error(
            "SciJo Derivative (Backward): Invalid absolute tolerance.\n"
            "  Expected: atol ≥ 0\n"
            "  Got: atol = " + String(atol) + "\n"
            "  Note: Absolute tolerance must be non-negative for convergence testing."
        )
    if rtol < 0:
        raise Error(
            "SciJo Derivative (Backward): Invalid relative tolerance.\n"
            "  Expected: rtol ≥ 0\n"
            "  Got: rtol = " + String(rtol) + "\n"
            "  Note: Relative tolerance must be non-negative for convergence testing."
        )
    
    # Validate step parameters
    if initial_step <= 0:
        raise Error(
            "SciJo Derivative (Backward): Invalid initial step size.\n"
            "  Expected: initial_step > 0\n"
            "  Got: initial_step = " + String(initial_step) + "\n"
            "  Note: Step size must be positive for finite difference computation."
        )
    if step_factor <= 1:
        raise Error(
            "SciJo Derivative (Backward): Invalid step reduction factor.\n"
            "  Expected: step_factor > 1\n"
            "  Got: step_factor = " + String(step_factor) + "\n"
            "  Note: Step factor must be > 1 for Richardson extrapolation convergence."
        )
    
    # Validate iteration parameters
    if max_iter <= 0:
        raise Error(
            "SciJo Derivative (Backward): Invalid maximum iterations.\n"
            "  Expected: max_iter > 0\n"
            "  Got: max_iter = " + String(max_iter) + "\n"
            "  Note: At least one iteration is required for derivative computation."
        )
    
    if order in (1, 2, 3, 4, 5, 6):
        coefficients = first_order_coefficients[order]
    else:
        raise Error(
            "SciJo Derivative (Backward): Invalid accuracy order specified.\n"
            "  Expected: order ∈ {1, 2, 3, 4, 5, 6}\n"
            "  Got: order = " + String(order) + "\n"
            "  Note: Backward differences are limited to order 6 due to numerical stability.\n"
            "  Available orders map to truncation errors: 1→O(h), 2→O(h²), ..., 6→O(h⁶)"
        )
    var step: Scalar[dtype] = initial_step

    for i in range(max_iter):
        central_diff = 0.0  # Reset for each iteration
        var j: Int = 0
        for ref coeff in coefficients:
            central_diff += coeff * func(x0 + step * (-j), args)
            j += 1
        central_diff /= step
        if i > 0:
            var diff_change = abs(central_diff - prev_diff)
            var tolerance_threshold = atol + rtol * abs(central_diff)
            if diff_change < tolerance_threshold:
                return Result[dtype](
                    success=True,
                    df=central_diff,
                    error=diff_change,
                    nit=i + 1,
                    nfev=(i + 1) * len(coefficients),
                )

        prev_diff = central_diff
        step /= step_factor

    return Result[dtype](
        success=False,
        df=central_diff,
        error=0.0,
        nit=max_iter,
        nfev=max_iter * len(coefficients),
    )


# first derivative only -> Central finite difference method
fn derivative[
    dtype: DType,
    func: fn (x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[
        dtype
    ],
    *,
    step_direction: Int = 0,
](
    x0: Scalar[dtype],
    args: Optional[List[Scalar[dtype]]],
    tolerance: Dict[String, Scalar[dtype]] = {"atol": 1e-6, "rtol": 1e-6},
    order: Int = 8,
    initial_step: Scalar[dtype] = 0.5,
    step_factor: Scalar[dtype] = 2.0,
    max_iter: Int = 10,
) raises -> Result[dtype]:
    """Computes the first derivative of a scalar function using finite differences.

    This function provides a unified interface for computing first-order derivatives
    using different finite difference methods. It automatically selects the appropriate
    method based on the step_direction parameter and uses Richardson extrapolation
    with adaptive step sizing for improved accuracy.

    The implementation is designed to be compatible with SciPy's derivative function
    while providing additional control over the finite difference method and
    convergence parameters. It supports central, forward, and backward differences
    with various accuracy orders.

    Parameters:
        dtype: The floating-point data type (e.g., DType.float32, DType.float64).
        func: Function to differentiate with signature fn(x, args) -> Scalar[dtype].
        step_direction: Direction of finite difference stencil (keyword-only parameter).

    Args:
        x0: Point at which to evaluate the derivative.
        args: Optional arguments to pass to the function.
        tolerance: Convergence tolerances with "atol" (absolute) and "rtol" (relative) keys.
        order: Accuracy order for finite differences. Valid ranges depend on method.
               Central differences support orders: 2, 4, 6, 8.
               Forward/Backward differences support orders: 1, 2, 3, 4, 5, 6.
        initial_step: Initial step size for finite differences.
        step_factor: Factor by which to reduce step size in each iteration (typically 2.0).
        max_iter: Maximum number of Richardson extrapolation iterations.

    Returns:
        Result[dtype] containing:
        - success: Whether the computation converged
        - df: The computed derivative value
        - error: Estimated error or final convergence criterion
        - nit: Number of iterations performed
        - nfev: Total number of function evaluations

    Raises:
        Error: If dtype is not a floating-point type.
        Error: If step_direction is not 0, 1, or -1.
        Error: If the specified order is not supported for the chosen method.

    Example:
        ```mojo
        fn square(x: Scalar[DType.float64], args: Optional[List[Scalar[DType.float64]]]) -> Scalar[DType.float64]:
            return x * x  # f(x) = x²
        
        # Compute derivative f'(2) = 2*2 = 4
        var result = derivative[DType.float64, square](
            x0=2.0,
            args=None,
            order=4,
            step_direction=0  # Use central differences
        )
        # result.df should be approximately 4.0
        ```

    Note:
        - step_direction=0: Central differences (default, highest accuracy)
        - step_direction=1: Forward differences (for left boundaries)
        - step_direction=-1: Backward differences (for right boundaries)
        
        Central differences generally provide the best accuracy but require
        function evaluations on both sides of x0. Use forward/backward
        differences at domain boundaries or when directional constraints exist.
    """
    constrained[
        dtype.is_floating_point(),
        msg="SciJo Derivative: Data type must be floating-point for numerical differentiation.\n"
        "  Floating-point types are required to handle fractional step sizes and avoid\n"
        "  precision loss during finite difference computations.\n"
        "  Supported types: DType.float16, DType.float32, DType.float64, DType.float128",
    ]()
    alias first_order_coefficients = generate_central_finite_difference_table[
        dtype
    ]()

    @parameter
    if step_direction == 0:
        alias first_order_coefficients = generate_central_finite_difference_table[
            dtype
        ]()
        return _derivative_central_difference[dtype, func](
            x0,
            args,
            first_order_coefficients,
            tolerance,
            order,
            initial_step,
            step_factor,
            max_iter,
        )
    elif step_direction == 1:
        alias first_order_coefficients = generate_forward_finite_difference_table[
            dtype
        ]()
        return _derivative_forward_difference[dtype, func](
            x0,
            args,
            first_order_coefficients,
            tolerance,
            order,
            initial_step,
            step_factor,
            max_iter,
        )
    elif step_direction == -1:
        alias first_order_coefficients = generate_backward_finite_difference_table[
            dtype
        ]()
        return _derivative_backward_difference[dtype, func](
            x0,
            args,
            first_order_coefficients,
            tolerance,
            order,
            initial_step,
            step_factor,
            max_iter,
        )
    else:
        raise Error(
            "SciJo Derivative: Invalid step direction parameter.\n"
            "  Expected: step_direction ∈ {-1, 0, 1}\n"
            "  Got: step_direction = " + String(step_direction) + "\n"
            "  Valid options:\n"
            "    • step_direction = 0:  Central differences (highest accuracy, requires f(x±h))\n"
            "    • step_direction = 1:  Forward differences (for left boundaries, uses f(x+h))\n"
            "    • step_direction = -1: Backward differences (for right boundaries, uses f(x-h))"
        )
