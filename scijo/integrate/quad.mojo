"""
Adaptive Quadrature Integration using Gauss-Kronrod Rules (QAGS Algorithm)

This module implements the QUADPACK QAGS algorithm used by SciPy's quad function.
The algorithm uses adaptive subdivision with 21-point Gauss-Kronrod quadrature rules
to accurately compute definite integrals.

Key features:
- 21-point Gauss-Kronrod quadrature for high accuracy
- Adaptive subdivision for handling difficult integrands
- Error estimation using difference between Gauss and Kronrod results
- Similar interface and behavior to SciPy's integrate.quad
- High-precision G10K21 coefficients from Advanpix (34 decimal places)

References:
- Piessens, R., de Doncker-Kapenga, E., Überhuber, C. W., & Kahaner, D. K. (1983).
  QUADPACK: A subroutine package for automatic integration. Springer-Verlag.
- SciPy documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
- Advanpix G10K21 coefficients: https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/
"""

from math import sqrt

from scijo.integrate.utility import (
    IntegralResult,
    kronrod_nodes,
    kronrod_weights,
    gauss_weights,
)

# QUAD function with adaptive Gauss-Kronrod integration        
fn quad[
    dtype: DType,
](
    func: fn (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]] = None
    ) -> Scalar[dtype],
    a: Scalar[dtype],
    b: Scalar[dtype],
    args: Optional[List[Scalar[dtype]]],
    epsabs: Scalar[dtype] = 1.49e-8,
    epsrel: Scalar[dtype] = 1.49e-8,
    limit: Int = 50,
) -> IntegralResult[dtype]:
    """`
    Adaptive quadrature using Gauss-Kronrod rules (QAGS algorithm).

    This is a Mojo implementation of the QUADPACK QAGS algorithm used by SciPy.
    It uses adaptive subdivision with 21-point Gauss-Kronrod quadrature rules.

    The algorithm works by:
    1. Evaluating the integral using both 10-point Gauss and 21-point Kronrod rules
    2. Using the difference as an error estimate
    3. If the error is too large, subdividing the interval and repeating
    4. Continuing until convergence or maximum subdivisions reached

    Args:
        func: The integrand function to integrate.
        a: Lower integration limit.
        b: Upper integration limit.
        args: Additional arguments for integrand.
        epsabs: Absolute error tolerance.
        epsrel: Relative error tolerance.
        limit: Maximum number of subintervals.

    Returns:
        Tuple of (integral_value, error_estimate, status_message).
    """
    if a == b:
        return IntegralResult(
            integral=Scalar[dtype](0),
            abserr=Scalar[dtype](0),
            neval=0,
            ier=0,
            last=0,
        )

    var sign = Scalar[dtype](1)
    var lower = a
    var upper = b
    if a > b:
        sign = Scalar[dtype](-1)
        lower = b
        upper = a

    var initial_result = _gauss_kronrod_21[dtype](func, lower, upper, args)
    var result_k = initial_result[0]
    var result_g = initial_result[1]  # Gauss result
    var neval = initial_result[2]

    var error_est = abs(result_k - result_g)

    var tolerance = max(epsabs, epsrel * abs(result_k))
    if error_est <= tolerance:
        return IntegralResult(
            integral=sign * result_k,
            abserr=error_est,
            neval=neval,
            ier=0,  # Success (QUADPACK convention: ier=0 means success)
            last=0,
        )

    var current_result = result_k
    var current_error = error_est
    var subdivisions = 0

    while current_error > tolerance and subdivisions < limit:
        var mid = (lower + upper) / 2

        var left_result = _gauss_kronrod_21[dtype](func, lower, mid, args)
        var left_k = left_result[0]
        var left_g = left_result[1]
        var left_error = abs(left_k - left_g)
        neval += left_result[2]

        var right_result = _gauss_kronrod_21[dtype](func, mid, upper, args)
        var right_k = right_result[0]
        var right_g = right_result[1]
        var right_error = abs(right_k - right_g)
        neval += right_result[2]

        current_result = left_k + right_k
        current_error = left_error + right_error
        tolerance = max(epsabs, epsrel * abs(current_result))

        if left_error > right_error:
            upper = mid
        else:
            lower = mid

        subdivisions += 1

    if current_error <= tolerance:
        return IntegralResult(
            integral=sign * current_result,
            abserr=current_error,
            neval=neval,
            ier=0,  # Success (QUADPACK convention: ier=0 means success)
            last=subdivisions,
        )
    else:  
        try:
            print(String("Maximum subdivisions {} reached in QUAD, returning current result").format(limit))
        except e:
            print("Error printing message:", e)
        return IntegralResult(
            integral=sign * current_result,
            abserr=current_error,
            neval=neval,
            ier=1,  # Maximum subdivisions reached (failure)
            last=subdivisions,
        )

# same quad function with compile time function
fn quad[
    dtype: DType,
    integrand: fn (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
](
    a: Scalar[dtype],
    b: Scalar[dtype],
    args: Optional[List[Scalar[dtype]]],
    epsabs: Scalar[dtype] = 1.49e-8,
    epsrel: Scalar[dtype] = 1.49e-8,
    limit: Int = 50,
) -> IntegralResult[dtype]:
    """
    Adaptive quadrature using Gauss-Kronrod rules (QAGS algorithm).

    This is a Mojo implementation of the QUADPACK QAGS algorithm used by SciPy.
    It uses adaptive subdivision with 21-point Gauss-Kronrod quadrature rules.

    The algorithm works by:
    1. Evaluating the integral using both 10-point Gauss and 21-point Kronrod rules
    2. Using the difference as an error estimate
    3. If the error is too large, subdividing the interval and repeating
    4. Continuing until convergence or maximum subdivisions reached

    Args:
        a: Lower integration limit.
        b: Upper integration limit.
        args: Additional arguments for integrand.
        epsabs: Absolute error tolerance.
        epsrel: Relative error tolerance.
        limit: Maximum number of subintervals.
    Returns:
        Tuple of (integral_value, error_estimate, status_message).
    """

    if a == b:
        return IntegralResult(
            integral=Scalar[dtype](0),
            abserr=Scalar[dtype](0),
            neval=0,
            ier=0,  
            last=0,
        )

    var sign = Scalar[dtype](1)
    var lower = a
    var upper = b
    if a > b:
        sign = Scalar[dtype](-1)
        lower = b
        upper = a

    var initial_result = _gauss_kronrod_21[dtype, integrand](lower, upper, args)
    var result_k = initial_result[0]
    var result_g = initial_result[1]  
    var neval = initial_result[2]

    var error_est = abs(result_k - result_g)

    var tolerance = max(epsabs, epsrel * abs(result_k))
    if error_est <= tolerance:
        return IntegralResult(
            integral=sign * result_k,
            abserr=error_est,
            neval=neval,
            ier=0,  # Success (QUADPACK convention: ier=0 means success)
            last=0,
        )

    var current_result = result_k
    var current_error = error_est
    var subdivisions = 0

    # look up references for better adaptive algorithms. 
    while current_error > tolerance and subdivisions < limit:
        var mid = (lower + upper) / 2

        var left_result = _gauss_kronrod_21[dtype, integrand](lower, mid, args)
        var left_k = left_result[0]
        var left_g = left_result[1]
        var left_error = abs(left_k - left_g)
        neval += left_result[2]

        var right_result = _gauss_kronrod_21[dtype, integrand](mid, upper, args)
        var right_k = right_result[0]
        var right_g = right_result[1]
        var right_error = abs(right_k - right_g)
        neval += right_result[2]

        current_result = left_k + right_k
        current_error = left_error + right_error

        tolerance = max(epsabs, epsrel * abs(current_result))

        if left_error > right_error:
            upper = mid
        else:
            lower = mid

        subdivisions += 1

    if current_error <= tolerance:
        return IntegralResult(
            integral=sign * current_result,
            abserr=current_error,
            neval=neval,
            ier=0,  # Success (QUADPACK convention: ier=0 means success)
            last=subdivisions,
        )
    else:  # subdivisions >= limit
        try:
            print(String("Maximum subdivisions {} reached in QUAD, returning current result").format(limit))
        except e:
            print("Error printing message:", e)
        return IntegralResult(
            integral=sign * current_result,
            abserr=current_error,
            neval=neval,
            ier=1,  # Maximum subdivisions reached (failure)
            last=subdivisions,
        )

fn _gauss_kronrod_21[
    dtype: DType,
    func: fn (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
](
    a: Scalar[dtype], b: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
) -> Tuple[Scalar[dtype], Scalar[dtype], Int]:
    """Compute integral using 21-point Gauss-Kronrod rule."""
    var center = (a + b) / 2
    var half_length = (b - a) / 2

    # Center point evaluation (node 0 = 0.0)
    var fc = func(center, args)
    var result_kronrod = Scalar[dtype](kronrod_weights[0]) * fc
    var result_gauss = Scalar[dtype](gauss_weights[0]) * fc
    var neval = 1

    # Process all pairs of symmetric points
    # For G10K21: Gauss points are at indices 1,3,5,7,9 in the Kronrod sequence
    var gauss_idx = 1

    for i in range(1, 11):
        # Transform node from [-1,1] to [a,b]: x_scaled = x * (b-a)/2 + (a+b)/2
        var node = Scalar[dtype](kronrod_nodes[i])
        var x1 = (
            center - half_length * node
        )  # Left point: (a+b)/2 - (b-a)/2 * node
        var x2 = (
            center + half_length * node
        )  # Right point: (a+b)/2 + (b-a)/2 * node

        var f1 = func(x1, args)
        var f2 = func(x2, args)
        var fsum = f1 + f2

        # Add to Kronrod result (all points contribute)
        result_kronrod += Scalar[dtype](kronrod_weights[i]) * fsum

        # Add to Gauss result only for Gauss points (odd indices: 1,3,5,7,9)
        if i % 2 == 1 and gauss_idx < 6:  # Guard against index overflow
            result_gauss += Scalar[dtype](gauss_weights[gauss_idx]) * fsum
            gauss_idx += 1

        neval += 2

    # Scale by interval length: w_scaled = w * (b-a)/2
    result_kronrod *= half_length
    result_gauss *= half_length

    return (result_kronrod, result_gauss, neval)

fn _gauss_kronrod_21[
    dtype: DType,
](
    func: fn (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    a: Scalar[dtype], b: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
) -> Tuple[Scalar[dtype], Scalar[dtype], Int]:
    """Compute integral using 21-point Gauss-Kronrod rule."""
    var center = (a + b) / 2
    var half_length = (b - a) / 2

    # Center point evaluation (node 0 = 0.0)
    var fc = func(center, args)
    var result_kronrod = Scalar[dtype](kronrod_weights[0]) * fc
    var result_gauss = Scalar[dtype](gauss_weights[0]) * fc
    var neval = 1

    # Process all pairs of symmetric points
    # For G10K21: Gauss points are at indices 1,3,5,7,9 in the Kronrod sequence
    var gauss_idx = 1

    for i in range(1, 11):
        # Transform node from [-1,1] to [a,b]: x_scaled = x * (b-a)/2 + (a+b)/2
        var node = Scalar[dtype](kronrod_nodes[i])
        var x1 = (
            center - half_length * node
        )  # Left point: (a+b)/2 - (b-a)/2 * node
        var x2 = (
            center + half_length * node
        )  # Right point: (a+b)/2 + (b-a)/2 * node

        var f1 = func(x1, args)
        var f2 = func(x2, args)
        var fsum = f1 + f2

        # Add to Kronrod result (all points contribute)
        result_kronrod += Scalar[dtype](kronrod_weights[i]) * fsum

        # Add to Gauss result only for Gauss points (odd indices: 1,3,5,7,9)
        if i % 2 == 1 and gauss_idx < 6:  # Guard against index overflow
            result_gauss += Scalar[dtype](gauss_weights[gauss_idx]) * fsum
            gauss_idx += 1

        neval += 2

    # Scale by interval length: w_scaled = w * (b-a)/2
    result_kronrod *= half_length
    result_gauss *= half_length

    return (result_kronrod, result_gauss, neval)
