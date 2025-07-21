"""
SciJo - Scientific Computing Library for Mojo
==============================================

Differentiate Module - Utility Functions
-----------------------------------------

This module provides utility functions and data structures for numerical differentiation
using finite difference methods. It implements central, forward, and backward finite 
difference coefficient tables based on the formulas from Wikipedia's finite difference 
coefficient reference.

Author: Shivasankar K.A
Version: 0.1.0
Date: July 2025

References:
- Wikipedia: Finite difference coefficient 
  https://en.wikipedia.org/wiki/Finite_difference_coefficient
- Fornberg, B. (1988). Generation of Finite Difference Formulas on Arbitrarily
  Spaced Grids. Mathematics of Computation, 51(184), 699-706.
"""


struct Result[dtype: DType]():
    """Result structure for numerical differentiation operations.

    This structure encapsulates the results of derivative computations, including
    the computed derivative value, convergence information, and diagnostic data
    for numerical analysis.

    Type Parameters:
        dtype: The floating-point data type (DType.float32, DType.float64, etc.)

    Fields:
        success: Whether the computation converged successfully
        df: The computed derivative value
        error: Estimated error or convergence tolerance achieved
        nit: Number of iterations performed
        nfev: Number of function evaluations used
        x: The point at which the derivative was evaluated

    Usage:
        Creates a result object containing derivative computation results
        with convergence and diagnostic information.
    """

    var success: Bool
    var df: Scalar[dtype]
    var error: Scalar[dtype]
    var nit: Int
    var nfev: Int
    var x: Scalar[dtype]

    fn __init__(
        out self,
        success: Bool,
        df: Scalar[dtype],
        error: Scalar[dtype],
        nit: Int,
        nfev: Int,
        x: Scalar[dtype] = 0.0,
    ):
        self.success = success
        self.df = df
        self.error = error
        self.nit = nit
        self.nfev = nfev
        self.x = x

    fn __str__(self) raises -> String:
        return String(
            "DerivativeResult"
            + "\n"
            + "success={}"
            + "\n"
            + "df={}"
            + "\n"
            + "error={}"
            + "\n"
            + "nit={}"
            + "\n"
            + "nfev={}"
            + "\n"
            + "x={})"
        ).format(self.success, self.df, self.error, self.nit, self.nfev, self.x)

    fn write_to[W: Writer](self, mut writer: W):
        try:
            writer.write(
                String(
                    "\n"
                    + "DerivativeResult:"
                    + "\n"
                    + "\t success = {}"
                    + "\n"
                    + "\t df = {}"
                    + "\n"
                    + "\t error = {}"
                    + "\n"
                    + "\t nit = {}"
                    + "\n"
                    + "\t nfev = {}"
                    + "\n"
                    + "\t x = {}"
                    + "\n"
                ).format(
                    self.success,
                    self.df,
                    self.error,
                    self.nit,
                    self.nfev,
                    self.x,
                )
            )
        except e:
            writer.write(
                "Cannot convert DerivativeResult to string.\n" + String(e)
            )


@parameter
fn generate_central_finite_difference_table[
    dtype: DType
]() -> Dict[Int, List[Scalar[dtype]]]:
    """Generates central finite difference coefficients for first-order derivatives.

    This function creates a lookup table of finite difference coefficients for
    central difference approximations of first derivatives. The coefficients are
    based on the standard central difference formulas that provide symmetric
    approximations around the evaluation point.

    Central difference methods offer the best accuracy for smooth functions when
    function evaluations are available on both sides of the evaluation point.
    The method uses symmetric stencils: f'(x) ≈ Σ(c_i * f(x + i*h)) / h.

    Parameters:
        dtype: The floating-point data type for the coefficients
               (e.g., DType.float32, DType.float64).

    Returns:
        Dict[Int, List[Scalar[dtype]]] containing coefficient arrays indexed by
        accuracy order. Available orders: 2, 4, 6, 8 with corresponding
        truncation errors O(h²), O(h⁴), O(h⁶), O(h⁸).

    Example:
        ```mojo
        alias coeffs = generate_central_finite_difference_table[DType.float64]()
        var order2 = coeffs[2]  # [-0.5, 0.0, 0.5] for 3-point stencil
        ```

    Note:
        Higher-order methods require more function evaluations but provide
        exponentially better accuracy for smooth functions. Order 8 uses
        9 function evaluation points.
    """
    var coefficients = Dict[Int, List[Scalar[dtype]]]()

    # Order 2 central difference: [-1/2, 0, 1/2] at points [-1, 0, 1]
    coefficients[2] = List[Scalar[dtype]](
        Scalar[dtype](-0.5), Scalar[dtype](0.0), Scalar[dtype](0.5)
    )

    # Order 4 central difference: [1/12, -2/3, 0, 2/3, -1/12] at points [-2, -1, 0, 1, 2]
    coefficients[4] = List[Scalar[dtype]](
        Scalar[dtype](1.0 / 12.0),
        Scalar[dtype](-2.0 / 3.0),
        Scalar[dtype](0.0),
        Scalar[dtype](2.0 / 3.0),
        Scalar[dtype](-1.0 / 12.0),
    )

    # Order 6 central difference: [-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60] at points [-3, -2, -1, 0, 1, 2, 3]
    coefficients[6] = List[Scalar[dtype]](
        Scalar[dtype](-1.0 / 60.0),
        Scalar[dtype](3.0 / 20.0),
        Scalar[dtype](-3.0 / 4.0),
        Scalar[dtype](0.0),
        Scalar[dtype](3.0 / 4.0),
        Scalar[dtype](-3.0 / 20.0),
        Scalar[dtype](1.0 / 60.0),
    )

    # Order 8 central difference: [1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280] at points [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    coefficients[8] = List[Scalar[dtype]](
        Scalar[dtype](1.0 / 280.0),
        Scalar[dtype](-4.0 / 105.0),
        Scalar[dtype](1.0 / 5.0),
        Scalar[dtype](-4.0 / 5.0),
        Scalar[dtype](0.0),
        Scalar[dtype](4.0 / 5.0),
        Scalar[dtype](-1.0 / 5.0),
        Scalar[dtype](4.0 / 105.0),
        Scalar[dtype](-1.0 / 280.0),
    )

    return coefficients


@parameter
fn generate_forward_finite_difference_table[
    dtype: DType
]() -> Dict[Int, List[Scalar[dtype]]]:
    """Generates forward finite difference coefficients for first-order derivatives.

    This function creates a lookup table of finite difference coefficients for
    forward difference approximations of first derivatives. Forward differences
    use only function values at the evaluation point and points ahead in the
    positive direction.

    Forward difference methods are essential when function evaluations are only
    available in the forward direction (e.g., at domain boundaries, or when
    backward evaluations are not feasible). The method uses: f'(x) ≈ Σ(c_i * f(x + i*h)) / h.

    Parameters:
        dtype: The floating-point data type for the coefficients
               (e.g., DType.float32, DType.float64).

    Returns:
        Dict[Int, List[Scalar[dtype]]] containing coefficient arrays indexed by
        accuracy order. Available orders: 1, 2, 3, 4, 5, 6 with corresponding
        truncation errors O(h), O(h²), O(h³), O(h⁴), O(h⁵), O(h⁶).

    Example:
        ```mojo
        alias coeffs = generate_forward_finite_difference_table[DType.float64]()
        var order1 = coeffs[1]  # [-1, 1] for simple 2-point forward difference
        ```

    Note:
        Forward differences are generally less accurate than central differences
        for the same number of points, but are necessary at domain boundaries
        or when only forward evaluations are possible.
    """
    var coefficients = Dict[Int, List[Scalar[dtype]]]()

    # Order 1 forward difference: [-1, 1] at points [0, 1]
    coefficients[1] = List[Scalar[dtype]](
        Scalar[dtype](-1.0), Scalar[dtype](1.0)
    )

    # Order 2 forward difference: [-3/2, 2, -1/2] at points [0, 1, 2]
    coefficients[2] = List[Scalar[dtype]](
        Scalar[dtype](-3.0 / 2.0), Scalar[dtype](2.0), Scalar[dtype](-1.0 / 2.0)
    )

    # Order 3 forward difference: [-11/6, 3, -3/2, 1/3] at points [0, 1, 2, 3]
    coefficients[3] = List[Scalar[dtype]](
        Scalar[dtype](-11.0 / 6.0),
        Scalar[dtype](3.0),
        Scalar[dtype](-3.0 / 2.0),
        Scalar[dtype](1.0 / 3.0),
    )

    # Order 4 forward difference: [-25/12, 4, -3, 4/3, -1/4] at points [0, 1, 2, 3, 4]
    coefficients[4] = List[Scalar[dtype]](
        Scalar[dtype](-25.0 / 12.0),
        Scalar[dtype](4.0),
        Scalar[dtype](-3.0),
        Scalar[dtype](4.0 / 3.0),
        Scalar[dtype](-1.0 / 4.0),
    )

    # Order 5 forward difference: [-137/60, 5, -5, 10/3, -5/4, 1/5] at points [0, 1, 2, 3, 4, 5]
    coefficients[5] = List[Scalar[dtype]](
        Scalar[dtype](-137.0 / 60.0),
        Scalar[dtype](5.0),
        Scalar[dtype](-5.0),
        Scalar[dtype](10.0 / 3.0),
        Scalar[dtype](-5.0 / 4.0),
        Scalar[dtype](1.0 / 5.0),
    )

    # Order 6 forward difference: [-49/20, 6, -15/2, 20/3, -15/4, 6/5, -1/6] at points [0, 1, 2, 3, 4, 5, 6]
    coefficients[6] = List[Scalar[dtype]](
        Scalar[dtype](-49.0 / 20.0),
        Scalar[dtype](6.0),
        Scalar[dtype](-15.0 / 2.0),
        Scalar[dtype](20.0 / 3.0),
        Scalar[dtype](-15.0 / 4.0),
        Scalar[dtype](6.0 / 5.0),
        Scalar[dtype](-1.0 / 6.0),
    )

    return coefficients


@parameter
fn generate_backward_finite_difference_table[
    dtype: DType
]() -> Dict[Int, List[Scalar[dtype]]]:
    """Generates backward finite difference coefficients for first-order derivatives.

    This function creates a lookup table of finite difference coefficients for
    backward difference approximations of first derivatives. Backward differences
    use only function values at the evaluation point and points behind in the
    negative direction.

    Backward difference methods are essential when function evaluations are only
    available in the backward direction (e.g., at domain boundaries, or when
    forward evaluations are not feasible). The method uses: f'(x) ≈ Σ(c_i * f(x - i*h)) / h.

    Parameters:
        dtype: The floating-point data type for the coefficients
               (e.g., DType.float32, DType.float64).

    Returns:
        Dict[Int, List[Scalar[dtype]]] containing coefficient arrays indexed by
        accuracy order. Available orders: 1, 2, 3, 4, 5, 6 with corresponding
        truncation errors O(h), O(h²), O(h³), O(h⁴), O(h⁵), O(h⁶).

    Example:
        ```mojo
        alias coeffs = generate_backward_finite_difference_table[DType.float64]()
        var order2 = coeffs[2]  # [0.5, -2, 1.5] for 3-point backward difference
        ```

    Note:
        Backward differences are derived from forward differences by reversing
        the stencil and adjusting signs for odd derivatives. They provide the
        same accuracy as forward differences but in the opposite direction.
    """
    var coefficients = Dict[Int, List[Scalar[dtype]]]()

    # Order 1 backward difference: [-1, 1] at points [-1, 0] (reversed from forward)
    coefficients[1] = List[Scalar[dtype]](
        Scalar[dtype](-1.0), Scalar[dtype](1.0)
    )

    # Order 2 backward difference: [1/2, -2, 3/2] at points [-2, -1, 0] (reversed from forward with sign change for odd derivative)
    coefficients[2] = List[Scalar[dtype]](
        Scalar[dtype](1.0 / 2.0), Scalar[dtype](-2.0), Scalar[dtype](3.0 / 2.0)
    )

    # Order 3 backward difference: [-1/3, 3/2, -3, 11/6] at points [-3, -2, -1, 0] (reversed from forward with sign change)
    var coeff3 = List[Scalar[dtype]]()
    coeff3.append(Scalar[dtype](-1.0 / 3.0))
    coeff3.append(Scalar[dtype](3.0 / 2.0))
    coeff3.append(Scalar[dtype](-3.0))
    coeff3.append(Scalar[dtype](11.0 / 6.0))
    coefficients[3] = coeff3

    # Order 4 backward difference: [1/4, -4/3, 3, -4, 25/12] at points [-4, -3, -2, -1, 0] (reversed from forward with sign change)
    coefficients[4] = List[Scalar[dtype]](
        Scalar[dtype](1.0 / 4.0),
        Scalar[dtype](-4.0 / 3.0),
        Scalar[dtype](3.0),
        Scalar[dtype](-4.0),
        Scalar[dtype](25.0 / 12.0),
    )

    # Order 5 backward difference: [-1/5, 5/4, -10/3, 5, -5, 137/60] at points [-5, -4, -3, -2, -1, 0] (reversed from forward with sign change)
    coefficients[5] = List[Scalar[dtype]](
        Scalar[dtype](-1.0 / 5.0),
        Scalar[dtype](5.0 / 4.0),
        Scalar[dtype](-10.0 / 3.0),
        Scalar[dtype](5.0),
        Scalar[dtype](-5.0),
        Scalar[dtype](137.0 / 60.0),
    )

    # Order 6 backward difference: [1/6, -6/5, 15/4, -20/3, 15/2, -6, 49/20] at points [-6, -5, -4, -3, -2, -1, 0] (reversed from forward with sign change)
    coefficients[6] = List[Scalar[dtype]](
        Scalar[dtype](1.0 / 6.0),
        Scalar[dtype](-6.0 / 5.0),
        Scalar[dtype](15.0 / 4.0),
        Scalar[dtype](-20.0 / 3.0),
        Scalar[dtype](15.0 / 2.0),
        Scalar[dtype](-6.0),
        Scalar[dtype](49.0 / 20.0),
    )

    return coefficients
