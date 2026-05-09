# ===----------------------------------------------------------------------=== #
# SciJo: Optimize module for Mojo
# Distributed under the Apache 2.0 License.
# ===----------------------------------------------------------------------=== #
"""Scalar Root-Finding (`scijo.optimize.root_scalar`)
=====================================================
Scalar root-finding methods for nonlinear equations. Includes bracketing methods
(bisection), derivative-based methods (Newton-Raphson), and derivative-free
methods (secant).

Examples
--------
    ```mojo
    from scijo.optimize import root_scalar

    def f[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
        return x * x - 2

    var result = root_scalar[f64, f](bracket=(1.0, 2.0), method="bisect")
    print(result.root)       # ≈ 1.41421356
    print(result.converged)  # True
    ```
"""

from scijo.optimize.utility import RootResult

# TODO: check if we are using the right tolerance conditions in all methods.

# ===----------------------------------------------------------------------=== #
# Root scalar
# ===----------------------------------------------------------------------=== #


def root_scalar[
    dtype: DType,
    f: def[dtype: DType](
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    fprime: Optional[
        def[
            dtype: DType
        ](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[
            dtype
        ]
    ] = None,
    *,
    method: String = "bisect",
](
    args: Optional[List[Scalar[dtype]]] = None,
    x0: Optional[Scalar[dtype]] = None,
    x1: Optional[Scalar[dtype]] = None,
    bracket: Optional[Tuple[Scalar[dtype], Scalar[dtype]]] = None,
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> RootResult[dtype]:
    """Finds a root of a scalar function using the specified method.

    Parameters:
        dtype: The floating-point data type.
        f: Function f(x, args) -> Scalar[dtype] for which to find a root.
        fprime: Derivative function f'(x, args). Required for Newton's method.
        method: Root-finding algorithm: "bisect", "newton", or "secant".

    Args:
        args: Optional arguments forwarded to f and fprime.
        x0: Initial guess. Required for Newton and secant methods.
        x1: Second initial guess. Required for secant method.
        bracket: (a, b) tuple where f(a) and f(b) have opposite signs. Required for bisection.
        xtol: Absolute tolerance for convergence.
        rtol: Relative tolerance for convergence.
        maxiter: Maximum number of iterations.

    Raises:
        Error: If required inputs for the chosen method are missing or invalid.

    Returns:
        RootResult[dtype] containing the root, convergence status, iteration
        count, function evaluation count, and method name.
    """

    comptime if method == "newton":
        if not fprime:
            raise Error(
                "Scijo [root_scalar]: Derivative fprime must be provided for"
                " Newton's method."
            )
        return newton[dtype, f, fprime.value()](args, x0, xtol, rtol, maxiter)
    elif method == "bisect":
        if not bracket:
            raise Error(
                "Scijo [root_scalar]: Bracket must be provided for bisection"
                " method."
            )
        return bisect[dtype, f](args, bracket.value(), xtol, rtol, maxiter)
    elif method == "secant":
        if not (x0 and x1):
            raise Error(
                "Scijo [root_scalar]: Initial guesses x0 and x1 must be"
                " provided for secant method."
            )
        return secant[dtype, f](
            args, x0.value(), x1.value(), xtol, rtol, maxiter
        )
    else:
        raise Error(
            "Scijo [root_scalar]: Unsupported method: " + String(method)
        )


# ===----------------------------------------------------------------------=== #
# Root scalar methods
# ===----------------------------------------------------------------------=== #


def newton[
    dtype: DType,
    f: def[dtype: DType](
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    fprime: def[dtype: DType](
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    method: String = "newton",
](
    args: Optional[List[Scalar[dtype]]],
    x0: Optional[Scalar[dtype]] = None,
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> RootResult[dtype]:
    """Finds a root using the Newton-Raphson method.

    Terminates when the step size or function value falls below
    max(xtol, rtol * |x|), or when maxiter is reached.

    Parameters:
        dtype: The floating-point data type.
        f: Function for which to find a root.
        fprime: Derivative of f.
        method: Method name (unused, for signature compatibility).

    Args:
        args: Optional arguments forwarded to f and fprime.
        x0: Initial guess. Required.
        xtol: Absolute tolerance for convergence.
        rtol: Relative tolerance for convergence.
        maxiter: Maximum number of iterations.

    Raises:
        Error: If x0 is not provided or the derivative is zero at any step.

    Returns:
        RootResult[dtype] containing the root and convergence information.
    """
    var xn: Scalar[dtype]
    if x0:
        xn = x0.value()
    else:
        raise Error(
            "Scijo [newton]: Initial guess x0 must be provided for Newton's"
            " method."
        )

    var nit: Int = 0
    var nfev: Int = 0

    for _ in range(maxiter):
        var fx = f(xn, args)
        var fpx = fprime(xn, args)
        nfev += 2
        nit += 1

        if fpx == 0:
            raise Error(
                "Scijo [newton]: Derivative is zero. Newton-Raphson step would"
                " divide by zero."
            )

        var delta = fx / fpx
        var xn_next = xn - delta

        var tol_x = max(xtol, rtol * abs(xn_next))
        var tol_f = max(xtol, rtol * abs(fx))

        if abs(fx) <= tol_f or abs(delta) <= tol_x:
            return RootResult[dtype](
                root=xn_next,
                nit=nit,
                nfev=nfev,
                success=True,
                message="converged",
                method="newton",
            )

        xn = xn_next

    return RootResult[dtype](
        root=xn,
        nit=nit,
        nfev=nfev,
        success=False,
        message="maximum iterations exceeded",
        method="newton",
    )


def bisect[
    dtype: DType,
    f: def[dtype: DType](
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
](
    args: Optional[List[Scalar[dtype]]],
    bracket: Tuple[Scalar[dtype], Scalar[dtype]],
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> RootResult[dtype]:
    """Finds a root using the bisection method over a bracket [a, b].

    Requires f(a) and f(b) to have opposite signs. Terminates when the
    interval half-width falls below max(xtol, rtol * |c|), f(c) == 0,
    or maxiter is reached.

    Parameters:
        dtype: The floating-point data type.
        f: Function for which to find a root.

    Args:
        args: Optional arguments forwarded to f.
        bracket: (a, b) tuple where f(a) and f(b) have opposite signs.
        xtol: Absolute tolerance for convergence.
        rtol: Relative tolerance for convergence.
        maxiter: Maximum number of iterations.

    Raises:
        Error: If f(a) and f(b) do not have opposite signs.

    Returns:
        RootResult[dtype] containing the root and convergence information.
    """
    var a: Scalar[dtype] = bracket[0]
    var b: Scalar[dtype] = bracket[1]

    var fa = f(a, args)
    var fb = f(b, args)
    var nfev: Int = 2
    var nit: Int = 0

    if fa == 0:
        return RootResult[dtype](
            root=a,
            nit=0,
            nfev=nfev,
            success=True,
            message="converged",
            method="bisect",
        )
    if fb == 0:
        return RootResult[dtype](
            root=b,
            nit=0,
            nfev=nfev,
            success=True,
            message="converged",
            method="bisect",
        )

    if fa * fb > 0:
        raise Error(
            "SciJo [bisect]: f(a) and f(b) must have opposite signs (bracket"
            " does not enclose a root)."
        )

    for _ in range(maxiter):
        var c = (a + b) / 2
        var fc = f(c, args)
        nfev += 1
        nit += 1

        if fc == 0:
            return RootResult[dtype](
                root=c,
                nit=nit,
                nfev=nfev,
                success=True,
                message="converged",
                method="bisect",
            )

        var tol_x = max(xtol, rtol * abs(c))
        var half_width = abs(b - a) / 2
        if half_width <= tol_x:
            return RootResult[dtype](
                root=c,
                nit=nit,
                nfev=nfev,
                success=True,
                message="converged",
                method="bisect",
            )

        if fa * fc < 0:
            b = c
        else:
            a = c
            fa = fc

    return RootResult[dtype](
        root=(a + b) / 2,
        nit=nit,
        nfev=nfev,
        success=False,
        message="maximum iterations exceeded",
        method="bisect",
    )


def secant[
    dtype: DType,
    f: def[dtype: DType](
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
](
    args: Optional[List[Scalar[dtype]]],
    x0: Scalar[dtype],
    x1: Scalar[dtype],
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> RootResult[dtype]:
    """Finds a root using the secant method.

    A derivative-free method that approximates the derivative using finite
    differences between the two most recent iterates.

    Parameters:
        dtype: The floating-point data type.
        f: Function for which to find a root.

    Args:
        args: Optional arguments forwarded to f.
        x0: First initial guess.
        x1: Second initial guess.
        xtol: Absolute tolerance for convergence.
        rtol: Relative tolerance for convergence.
        maxiter: Maximum number of iterations.

    Raises:
        Error: If zero slope is encountered.

    Returns:
        RootResult[dtype] containing the root and convergence information.
    """
    var a: Scalar[dtype] = x0
    var b: Scalar[dtype] = x1
    var nfev: Int = 0
    var nit: Int = 0

    for _ in range(maxiter):
        var f0 = f(a, args)
        var f1 = f(b, args)
        nfev += 2
        nit += 1

        var denom = f1 - f0
        if denom == 0:
            raise Error(
                "SciJo [secant]: Secant method encountered zero slope (f1 - f0"
                " == 0)."
            )

        var xn = b - (f1 * (b - a)) / denom

        var fxn = f(xn, args)
        nfev += 1
        var tol_x = max(xtol, rtol * abs(xn))
        var tol_f = max(xtol, rtol * abs(fxn))

        if abs(fxn) <= tol_f or abs(xn - b) <= tol_x:
            return RootResult[dtype](
                root=xn,
                nit=nit,
                nfev=nfev,
                success=True,
                message="converged",
                method="secant",
            )

        a = b
        b = xn

    return RootResult[dtype](
        root=b,
        nit=nit,
        nfev=nfev,
        success=False,
        message="maximum iterations exceeded",
        method="secant",
    )
