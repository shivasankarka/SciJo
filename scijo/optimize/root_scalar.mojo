# TODO: check if we are using the right tolerance conditions in all methods.

fn root_scalar[
    dtype: DType,
    f: fn[dtype: DType] (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    fprime: Optional[fn[
            dtype: DType
        ] (x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[
            dtype
        ]] = None,
    method: String = "bisect",
](
    args: Optional[List[Scalar[dtype]]] = None,
    x0: Optional[Scalar[dtype]] = None,
    x1: Optional[Scalar[dtype]] = None,
    bracket: Optional[Tuple[Scalar[dtype], Scalar[dtype]]] = None,
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
    # options: SolverOptions
) raises -> Scalar[dtype]:
    """
    Find a root of a scalar function.

    Parameters
    - f: function f(x, args) that returns a scalar of dtype.
    - fprime: derivative function f'(x, args). Required for Newton method.
    - method: "newton" or "bisect".
    - args: optional additional arguments forwarded to f and fprime.
    - x0: initial guess for Newton's method (required if method == "newton").
    - bracket: (a, b) tuple with a and b such that f(a) and f(b) have opposite signs for bisection.
    - xtol, rtol: absolute and relative tolerances for termination.
    - maxiter: maximum number of iterations.

    Raises
    - Error if required inputs for the chosen method are missing or invalid.
    """

    @parameter
    if method == "newton" and fprime:
            return newton[dtype, f, fprime.value()](args, x0, xtol, rtol, maxiter)
    elif method == "bisect":
        if not bracket:
            raise Error("Bracket must be provided for bisection method.")
        return bisect[dtype, f](args, bracket.value(), xtol, rtol, maxiter)
    elif method == "secant":
        if not (x0 and x1):
            raise Error("Initial guesses x0 and x1 must be provided for secant method.")
        return secant[dtype, f](args, x0.value(), x1.value(), xtol, rtol, maxiter)
    else:
        raise Error("Unsupported method: " + String(method))


fn newton[
    dtype: DType,
    f: fn[dtype: DType] (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    fprime: fn[dtype: DType] (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
    method: String = "newton",
](
    args: Optional[List[Scalar[dtype]]],
    x0: Optional[Scalar[dtype]] = None,
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> Scalar[dtype]:
    """
    Newton-Raphson method for finding a root.

    Terminates when either:
    - the change in x is below max(xtol, rtol * |x|), or
    - |f(x)| is below max(xtol, rtol * |x|),
    or when maxiter is reached.

    Raises Error if initial guess x0 is not provided or derivative is zero.
    """
    var xn: Scalar[dtype]
    if x0:
        xn = x0.value()
    else:
        raise Error("Initial guess x0 must be provided for Newton's method.")

    for _ in range(maxiter):
        var fx = f(xn, args)
        var fpx = fprime(xn, args)

        if fpx == 0:
            raise Error(
                "Derivative is zero. Newton-Raphson step would divide by zero."
            )

        var delta = fx / fpx
        var xn_next = xn - delta

        var tol_x = max(xtol, rtol * abs(xn_next))
        var tol_f = max(xtol, rtol * abs(fx))

        if abs(fx) <= tol_f:
            return xn_next
        if abs(delta) <= tol_x:
            return xn_next

        xn = xn_next

    return xn


fn bisect[
    dtype: DType,
    f: fn[dtype: DType] (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
](
    args: Optional[List[Scalar[dtype]]],
    bracket: Tuple[Scalar[dtype], Scalar[dtype]],
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> Scalar[dtype]:
    """
    Bisection method for finding a root in a bracket [a, b].

    Requirements:
    - f(a) and f(b) must have opposite signs (f(a) * f(b) < 0),
      otherwise there is no guarantee of a root in the interval.

    Terminates when:
    - the interval half-width (b - a) / 2 is below max(xtol, rtol * |c|), or
    - f(c) == 0, or when maxiter is reached.
    """
    var a: Scalar[dtype] = bracket[0]
    var b: Scalar[dtype] = bracket[1]

    var fa = f(a, args)
    var fb = f(b, args)

    if fa == 0:
        return a
    if fb == 0:
        return b

    if fa * fb > 0:
        raise Error(
            "f(a) and f(b) must have opposite signs (bracket does not enclose a"
            " root)."
        )

    for _ in range(maxiter):
        var c = (a + b) / 2
        var fc = f(c, args)

        if fc == 0:
            return c

        var tol_x = max(xtol, rtol * abs(c))
        var half_width = (b - a) / 2
        if half_width <= tol_x:
            return c

        if fa * fc < 0:
            # Root is in [a, c]
            b = c
        else:
            # Root is in [c, b]
            a = c

    return (a + b) / 2

fn secant[
    dtype: DType,
    f: fn[dtype: DType] (
        x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]
    ) -> Scalar[dtype],
](
    args: Optional[List[Scalar[dtype]]],
    x0: Scalar[dtype],
    x1: Scalar[dtype],
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) -> Scalar[dtype]:
    """
    Secant method for finding a root.

    Parameters:
        dtype: Data type of the scalar function.
        f: Function for which to find the root.

    Args:
        args: Optional additional arguments for the function.
        x0: First initial guess.
        x1: Second initial guess.
        xtol: Absolute tolerance for convergence.
        rtol: Relative tolerance for convergence.
        maxiter: Maximum number of iterations.

    Returns:
        Approximation of the root of the function.
    """
    var a: Scalar[dtype] = x0
    var b: Scalar[dtype] = x1

    for _ in range(maxiter):
        var f0 = f[dtype](a, args)
        var f1 = f[dtype](b, args)

        var xn = b - (f1 * (b - a)) / (f1 - f0)

        var fxn = f[dtype](xn, args)
        var tol_x = max(xtol, rtol * abs(xn))
        var tol_f = max(xtol, rtol * abs(fxn))

        if abs(fxn) <= tol_f:
            return xn

        if abs(xn - b) <= tol_x:
            return xn

        b = xn
        a = x1

    return x1
