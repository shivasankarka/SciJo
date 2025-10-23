fn root_scalar[
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
    x1: Optional[Scalar[dtype]] = None,
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
    # options: SolverOptions
) raises -> Scalar[dtype]:
    @parameter
    if method == "newton":
        return _newton_raphson[dtype, f, fprime](
            args, x0, x1, xtol, rtol, maxiter
        )
    else:
        raise Error("Unsupported method: " + String(method))


fn _newton_raphson[
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
    x1: Optional[Scalar[dtype]] = None,
    xtol: Scalar[dtype] = 1e-8,
    rtol: Scalar[dtype] = 1e-8,
    maxiter: Int = 100,
) raises -> Scalar[dtype]:
    var xn: Scalar[dtype]
    if x0:
        xn = x0.value()
    else:
        raise Error("Initial guess x0 must be provided for Newton's method.")

    for _ in range(maxiter):
        var fx = f(xn, args)
        var fpx = fprime(xn, args)
        xn = xn - (fx / fpx)

    return xn
