from utils import StaticTuple

fn get_quad_error_message(ier: Int) -> String:
    """
    Get error message for QUADPACK integration error codes.

    Arguments:
        ier: Integration error code (0 = success, >0 = error type)
    """
    if ier == 0:
        return String("The integral converged successfully.")
    elif ier == 1:
        return String("Maximum subdivisions reached. The integrand may have singularities or discontinuities. Try splitting the integration interval at problem points or using a specialized integrator.")
    elif ier == 2:
        return String("Roundoff error prevents achieving the requested tolerance. The error estimate may be inaccurate.")
    elif ier == 3:
        return String("The function has problematic behavior at some points in the integration interval.")
    elif ier == 4:
        return String("Algorithm failed to converge due to roundoff errors. The result may be the best possible approximation.")
    elif ier == 5:
        return String("The integral is likely divergent or converges very slowly.")
    elif ier == 6:
        return String("Invalid input parameters: tolerance values or limits are out of valid range.")
    else:
        return String("Unknown error code.")

struct IntegralResult[dtype: DType](Copyable, Movable, Writable):
    """Result structure for numerical integration operations.

    This structure encapsulates the results of integration computations, including
    the computed integral value, error estimates, and diagnostic data for numerical
    analysis following the QUADPACK convention.

    Type Parameters:
        dtype: The floating-point data type (DType.float32, DType.float64, etc.)

    Fields:
        integral: The computed integral value
        abserr: Absolute error estimate
        neval: Number of function evaluations used
        ier: Integration error code (0 = success, >0 = error type)
        last: Number of subintervals used in adaptive algorithm

    Usage:
        Creates a result object containing integration computation results
        with error estimates and diagnostic information.
    """

    var integral: Scalar[dtype]
    var abserr: Scalar[dtype]
    var neval: Int
    var ier: Int  # Changed from Bool to Int for QUADPACK error codes
    var last: Int

    fn __init__(
        out self,
        integral: Scalar[dtype] = Scalar[dtype](0),
        abserr: Scalar[dtype] = Scalar[dtype](0),
        neval: Int = 0,
        ier: Int = 0,
        last: Int = 0,
    ):
        self.integral = integral
        self.abserr = abserr
        self.neval = neval
        self.ier = ier
        self.last = last

    fn success(self) -> Bool:
        """Check if integration was successful."""
        return self.ier == 0

    fn message(self) -> String:
        """Get SciPy-style error message."""
        return get_quad_error_message(self.ier)

    fn __str__(self) raises -> String:
        var status = "SUCCESS" if self.success() else "ERROR"
        return String(
            "QuadResult(status={}, message='{}', integral={}, abserr={:.2e},"
            " neval={}, last={})"
        ).format(
            status,
            self.message(),
            self.integral,
            self.abserr,
            self.neval,
            self.last,
        )

    fn write_to[W: Writer](self, mut writer: W):
        try:
            var status = "SUCCESS" if self.success() else "ERROR"
            writer.write(
                String(
                    "Numerical Integration Result\n"
                    + "============================\n"
                    + "Status      : {}\n"
                    + "Message     : {}\n"
                    + "Integral    : {}\n"
                    + "Error Est.  : {}\n"
                    + "Func Evals  : {}\n"
                    + "Subintervals: {}\n"
                ).format(
                    status,
                    self.message(),
                    self.integral,
                    self.abserr,
                    self.neval,
                    self.last,
                )
            )
        except e:
            writer.write("Error displaying QuadResult: " + String(e) + "\n")


# Nodes and weights for G10K21 rule (from Advanpix)
# Structure: 21 Kronrod points total, with 11 Gauss points embedded within
alias kronrod_nodes: StaticTuple[Float64, 11] = StaticTuple[Float64, 11](
    # All 21 Kronrod nodes in ascending order (symmetric about 0)
    0.0000000000000000000000000000000000,  # 0: center (Gauss)
    0.1488743389816312108848260011297200,  # 1: ± (Gauss)
    0.2943928627014601981311266031038656,  # 2: ± (Kronrod only)
    0.4333953941292471907992659431657842,  # 3: ± (Gauss)
    0.5627571346686046833390000992726941,  # 4: ± (Kronrod only)
    0.6794095682990244062343273651148736,  # 5: ± (Gauss)
    0.7808177265864168970637175783450424,  # 6: ± (Kronrod only)
    0.8650633666889845107320966884234930,  # 7: ± (Gauss)
    0.9301574913557082260012071800595083,  # 8: ± (Kronrod only)
    0.9739065285171717200779640120844521,  # 9: ± (Gauss)
    0.9956571630258080807355272806890028,  # 10: ± (Kronrod only)
)

alias kronrod_weights: StaticTuple[Float64, 11] = StaticTuple[Float64, 11](
    0.1494455540029169056649364683898212,  # 0: center
    0.1477391049013384913748415159720680,  # 1: ± (Gauss node)
    0.1427759385770600807970942731387171,  # 2: ± (Kronrod only)
    0.1347092173114733259280540017717068,  # 3: ± (Gauss node)
    0.1234919762620658510779581098310742,  # 4: ± (Kronrod only)
    0.1093871588022976418992105903258050,  # 5: ± (Gauss node)
    0.0931254545836976055350654634063634,  # 6: ± (Kronrod only)
    0.0750396748109199527670431409161900,  # 7: ± (Gauss node)
    0.0547558965743519960313813024458018,  # 8: ± (Kronrod only)
    0.0325581623079647274788190031921602,  # 9: ± (Gauss node)
    0.0116946388673718742780643960621925,  # 10: ± (Kronrod only)
)

alias gauss_weights: StaticTuple[Float64, 5] = StaticTuple[Float64, 5](
    # Gauss-Legendre weights for 10-point rule embedded in G10K21 (from official table n=10)
    0.2955242247147529,  # ±0.1488743389816312
    0.2692667193099963,  # ±0.4333953941292472
    0.2190863625159820,  # ±0.6794095682990244
    0.1494513491505806,  # ±0.8650633666889845
    0.0666713443086881,  # ±0.9739065285171717
)
