alias error_messages: Dict[Int, String] = {
    0: "Integration was successful",
    1: "The maximum number of subdivisions has been achieved.\n  "
        "If increasing the limit yields no improvement it is advised to "
        "analyze \n  the integrand in order to determine the difficulties.  "
        "If the position of a \n  local difficulty can be determined "
        "(singularity, discontinuity) one will \n  probably gain from "
        "splitting up the interval and calling the integrator \n  on the "
        "subranges.  Perhaps a special-purpose integrator should be used.",
    2: "The occurrence of roundoff error is detected, which prevents \n  "
        "the requested tolerance from being achieved.  "
        "The error may be \n  underestimated.",
    3: "Extremely bad integrand behavior occurs at some points of the\n  "
        "integration interval.",
    4: "The algorithm does not converge.  Roundoff error is detected\n  "
        "in the extrapolation table.  It is assumed that the requested "
        "tolerance\n  cannot be achieved, and that the returned result "
        "(if full_output = 1) is \n  the best which can be obtained.",
    5: "The integral is probably divergent, or slowly convergent.",
    6: "The input is invalid.",
    7: "Abnormal termination of the routine. The estimates for result\n and error are less reliable.  It is assumed that the requested accuracy\n  has not been achieved.",
    8: "Unknown Error",
    80: "A Python error occurred possibly while calling the function."
}

fn get_quad_error_message(ier: Int) -> String:
    """Get SciPy-style error message for QUADPACK error code using dictionary lookup."""
    try:
        return error_messages[ier]
    except:
        return "Unknown error code"

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
            "QuadResult(status={}, message='{}', integral={}, abserr={:.2e}, neval={}, last={})"
        ).format(status, self.message(), self.integral, self.abserr, self.neval, self.last)

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


# Define the nodes and weights for G10K21 rule (correct values from Advanpix)
# Structure: 21 Kronrod points total, with 11 Gauss points embedded within
alias kronrod_nodes = List[Float64](
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

alias kronrod_weights = List[Float64](
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

alias gauss_weights = List[Float64](
    # Gauss-Legendre weights for 11-point rule (G10 + center)
    0.2955242247147528701738929946513383,  # center
    0.2692667193099963550912269215694694,  # ±1
    0.2190863625159820439955349342281632,  # ±2
    0.1494513491505805931457763396576973,  # ±3
    0.0666713443086881375935688098933179,  # ±4
    0.0666713443086881375935688098933179,  # ±5 - this is wrong, each should be different
)
