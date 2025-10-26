from scijo.differentiate.derivative import derivative
import scijo as sj
from utils import StaticTuple
import benchmark
from benchmark.compiler import keep
from numojo.prelude import *
from algorithm.functional import vectorize
from sys.info import simd_width_of

alias dtype: DType = DType.float64
alias SIZE: Int = 10
alias coeffs: List[Scalar[dtype]] = List[Scalar[dtype]](
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
    31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0,
    41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0
)
alias simd_width: Int = simd_width_of[dtype]()

fn poly_func[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
    var result: Scalar[dtype] = 0.0
    var coeff = materialize[coeffs]()
    var coeff_ptr = coeff.unsafe_ptr()
    # test only for 10 coeffs.
    @parameter
    fn closure[wid: Int](i: Int):
        result += SIMD[dtype, wid](coeff_ptr.load[width=wid](i)).reduce_add() * (x ** i)
    vectorize[closure, simd_width](SIZE)

    # the speed is not so different with native loop.
    # for i in range(SIZE):
    #     result = result + Scalar[dtype](coeff[i]) * x ** i
    return result

alias x0 = 1.0
def numeric_derivative_1point[order: Int]():
    var result = derivative[dtype, poly_func, step_direction=0](
        x0=x0,
        args=None,
        order=order,
    )
    keep(result.df)

def numeric_derivative_npoints[order: Int]():
    var xs: NDArray[dtype] = nm.linspace[dtype](-1.0, 1.0, 20)
    for i in range(20):
        var result = derivative[dtype, poly_func, step_direction=0](
            x0=xs._buf.ptr[i],
            args=None,
            order=order,
        )
        keep(result.df)

fn main() raises:
    # var result = derivative[dtype, simple_func, step_direction=0](
    #     x0=1.0,
    #     args=None,
    #     order=8,
    # )
    # print("Computed derivative: ", result)

    # testing only order 8 for now.
    alias ord2 = 8
    var report2 = benchmark.run[numeric_derivative_1point[ord2]]()
    print("Derivative Benchmark - 1 Point Evaluation")
    print("order=: ", ord2)
    report2.print()

    var report3 = benchmark.run[numeric_derivative_npoints[ord2]]()
    print("Derivative Benchmark - N Points Evaluation")
    print("order=: ", ord2)
    report3.print()
