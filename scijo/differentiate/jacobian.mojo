from numojo.routines.creation import zeros, full
from numojo.core import NDArray, Shape

from algorithm.functional import parallelize

fn jacobian[
    dtype: DType,
    f: fn [dtype: DType] (x: NDArray[dtype], args: Optional[List[Scalar[dtype]]]) raises -> NDArray[dtype],
](
    x: NDArray[dtype],
    args: Optional[List[Scalar[dtype]]] = None,
    tolerances: Dict[String, Scalar[dtype]] = {"abs": 1e-5, "rel": 1e-3},
    maxiter: Int = 10
) raises -> NDArray[dtype]:
    var n: Int = len(x)
    var f0: NDArray[dtype] = f(x, args)
    var m: Int = len(f0)

    var jacob: NDArray[dtype] = zeros[dtype](Shape(m, n))
    var step: NDArray[dtype] = full[dtype](Shape(n), fill_value=0.5)

    # For each input dimension j, perturb x at index j and compute the
    # corresponding column of the Jacobian: (f(x+h e_j) - f(x-h e_j)) / (2 h_j)
    @parameter
    fn closure(j: Int):
        try:
            var x_plus  = x.copy()
            var x_minus = x.copy()
            var hj = step.load(j)
            x_plus.store(j, val=x.load(j) + hj)
            x_minus.store(j, val=x.load(j) - hj)
            var f_plus  = f(x_plus, args)
            var f_minus = f(x_minus, args)

            var col: NDArray[dtype] = (f_plus - f_minus) / (2.0 * hj)

            for i in range(m):
                var val: Scalar[dtype] = col.load(i)
                var flat_idx: Int = i * n + j
                jacob.store(flat_idx, val=val)
        except:
            print("Error in computing Jacobian column for j =", j)

    parallelize[closure](n)

    # naive algorithm
    # for j in range(m):
    #     var x_plus  = x.copy()
    #     var x_minus = x.copy()
    #     var hj = step.load(j)
    #     x_plus.store(j, val=x.load(j) + hj)
    #     x_minus.store(j, val=x.load(j) - hj)
    #     var f_plus  = f(x_plus)
    #     var f_minus = f(x_minus)

    #     for i in range(n):
    #         var val: Scalar[dtype] = (f_plus.item(i) - f_minus.item(i)) / (2.0 * hj)
    #         var flat_idx: Int = i * n + j
    #         jacob.store(flat_idx, val=val)
    #
    return jacob^
