from python import Python, PythonObject
from testing.testing import assert_true
import numojo as nm

fn check_is_close[
    dtype: DType
](
    scalar_1: Scalar[dtype],
    scalar_2: Scalar[dtype],
    st: String,
    rtol: Scalar[dtype] = 1e-6,
) raises:
    assert_true(scalar_1 - scalar_2 <= rtol, st)
