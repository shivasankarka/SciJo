# SciJo

<div align="center">
    <img src="./assets/scijo.png" alt="SciJo Logo" width="200" style="border-radius: 32px; margin-bottom: 200px; display: block; border: 3px solid rgba(0, 0, 0, 0.15); box-shadow: 0 20px 45px rgba(0, 0, 0, 0.25); background: #fff;"/>
  <p style="font-size: 1.2em; color: #666; margin: 0; padding: 10px 20px; line-height: 1.5;">
    <em>High-performance scientific computing library for Mojo, written in pure Mojo, inspired by SciPy</em>
  </p>
</div>

**[Changelog»](https://github.com/shivasankarka/SciJo/tree/main/docs/changelog.md)**

## Overview

SciJo is a high-performance scientific computing library for Mojo that brings the power and familiarity of SciPy to the Mojo ecosystem. Written in pure Mojo and built on top of **[NuMojo](https://github.com/Mojo-Numerics-and-Algorithms-group/NuMojo)**, SciJo combines the performance benefits of native compilation with the type safety guarantees of Mojo's advanced type system.

## Features

- **Pure Mojo**: Native implementation for maximum performance
- **Familiar APIs**: SciPy-inspired interfaces for easy adoption
- **Type Safe**: Compile-time guarantees with Mojo's type system
- **NuMojo Backend**: Efficient array operations and complex number support

## Current Modules

### Numerical Differentiation (`scijo.differentiate`)
Accurate derivatives using finite difference methods:
- **`derivative`**: Central, forward, and backward differences
- **`jacobian`**: Jacobian matrix computation for vector-valued functions
- **Order control**: Specify accuracy order (1-6 for forward/backward, 2-8 for central)
- **Adaptive stepping**: Automatic step size refinement with Richardson extrapolation
- **Error estimation**: Built-in convergence tracking

### Integration (`scijo.integrate`)
Numerical integration with adaptive algorithms:
- **`quad`**: (QUADPACK QNG algorithm)
  - Successively increasing precision levels (10, 21, 43, 87 point rules)
- **`trapezoid`**: Composite trapezoidal rule for uniform or non-uniform grids
- **`simpson`**: Simpson's rule for discrete data
- **`romb`**: Romberg integration with Richardson extrapolation

### Interpolation (`scijo.interpolate`)
1D data interpolation:
- **`interp1d`**: Linear interpolation (functional and callable interfaces)
- **`LinearInterpolator`**: Reusable callable interpolation object
- Handles both extrapolation and boundary fill methods
- Compatible with NuMojo arrays

### FFT (`scijo.fft`)
Fast Fourier Transform operations:
- **`fft`**: Forward FFT using Cooley-Tukey algorithm
- **`ifft`**: Inverse FFT with proper normalization
- Supports complex arrays (power-of-2 sizes)
- Compatible with NumPy's FFT conventions

### Physical Constants (`scijo.constants`)
Access fundamental physical constants from CODATA 2022:
- Comprehensive physical constants with values, units, and uncertainties
- Mathematical constants (pi, golden ratio, etc.)
- SI prefixes, binary prefixes, and unit conversions
- Helper functions: `value()`, `unit()`, `precision()`, `find()`
- Temperature conversion utilities
- Compatible with `scipy.constants` structure

### Optimization (`scijo.optimize`)
Scalar root-finding and minimization:
- **`root_scalar`**: Unified interface for root finding
  - **`newton`**: Newton-Raphson method
  - **`bisect`**: Bisection method
  - **`secant`**: Secant method (derivative-free)
- **`minimize_scalar`**: Scalar function minimization
  - Brent's method, golden section search, bounded minimization

## Installation

### Method 1:
1) Add to pixi.toml
```toml
[workspace]
preview = ["pixi-build"]

[dependencies]
modular = ">=25.6.1,<26"
scijo = { git = "https://github.com/shivasankarka/SciJo.git", branch = "main"}
```
Note that SciJo and NuMojo require the `modular` package. We will move to `mojo` only package in future if possible.

2) Install in pixi
```bash
pixi install
```

### Method 2: Build from Source
```bash
# Clone and build
git clone https://github.com/shivasankarka/SciJo.git
cd SciJo
mojo build scijo

# Move package to your project
mv build/scijo.mojopkg /path/to/your/project
```

## Quick Start

### Numerical Differentiation
```mojo
import scijo as sj
from scijo.differentiate import derivative

fn simple_function[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]] = None) -> Scalar[dtype]:
    var a = args.value()[0]
    return a * x * x + 2.0 * x + 1.0

fn main() raises:
    var result = derivative[sj.f64, simple_function, step_direction=0](
        x0=1.0,
        args=List[Scalar[sj.f64]](2.0),
        tolerance={"atol": 1e-8, "rtol": 1e-8},
        order=6
    )
    print("Derivative result:", result)
```

### Integration
```mojo
from scijo.integrate.quad import quad

fn simple_function[
    dtype: DType
](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]] = None) -> Scalar[
    dtype
]:
    """A simple function for testing."""
    var a = args.value()[0]
    return a * x * x + 2.0 * x + 1.0

fn main():
    var result = quad[sj.f64, simple_function](
        a=0.0,
        b=1.0,
        args=List[Scalar[sj.f64]](2.0),
        epsabs=1e-6,
        epsrel=1e-6,
    )
    print("Integral value:", result.integral)
```

### Interpolation
```mojo
from scijo.interpolate.interp1d import interp1d
import numojo as nm

fn main() raises:
    var x = nm.arange[nm.f64](0, 5, 1)
    var y = x * x
    var xi = nm.linspace[nm.f64](0.5, 3.5, 4)

    var yi = interp1d[nm.f64, type="linear", fill_method="interpolate"](xi, x, y)
    print("Interpolated values:", yi)
```

### FFT
```mojo
from scijo.fft import fft, ifft
import numojo as nm

fn main() raises:
    # Create complex array
    var arr = nm.arange[nm.cf64](nm.CScalar[nm.cf64](0), nm.CScalar[nm.cf64](8))

    # Forward FFT
    var y_fft = fft[nm.cf64](arr)
    print("FFT result:", y_fft)

    # Inverse FFT
    var y_ifft = ifft[nm.cf64](y_fft)
    print("IFFT result:", y_ifft)
```

### Physical Constants
```mojo
from scijo.constants import physical_constants, value, unit

fn main() raises:
    print("Speed of light:", value("speed_of_light_in_vacuum"), "m/s")
    print("Planck constant:", value("Planck_constant"), unit("Planck_constant"))
```

### Optimization
```mojo
from scijo.optimize import root_scalar, minimize_scalar

fn f[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
    return x * x - 2

fn objective[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
    return (x - 2) * (x - 2) + 1

fn main() raises:
    # Root finding
    var root = root_scalar[f64, f](bracket=(1.0, 2.0), method="bisect")
    print("Root:", root)

    # Minimization
    var result = minimize_scalar[Float64, objective, method="Brent"](
        Bracket=(0.0, 4.0),
        tol=1e-8,
        maxiter=100
    )
    print("Minimum at:", result.x)
```


## Roadmap

### Near Term
- More integration algorithms (QAGSE etc)
- Real FFT (`rfft`, `irfft`) and 2D FFT support
- Additional interpolation methods (cubic, spline)
- Expand differentiation module

### Future
- **Optimization**: Multi-dimensional minimization, curve fitting
- **Signal Processing**: Filtering, windowing, convolution
- **Linear Algebra**: Matrix decompositions (SVD, QR, Cholesky)
- **Sparse Matrices**: Efficient storage and operations

## Contributing

Contributions are most welcome! Feel free to add a functionality and open a PR!

Priority areas:
- Algorithm implementations (see Roadmap)
- Performance benchmarks and optimization
- Tests and documentation
- Bug reports and feature requests

## License

Distributed under the Apache 2.0 License with LLVM Exceptions. See [LICENSE](LICENSE) for more information.

## Citation
Feel free to cite SciJo in your work, helps with visibility :)
```bibtex
@software{scijo,
  author = {Shivasankar K.A. and SciJo Contributors},
  title = {SciJo: High-Performance Scientific Computing in Mojo},
  url = {https://github.com/shivasankarka/SciJo},
  year = {2026}
}
```

---

⚠️ **Note**: This library is in early development and may introduce breaking changes between versions.
