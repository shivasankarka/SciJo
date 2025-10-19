# SciJo

<div align="center">
  <img src="./assets/scijo.png" alt="SciJo Logo" width="200" style="border-radius: 50%;"/>
  <p><em>A high-performance scientific computing library for Mojo, providing SciPy-like functionality with native Mojo performance.</em></p>
</div>

## Overview

SciJo is a scientific computing library written in pure Mojo that provides essential features similar to SciPy. By leveraging Mojo's performance capabilities and type safety, SciJo delivers familiar scientific computing tools optimized for high-performance applications.

## Features

- **High Performance**: Written in pure Mojo for maximum speed and efficiency
- **SciPy Compatible**: Familiar APIs and function signatures for easy migration
- **Type Safe**: Leverages Mojo's strong type system for reliability
- **Zero Dependencies**: Pure Mojo implementation with NumoJo backend

## Current Modules

### Physical Constants (`scijo.constants`)
- **CODATA 2022 Values**: Complete set of fundamental physical constants
- **100+ Constants**: Atomic units, particle masses, electromagnetic constants
- **SciPy Compatible**: Same structure as `scipy.constants` module

### Numerical Differentiation (`scijo.differentiate`)
- **Finite Difference Methods**: Central, forward, and backward differences
- **Multiple Orders**: Supports orders 1-6 for forward/backward, 2,4,6,8 for central differences

### Integration (`scijo.integrate`)
- **Adaptive Quadrature** (`quad`): QUADPACK QAGS algorithm with 21-point Gauss-Kronrod rules
- **Trapezoidal Rule** (`trapezoid`): Basic numerical integration

### Interpolation (`scijo.interpolate`)
- **Linear Interpolation** (`interp1d`): Basic 1D interpolation functionality

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
Note that SciJo and NuMojo require the `modular` package. We will move to `mojo` only package in future if possible to reduce dependancy size.

2) Install in pixi
```bash
pixi install
```

### Method 2:
1) Clone Repository
```bash
git clone https://github.com/shivasankarka/SciJo.git
cd SciJo
```
2) Build the package
```bash
mojo build scijo
```

3) Move .mojopkg to your project directory
```bash
mv build/scijo.mojopkg /path/to/your/project
```

## Quick Start

### Physical Constants
```mojo
from scijo.constants import physical_constants, value, unit

fn main() raises:
    print("Speed of light:", value("speed_of_light_in_vacuum"), "m/s")
    print("Planck constant:", value("Planck_constant"), unit("Planck_constant"))
```

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

fn my_function(x: Float64, args: List[Float64]) -> Float64:
    return x * x  // f(x) = x²

fn main():
    var args = List[Float64]()
    var result = quad[DType.float64, my_function](0.0, 4.0, args)
    print("Integral value:", result.integral)
```

## Roadmap

- **Linear Algebra**: Matrix operations and decompositions
- **Optimization**: Function minimization and root finding
- **Statistics**: Statistical functions and distributions
- **Signal Processing**: FFT, filtering, and transforms
- **Sparse Matrices**: Efficient sparse matrix operations

## Contributing

We welcome contributions! Areas where help is needed:
- Algorithm implementation and optimization
- Testing and benchmarks
- Documentation and examples

## License

Distributed under the Apache 2.0 License with LLVM Exceptions. See [LICENSE](LICENSE) for more information.

## Citation

```bibtex
@software{scijo,
  author = {Shivasankar K.A. and SciJo Contributors},
  title = {SciJo: High-Performance Scientific Computing in Mojo},
  url = {https://github.com/shivasankarka/SciJo},
  year = {2025}
}
```

---

⚠️ **Note**: This library is in early development and may introduce breaking changes between versions.
