# SciJo

<div align="center">
  <img src="./assets/scijo.png" alt="SciJo Logo" width="200" style="border-radius: 50%;"/>
  
  <p><em>A high-performance scientific computing library for Mojo, providing SciPy-like functionality with the speed and efficiency of native Mojo code.</em></p>
</div>

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Current Modules](#current-modules)
  - [Numerical Differentiation](#numerical-differentiation-scijodifferentiate)
  - [Physical Constants](#physical-constants-scijoconstants)
  - [Integration](#integration-scijointegrate---work-in-progress)
  - [Interpolation](#interpolation-scijointerpolate---early-development)
  - [Coming Soon](#coming-soon)
  - [Under Development](#under-development)
- [Installation](#installation)
- [Quick Start](#quick-start)
  - [Physical Constants Example](#physical-constants-example)
  - [Numerical Differentiation Example](#numerical-differentiation-example)
  - [Integration Example](#integration-example-work-in-progress)
- [Design Philosophy](#design-philosophy)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [Performance Benchmarks](#performance-benchmarks)
- [License](#license)
- [Acknowledgments](#acknowledgments)
- [Citation](#citation)

## Overview

SciJo is a comprehensive scientific computing library written in pure Mojo that aims to provide all the essential features of SciPy. By leveraging Mojo's performance capabilities and zero-cost abstractions, SciJo delivers scientific computing tools that are both familiar to Python users and optimized for high-performance applications.

## Key Features

- 🚀 **High Performance**: Written in pure Mojo for maximum speed and efficiency
- 🔧 **SciPy Compatible**: Familiar APIs and function signatures for easy migration
- 📊 **Comprehensive**: Covers major scientific computing domains
- 🎯 **Type Safe**: Leverages Mojo's strong type system for reliability
- 🔄 **Zero Dependencies**: Pure Mojo implementation with no external dependencies

## Current Modules

### Numerical Differentiation (`scijo.differentiate`)
- **Finite Difference Methods** (`derivative`): Central, forward, and backward finite differences
- **Richardson Extrapolation**: Adaptive step sizing for improved accuracy  
- **Multiple Orders**: Supports orders 1-6 for forward/backward, 2,4,6,8 for central differences
- **SciPy Compatible**: Similar interface to `scipy.derivative`
- **Comprehensive Error Handling**: Detailed validation and informative error messages

### Physical Constants (`scijo.constants`)
- **CODATA 2022 Values**: Complete set of fundamental physical constants
- **SI Units and Prefixes**: Mathematical constants, unit conversions, and SI prefixes
- **SciPy Compatible**: Same structure as `scipy.constants` module
- **100+ Constants**: Atomic units, particle masses, electromagnetic constants, and more
- **Type Safety**: Structured constants with values, units, and uncertainties

### Integration (`scijo.integrate`) - *Work in Progress*
- **Adaptive Quadrature** (`quad`): QUADPACK QAGS algorithm with 21-point Gauss-Kronrod rules
- **Trapezoidal Rule** (`trapezoid`): Basic numerical integration
- *Note: Integration module is still under development and testing*

### Interpolation (`scijo.interpolate`) - *Early Development*
- **Linear Interpolation** (`interp1d`): Basic 1D interpolation functionality
- *Note: Module structure exists but requires NumoJo dependency*

### Coming Soon
- **Linear Algebra** (`scijo.linalg`): Matrix operations, decompositions, and solvers
- **Optimization** (`scijo.optimize`): Function minimization and root finding  
- **Statistics** (`scijo.stats`): Statistical functions and distributions
- **Signal Processing** (`scijo.signal`): Filtering, transforms, and signal analysis
- **Sparse Matrices** (`scijo.sparse`): Efficient sparse matrix operations

### Under Development
- **Integration** (`scijo.integrate`): Completing QUADPACK implementation and testing
- **Interpolation** (`scijo.interpolate`): Expanding beyond basic linear interpolation

## Installation

```bash
# Clone the repository
git clone https://github.com/your-username/scijo.git
cd scijo

# Add to your Mojo project
# (Specific installation instructions will be added as the project matures)
```

## Quick Start

### Physical Constants Example

```mojo
from scijo.constants import physical_constants, value, unit

fn main() raises:
    # Access fundamental constants
    print("Speed of light:", physical_constants["speed_of_light_in_vacuum"].value, "m/s")
    print("Planck constant:", value("Planck_constant"), unit("Planck_constant"))
    print("Elementary charge:", physical_constants["elementary_charge"].value, "C")
    
    # Use in calculations
    var c = physical_constants["speed_of_light_in_vacuum"].value
    var energy = 0.511e6 * physical_constants["electron_volt"].value  # electron rest energy
```

### Numerical Differentiation Example

```mojo
from scijo.differentiate import derivative
from math import sin, cos

fn my_function(x: Float64, args: List[Float64]) -> Float64:
    return sin(x)  # f(x) = sin(x), f'(x) should be cos(x)

fn main() raises:
    var args = List[Float64]()
    var x0 = 1.0  # Point to evaluate derivative
    
    # Compute derivative using central differences
    var result = derivative[DType.float64, my_function](
        x0, args, 
        step_direction="central", 
        order=4  # 4th order accuracy
    )
    
    print("f'(1.0) ≈", result.df)          # ≈ cos(1.0) ≈ 0.5403
    print("Error estimate:", result.error)  # Very small error
    print("Converged:", result.success)     # True if converged
```

### Integration Example (Work in Progress)

```mojo
// Note: Integration module is still being tested and refined
from scijo.integrate.quad import quad

fn my_function(x: Float64, args: List[Float64]) -> Float64:
    return x * x  // f(x) = x²

fn main():
    var args = List[Float64]()
    var result = quad[DType.float64, my_function](0.0, 4.0, args)
    
    print("Integral value:", result.integral)    // ≈ 21.333333  
    print("Error estimate:", result.abserr)      // Small error
    print("Success:", result.ier)                // Integration status
```

## Design Philosophy

### Performance First
SciJo is designed from the ground up for performance, taking advantage of Mojo's:
- Zero-cost abstractions
- Compile-time optimizations
- Memory-efficient data structures
- Parallel execution capabilities

### SciPy Compatibility
While optimized for performance, SciJo maintains familiar APIs:
- Similar function signatures and parameter names
- Consistent return value formats
- Compatible default values and behavior

### Type Safety
Leverages Mojo's type system for:
- Compile-time error detection
- Generic programming with parametric types
- Memory safety without runtime overhead

## Roadmap

### Phase 1 (Current Progress)
- [x] **Physical Constants**: Complete CODATA 2022 dataset with 100+ constants
- [x] **Numerical Differentiation**: Finite differences
- [x] **Core Infrastructure**: Module structure and type-safe design
- [ ] **Integration**: Finalizing QUADPACK QAGS implementation and testing
- [ ] **Documentation**: API docs and comprehensive examples

### Phase 2 (Next Steps)
- [ ] **Linear Algebra**: Basic matrix operations and decompositions
- [ ] **Optimization**: Root finding and function minimization
- [ ] **Advanced Integration**: Multi-dimensional quadrature
- [ ] **Interpolation**: Spline and polynomial methods

### Phase 3 (Future)
- [ ] **Statistics**: Probability distributions and statistical functions
- [ ] **Signal Processing**: FFT, filtering, and transforms
- [ ] **Sparse Matrices**: Efficient sparse matrix operations
- [ ] **Performance Optimization**: SIMD and parallelization

### Phase 4 (Advanced Features)
- [ ] **Machine Learning**: Utilities for ML algorithms
- [ ] **Advanced Linear Algebra**: SVD, eigendecomposition
- [ ] **Specialized Integration**: Oscillatory and singular integrals
- [ ] **High-Performance Computing**: GPU acceleration

## Contributing

We welcome contributions to SciJo! Areas where help is needed:

- **Algorithm Implementation**: Porting SciPy algorithms to Mojo
- **Performance Optimization**: Leveraging Mojo's performance features
- **Testing**: Comprehensive test suites and benchmarks
- **Documentation**: Examples, tutorials, and API documentation

### Getting Started with Development

1. Fork the repository
2. Create a feature branch
3. Implement your changes with tests
4. Ensure compatibility with SciPy interfaces
5. Submit a pull request

## Performance Benchmarks

*Note: Comprehensive benchmarks are in development. Initial tests show promising results:*

| Operation | Status | Notes |
|-----------|--------|-------|
| Physical Constants Access | ✅ Implemented | Zero-cost compile-time constants |
| Numerical Differentiation | ✅ Implemented | Richardson extrapolation with adaptive stepping |
| Quadrature Integration | 🔄 Testing | QUADPACK implementation under validation |
| *More benchmarks coming soon* | | |

*Detailed performance comparisons with SciPy will be published once modules are fully validated.*

## License

SciJo is released under the [MIT License](LICENSE).

## Acknowledgments

- **SciPy Team**: For the incredible foundation and API design
- **Modular Team**: For creating the Mojo programming language

## Citation

If you use SciJo in your research, please cite:

```bibtex
@software{scijo,
  author = {Shivasankar K.A. and SciJo Contributors},
  title = {SciJo: High-Performance Scientific Computing in Mojo},
  url = {https://github.com/your-username/scijo},
  year = {2025}
}
```

---

**Note**: SciJo is under active development. APIs may change as the library matures. Please check the documentation for the latest updates.
