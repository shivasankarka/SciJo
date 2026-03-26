# SciJo

<div align="center">
    <img src="./assets/scijo.png" alt="SciJo Logo" width="200" style="border-radius: 32px; margin-bottom: 200px; display: block; border: 3px solid rgba(0, 0, 0, 0.15); box-shadow: 0 20px 45px rgba(0, 0, 0, 0.25); background: #fff;"/>
  <p style="font-size: 1.2em; color: #666; margin: 0; padding: 10px 20px; line-height: 1.5;">
    <em>High-performance scientific computing library for Mojo, written in pure Mojo, inspired by SciPy</em>
  </p>
</div>

**[Changelog»](https://github.com/shivasankarka/SciJo/tree/main/docs/changelog.md)**  |  **[NuMojo Docs»](https://numojo.readthedocs.io)**

## Overview

SciJo brings SciPy-like numerical tools to Mojo. It is built in pure Mojo on top of **[NuMojo](https://github.com/Mojo-Numerics-and-Algorithms-group/NuMojo)** and focuses on fast, type-safe scientific computing.

## Modules (current)

- **`scijo.differentiate`**: finite-difference derivatives (`derivative`, `jacobian`)
- **`scijo.integrate`**: quadrature and fixed-sample rules (`quad`, `trapezoid`, `simpson`, `romb`)
- **`scijo.interpolate`**: 1D interpolation (`interp1d`)
- **`scijo.fft`**: FFT/IFFT for complex arrays
- **`scijo.constants`**: CODATA physical constants
- **`scijo.optimize`**: Optimizers. 

## Installation

### Pixi (recommended)

1) Add to `pixi.toml`:

  ```toml
    [workspace]
    preview = ["pixi-build"]

    [dependencies]
    modular = ">=25.6.1,<26"
    scijo = { git = "https://github.com/mojomath/SciJo.git", branch = "main" }
  ```

2) Install:
  ```console
    pixi install
  ```

### Build from source
  ```console
    git clone https://github.com/mojomath/SciJo.git
    cd SciJo
    mojo build scijo
    mv build/scijo.mojopkg /path/to/your/project
  ```

## Quick Start

### Derivative
  ```mojo
    import scijo as sj
    from scijo.differentiate import derivative

    fn f[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
        return x * x

    fn main() raises:
        var res = derivative[sj.f64, f](1.0)
        print("df:", res.df)
  ```
  
### Integration (quad)
  ```mojo
    import scijo as sj
    from scijo.integrate import quad

    fn f[dtype: DType](x: Scalar[dtype], args: Optional[List[Scalar[dtype]]]) -> Scalar[dtype]:
        return x * x

    fn main() raises:
        var res = quad[sj.f64, f](0.0, 1.0, None)
        print("integral:", res.integral)
  ```
  
### FFT / IFFT
  ```mojo
    import numojo as nm
    from scijo.fft import fft, ifft
    from scijo.prelude import *

    fn main() raises:
        var arr = nm.arange[cf32](CScalar[f32](0), CScalar[f32](8))
        var y_fft = fft[cf32](arr)
        var y_ifft = ifft[cf32](y_fft)
        print(y_fft)
        print(y_ifft)
  ```

## Roadmap (short)

- More integrators (QAGSE, additional adaptive methods)
- Real FFT (`rfft`, `irfft`) and 2D FFT
- Additional interpolation methods (cubic, spline)
- Expanded optimization tools

## Contributing

Contributions are welcome. Focus areas:
- New algorithms
- Performance improvements
- Tests and docs

## License

Distributed under the Apache 2.0 License with LLVM Exceptions. See [LICENSE](LICENSE) for details.

## Citation

Feel free to cite SciJo in your work.
```tex
    @software{scijo,
      author = {Shivasankar K.A. and SciJo Contributors},
      title = {SciJo: High-Performance Scientific Computing in Mojo},
      url = {https://github.com/shivasankarka/SciJo},
      year = {2025}
    }
```

---

⚠️ **Note**: This library is in early development and may introduce breaking changes between versions.
