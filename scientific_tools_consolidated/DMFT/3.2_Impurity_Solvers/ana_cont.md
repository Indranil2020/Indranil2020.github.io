# ana_cont

## Official Resources
- Homepage: https://github.com/josefkaufmann/ana_cont
- Documentation: https://josefkaufmann.github.io/ana_cont/
- Source Repository: https://github.com/josefkaufmann/ana_cont
- License: MIT License
- Installation: `pip install ana_cont`

## Overview
`ana_cont` is a Python package dedicated to the analytic continuation of many-body Green's functions. It provides a user-friendly interface to standard continuation methods, specifically Padé approximants and the Maximum Entropy Method (MaxEnt), enabling the extraction of real-frequency spectral functions from imaginary-axis QMC data.

**Scientific domain**: Analytic Continuation, Many-body Physics
**Target user community**: Researchers converting Euclidean QMC data to real-frequency spectra

## Theoretical Methods
- **Maximum Entropy Method (MaxEnt)**: Probabilistic approach using Bayesian inference to find the most probable spectrum.
- **Padé Approximants**: Rational polynomial interpolation for continuation.
- **Analytic Continuation**: Inversion of $G(i\omega_n) = \int \frac{A(\omega)}{i\omega_n - \omega} d\omega$.

## Capabilities (CRITICAL)
- **MaxEnt Solver**: Robust MaxEnt implementation (`MaxEntSolver` class) for spectral reconstruction.
- **Padé Solver**: Fast rational approximation for lower-noise data (`PadeSolver` class).
- **Kernel Support**: Supports Fermionic and Bosonic kernels (Temperature dependence).
- **Preprocessing**: Tools for handling input data formats and preblur.

## Key Features

### Python Interface:
- Easy to use Python API (OOP design).
- Integrates seamlessly with NumPy/SciPy arrays.
- Scriptable for batch processing of many self-energies.

### Reliability:
- Implements standard, well-tested algorithms (Bryan's method / Historic MaxEnt).
- Configurable parameters (alpha, blur, error tolerance) for MaxEnt.
- GUI available for interactive parameter tuning.

## Inputs & Outputs
- **Inputs**:
  - Matsubara Green's functions $G(i\omega_n)$ (complex numpy array).
  - Error estimates (scalar or covariance).
  - Real frequency grid $\omega$.
- **Outputs**:
  - Real-frequency Spectral Function $A(\omega)$.
  - Trace of fit quality ($\chi^2$).

## Workflow and Usage
Standard post-processing tool.
1.  Load $G(i\omega_n)$ from a file (e.g., TRIQS or w2dynamics output).
2.  Instantiate `MaxEntSolver` or `PadeSolver`.
3.  Run the solver (e.g., `solver.solve()`).
4.  Plot the resulting spectrum.

## Comparison with Other Methods
| Method | Description | Key Characteristics |
| :--- | :--- | :--- |
| **ana_cont** | Library implementing MaxEnt and Pade | General purpose, standard MaxEnt/Pade algorithms |
| **SpM** | Sparse Modeling | Stable against noise, automatic basis selection |
| **SOM** | Stochastic Optimization Method | Sampling based, avoids default models, handles noise well |
| **Nevanlinna** | Nevanlinna analytical continuation | Preserves causality, potentially more rigorous |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/josefkaufmann/ana_cont
2. Documentation: https://josefkaufmann.github.io/ana_cont/
3. Dependencies: numpy, scipy, matplotlib.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Utility: Widely used Python tool for continuation
