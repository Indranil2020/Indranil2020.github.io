# Wannier.jl

## Official Resources
- **Repository**: https://github.com/qilauea/Wannier.jl
- **Documentation**: https://qilauea.github.io/Wannier.jl/dev/
- **License**: MIT License

## Overview
**Wannier.jl** is a modern, pure Julia package for generating **Maximally Localized Wannier Functions (MLWFs)**. Designed as a flexible and high-performance alternative to the standard Fortran **Wannier90** code, it implements core algorithms for disentanglement and localization while leveraging Julia's strengths in composability, automatic differentiation, and GPU acceleration. It allows for fast prototyping of new Wannierization methods and seamless integration with Julia-based density functional theory codes like **DFTK.jl**.

**Scientific domain**: Electronic Structure, Quantum Chemistry, Condensed Matter Physics
**Target user community**: Julia developers, researchers developing new Wannier algorithms, and users of DFTK.jl

## Theoretical Methods
- **Maximally Localized Wannier Functions (MLWF)**: Minimization of the spread functional $\Omega$.
- **Disentanglement**:
  - Implements the Souza-Marzari-Vanderbilt (SMV) method for entangled bands (metals).
  - Supports frozen and outer windows for subspace selection.
- **Manifold Optimization**: Uses Riemannian optimization techniques for unitary updates on the Stiefel manifold (optional algorithms).
- **Interpolation**: Fourier interpolation of operators (Hamiltonian, Position) to arbitrary k-points.
- **Auto-Differentiation**: Architecture compatible with Julia's AD ecosystem for gradient-based analysis (experimental).

## Capabilities
- **Wannierization**:
  - Calculation of $U(k)$ matrices to gauge-transform Bloch states.
  - Calculation of centers and spreads of Wannier functions.
- **Post-Processing**:
  - Band structure interpolation.
  - Fermi surface plotting.
  - Berry curvature and anomalous Hall conductivity (via interpolation).
- **Hardware Acceleration**:
  - Experimental support for GPU execution via **CUDA.jl**.
  - Multi-threading support built into Julia.

## Key Strengths
- **Language**: Pure Julia allows for easier inspection, modification, and extension than legacy Fortran.
- **Integration**: Direct in-memory coupling with **DFTK.jl** avoids disk I/O bottlenecks.
- **Algorithmic Experiments**: A playground for testing new minimization schemes (e.g., AD-based descent).
- **Reproducibility**: Strictly type-checked and tested against Wannier90 benchmarks.

## Inputs & Outputs
- **Inputs**:
  - Standard Wannier90 files (`.mmn`, `.amn`, `.eig`) via **WannierIO.jl**.
  - Direct data structures from Julia SCF codes.
- **Outputs**:
  - `.chk` checkpoint files compatible with Wannier90.
  - Interpolated physical quantities (Bands, Berry curvature).
  - Export to `_hr.dat` for other post-processors.

## Interfaces & Ecosystem
- **Upstream**:
  - **DFTK.jl**: Primary DFT engine in the Julia ecosystem.
  - **Wannier90**: Full file compatibility allows input from VASP, QE, etc.
- **I/O Backend**: Relies on **WannierIO.jl**.

## Performance Characteristics
- **Speed**: Comparable to the Fortran `wannier90.x` for serial execution on typical system sizes.
- **Accelerator Support**: Unlocks GPU potential which is currently limited in the standard Fortran version.
- **Scalability**: Currently limited to single-node shared memory (no MPI support yet), suitable for small to medium-sized systems (hundreds of atoms).

## Limitations & Known Constraints
- **MPI**: Lack of MPI support limits its use for extremely large-scale calculations compared to the MPI-parallelized `wannier90.x`.
- **Feature Parity**: While core features are present, some specialized modules of Wannier90 (e.g., transport, symmetry-adapted Wannier functions) may be less mature or missing.

## Comparison with Other Codes
- **vs. Wannier90**: Wannier.jl is more modular and amenable to algorithmic research/AD, while Wannier90 is the battle-tested, MPI-parallel industry standard for production runs on supercomputers.
- **vs. WannierTools**: WannierTools is a specific post-processor (Topological analysis); Wannier.jl generates the WFs themselves (like Wannier90).

## Application Areas
- **Method Development**: Testing new spread functionals or optimization strategies.
- **High-Throughput**: Scriptable workflows in Julia.
- **Education**: Easier to read and understand the mathematics of Wannierization.

## Community and Support
- **Development**: Maintained by the Qilauea organization (collaborators from various institutions).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/qilauea/Wannier.jl](https://github.com/qilauea/Wannier.jl)
- **Primary Website**: [https://wannierjl.org/](https://wannierjl.org/)
- **Verification status**: ✅ VERIFIED
  - Active and growing project.
  - Validated against exact results for standard benchmarks (Silicon, Copper, etc.).
