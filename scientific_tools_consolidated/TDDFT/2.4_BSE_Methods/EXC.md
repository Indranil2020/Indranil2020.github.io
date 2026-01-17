# EXC

## Official Resources
- **Website**: http://www.bethe-salpeter.org/
- **Documentation**: http://www.bethe-salpeter.org/ (Manuals and Tutorials)
- **License**: Research/Academic License

## Overview
EXC is a dedicated _ab initio_ code for calculating the dielectric and optical properties of materials by solving the Bethe-Salpeter Equation (BSE). Developed at Ecole Polytechnique (LSI), it operates in reciprocal space and the frequency domain, utilizing a plane-wave basis. It is designed to capture electron-hole interaction effects (excitons) in absorption and energy loss spectra.

**Scientific domain**: Solid state physics, optical spectroscopy, core-level spectroscopy
**Target user community**: Researchers studying optical properties of solids, surfaces, and clusters

## Theoretical Methods
- **Bethe-Salpeter Equation (BSE)**: Two-particle Green's function formalism
- **GW Approximation**: Uses quasiparticle energies as input
- **Reciprocal Space**: Solving the BSE matrix in k-space
- **Haydock Method**: Iterative solver for finding spectra without full diagonalization

## Capabilities
- **Optical Absorption**: Dielectric function with excitonic effects
- **EELS**: Electron Energy Loss Spectroscopy
- **IXS/XAS**: Inelastic X-ray Scattering and X-ray Absorption Spectroscopy
- **Supercells**: Treatment of finite systems (atoms, clusters) via supercell approximation
- **Scissor Operator**: Support for simple GW corrections via energy shifts

## Inputs & Outputs
- **Input files**:
  - `input`: Main configuration file (parameters like `exciton`, `omegai`, `broad`)
  - `file.kss`: Ground state Kohn-Sham structure (density/wavefunctions)
  - `file.gw`: Quasiparticle energy corrections
  - `file.em1`: Inverse dielectric matrix (screening)
- **Output data types**:
  - Dielectric function (real and imaginary parts)
  - Absorption spectra
  - Reflectivity and Refractive Index data
  - Exciton binding energies

## Performance Characteristics
- **Parallelization**: MPI parallelization for matrix construction and diagonalization.
- **Solvers**: Haydock iterative method allows handling larger matrices than direct diagonalization.

## Comparisons
- **vs Yambo**: Both are plane-wave BSE codes. Yambo has a larger modern community and more integrated workflow. EXC is a robust, specialized solver with a long history.
- **vs OCEAN**: OCEAN focuses on core-level spectroscopy (BSE), EXC covers both valence and core (check specific XAS features).

## Usage & Best Practices
- **Workflow**: Ground State (ABINIT/others) -> GW (ABINIT/others) -> Screening -> EXC.
- **Convergence**: Careful convergence of k-points and empty bands is crucial for BSE.

## Citations
- **Primary**: Reference L. Reining et al. and the http://www.bethe-salpeter.org/ website.
