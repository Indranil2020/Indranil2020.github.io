# DFTpy

## Official Resources
- Homepage: https://gitlab.com/pavanello-research-group/dftpy
- Documentation: https://dftpy.readthedocs.io/
- Source Repository: https://gitlab.com/pavanello-research-group/dftpy
- License: MIT License

## Overview
DFTpy is a modern, pure Python framework for Orbital-Free Density Functional Theory (OF-DFT) and Kohn-Sham DFT (KS-DFT). It utilizes a plane-wave basis set represented on a real-space grid via Fast Fourier Transforms (FFTs). Developed at Rutgers University (Pavanello Research Group), it is designed to be highly modular, enabling the rapid development and testing of new functionals (especially Kinetic Energy Functionals), embedding potentials, and non-standard SCF drivers.

**Scientific domain**: Orbital-Free DFT, Embedding Theory, Method Development
**Target user community**: Developers of new functionals, researchers in multi-scale embedding

## Theoretical Methods
- **Orbital-Free DFT (OF-DFT)**: Scales linearly $O(N)$ by evolving the density directly.
- **Kohn-Sham DFT (KS-DFT)**: Standard orbital-based approach available.
- **Basis Set**: Plane-waves (reciprocal space) / Uniform Grid (real space).
- **Potentials**: Local and Non-local Pseudopotentials (BLPS, PAW support in development).
- **Time-Dependent DFT**: Hydrodynamic formulations.

## Capabilities
- Ground-state energy minimization (Direct Minimization).
- Geometry optimization.
- Ab initio Molecular Dynamics (AIMD).
- Density-based embedding (subsystem DFT).
- Parallelization via MPI (domain decomposition).

## Key Strengths

### Innovation:
- One of the few modern codes focusing on **Orbital-Free DFT**, allowing simulation of very large systems (thousands of atoms) at linear scaling cost.
- **Embedding**: Excellent tools for density embedding protocols.

### Usability:
- **Python-First**: Entirely written in Python; easy to interface with ML libraries (PyTorch) or other Python tools.

## Inputs & Outputs
- **Input formats**:
  - Python scripts (standard mode).
  - Structure files (XYZ, POSCAR, CIF via ASE).
  - Pseudopotentials (UPF, standard formats).
  
- **Output data types**:
  - Density files (Cube, XSF).
  - Energy/Force logs.
  - Trajectory files.

## Interfaces & Ecosystem
- **ASE**: Full integration with the Atomic Simulation Environment.
- **Libxc**: Interfaces for exchange-correlation functionals.
- **PyTorch**: Experimental support for ML-based functionals.

## Computational Cost
- **Orbital-Free Mode**: extremely efficient ($O(N)$), capable of handling 10,000+ atoms on modest hardware if the KEDF is accurate enough for the system.
- **Kohn-Sham Mode**: Slower than optimized C++/Fortran codes (QE/VASP) due to Python overhead, though heavy FFTs are offloaded to efficient backends (FFTW via pyfftw).
- **Scaling**: Good MPI scaling for the grid operations.

## Best Practices
- **OF-DFT Usage**: Use for liquid metals or systems where Kinetic Energy Functionals are known to work well; avoid for covalent bonds unless using advanced non-local KEDFs.
- **Backend**: Ensure `pyfftw` or a GPU-accelerated FFT backend is configured for production performance.

## Comparison with Other Codes
- **vs eminus**: Both are Python-based, but DFTpy focuses more on Orbital-Free DFT and embedding, while eminus is strictly a Kohn-Sham educational prototype.
- **vs PROFESS**: PROFESS is a dedicated C++ OF-DFT code (faster); DFTpy is more flexible for development but slower in pure execution.

## Community and Support
- **Hosting**: GitLab.
- **Support**: Issue tracker and direct academic contact.
- **Status**: Active research code.

## Verification & Sources
**Primary sources**:
1. Official GitLab Repository
2. "DFTpy: An efficient and scalable software package..." (Journal of Chemical Physics)
3. Documentation (ReadTheDocs)

**Confidence**: VERIFIED - Active research software with peer-reviewed publications.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: OF-DFT/Python
- Key Feature: Linear Scaling (OF-DFT)
