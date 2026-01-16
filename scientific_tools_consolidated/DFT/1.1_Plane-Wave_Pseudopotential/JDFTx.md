# JDFTx

## Official Resources
- Homepage: https://jdftx.org/
- Documentation: https://jdftx.org/Documentation.html
- Source Repository: https://github.com/shankar1729/jdftx
- License: GNU General Public License v3.0

## Overview
JDFTx (Joint Density Functional Theory - extended) is a plane-wave DFT code designed for electronic structure calculations in molecules, surfaces, and liquids, with particular emphasis on solvation and electrochemistry. Developed by Ravishankar Sundararaman and collaborators, JDFTx features efficient joint density functional theory for implicit solvation, advanced models for electrochemical interfaces, and GPU acceleration. It is particularly strong in studying charged systems, solvated systems, and electrochemical processes.

**Scientific domain**: DFT, solvation, electrochemistry, charged systems, surfaces  
**Target user community**: Electrochemists, surface scientists, battery researchers, solvation specialists

## Theoretical Methods
- Kohn-Sham DFT (LDA, GGA, meta-GGA)
- Plane-wave basis with pseudopotentials
- Norm-conserving and ultrasoft pseudopotentials
- Joint Density Functional Theory (JDFT)
- Implicit solvation models (PCM, CANDLE, GLSSA13)
- Explicit fluid models
- Electrochemical interfaces
- Charged systems with background charge
- Grand canonical DFT
- van der Waals corrections (DFT-D2/D3)
- DFT+U for correlated systems
- Hybrid functionals
- Non-local correlation functionals
- Spin-orbit coupling
- Explicit solvent models
- Grand canonical DFT for charged systems

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Joint density-functional theory for solvation
- Implicit solvation models (PCM, CANDLE, SaLSA)
- Explicit solvent via classical density functionals
- Charged systems and grand canonical ensemble
- Electrochemical double layers
- Geometry optimization in solution
- Vibrational calculations with solvation
- Phonon calculations
- Band structure and DOS
- Wannier functions
- Berry phase polarization
- GPU acceleration for key operations
- Efficient handling of charged periodic systems
- Lattice optimization under constant potential

**Sources**: Official JDFTx documentation, cited in 5/7 source lists

## Inputs & Outputs
- **Input formats**:
  - Text input files (JDFTx format)
  - Coordinate files (XYZ, POSCAR-like)
  - Pseudopotential files (UPF, USPP)
  - Ionosphere files for electrolyte models
  
- **Output data types**:
  - Standard output with energies, forces, stresses
  - Electronic density (binary format)
  - Wavefunction files
  - DOS and band structure files
  - Solvation structure output
  - Wannier function output

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE - calculator interface available
  - Can interface with standard DFT workflows
  
- **GPU support**:
  - CUDA acceleration for FFTs and operations
  - Significant speedup on GPUs vs CPUs
  
- **Solvation models**:
  - NonlinearPCM - polarizable continuum
  - CANDLE - cavity and dispersion
  - SaLSA - nonlocal solvation
  - Classical DFT for explicit fluids

## Limitations & Known Constraints
- **Specialized focus**: Optimized for solvation; less feature-complete than general codes
- **Documentation**: Good but less extensive than major codes
- **Community**: Smaller user base than VASP/QE
- **Pseudopotentials**: Must use compatible pseudopotential formats
- **Solvation convergence**: Can be challenging for complex systems
- **GPU requirement**: Best performance requires CUDA-capable GPU
- **Input format**: Unique syntax requires learning curve
- **Post-processing**: Fewer established analysis tools
- **Platform support**: Primarily Linux/Unix; GPU drivers required

## Computational Cost
- **Scaling**: Excellent on GPUs; CPU scaling is standard $O(N^3)$.
- **Solvation**: Implicit solvation adds minimal overhead (~10-20%) compared to vacuum calculations, much cheaper than explicit solvent.
- **Memory**: GPU memory can be a bottleneck for large systems.

## Comparison with Other Codes
- **vs VASP/QE**: JDFTx is specialized for solvation/electrochemistry. If you need standard solid-state properties, VASP/QE are more feature-rich. If you need electrochemical interfaces, JDFTx is superior.
- **vs Qbox**: Both are C++ and scalable, but JDFTx focuses on solvation/embedding, Qbox on massive MD.

## Best Practices
- **GPU**: Use `jdw` (JDFTx wrapper) to manage GPU resources effectively.
- **Solvation**: Start with `PCM` or `CANDLE` before moving to more complex `SaLSA` models.
- **Minimization**: Use `Minimize` for electronic steps; it is often more robust than `SCF` mixing for difficult solvated systems.

## Community and Support
- **Support**: GitHub Issues are the primary support channel.
- **Development**: Active interaction with developers (Sundararaman group).

## Verification & Sources
**Primary sources**:
1. Official website: https://jdftx.org/
2. Documentation: https://jdftx.org/Documentation.html
3. GitHub repository: https://github.com/shankar1729/jdftx
4. R. Sundararaman et al., J. Chem. Phys. 146, 114104 (2017) - JDFTx code
5. K. A. Schwarz et al., Phys. Rev. B 85, 201102 (2012) - JDFT formalism

**Secondary sources**:
1. JDFTx tutorials and examples
2. Published electrochemistry applications
3. Solvation model benchmarks
4. Confirmed in 5/7 source lists (claude, g, gr, k, q)

**Confidence**: CONFIRMED - Appears in 5 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub)
- Community support: Active (GitHub issues, email)
- Academic citations: >100 (main papers)
- Active development: Regular updates, GPU optimization
