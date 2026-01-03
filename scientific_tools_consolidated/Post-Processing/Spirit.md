# Spirit

## Official Resources
- Homepage: https://spirit-code.github.io/
- Documentation: https://spirit-code.github.io/
- Source Repository: https://github.com/spirit-code/spirit
- License: MIT License

## Overview
Spirit is a framework for spin dynamics and calculating energy landscapes of magnetic systems. It allows for atomistic spin dynamics (ASD) simulations and the determination of minimum energy paths and transition states (e.g., for magnetic skyrmion collapse) using the Geodesic Nudged Elastic Band (GNEB) method. Spirit provides a desktop GUI with real-time visualization and a powerful Python API.

**Scientific domain**: Atomistic spin dynamics, magnetism, transition states, skyrmions  
**Target user community**: Magnetism researchers, spintronics engineers

## Theoretical Methods
- Atomistic Spin Dynamics (ASD) based on Landau-Lifshitz-Gilbert (LLG)
- Geodesic Nudged Elastic Band (GNEB) method
- Heisenberg Hamiltonian (Exchange, DMI, Anisotropy, Zeeman, Dipole-Dipole)
- Monte Carlo simulations
- Minimum Mode Following
- Calculation of lifetimes via Harmonic Transition State Theory (HTST)

## Capabilities (CRITICAL)
- Real-time visualization of spin dynamics
- Finding energy minima and saddle points (transition states)
- Calculating minimum energy paths for magnetic transitions
- Simulation of skyrmions, domain walls, and vortices
- Temperature-dependent simulations (Langevin dynamics)
- Interactive GUI and scriptable Python interface
- Support for various lattices and defects

**Sources**: Spirit website, Phys. Rev. B 99, 224414 (2019)

## Inputs & Outputs
- **Input formats**: Config files (cfg), structural data, interaction parameters (J, D, K)
- **Output data types**: Spin configurations (OVF, text), energy paths, images, properties

## Interfaces & Ecosystem
- **Python**: `spirit` python package (ctypes interface)
- **GUI**: Qt-based desktop application
- **Web**: Web-based version available
- **VASP/QE**: Can use parameters derived from DFT (e.g., via TB2J)

## Workflow and Usage
1. Define system: Lattice, magnetic moments, interactions (J_ij, D_ij).
2. Initialize state: Random, Skyrmion, Domain Wall.
3. Run simulation: LLG dynamics or GNEB minimization.
4. Visualize: Use the GUI or export to ParaView/POV-Ray.

## Performance Characteristics
- High-performance C++ core
- Parallelized (OpenMP, CUDA for GPU)
- Efficient for large systems (millions of spins)

## Application Areas
- Skyrmion stability and lifetimes
- Domain wall motion
- Magnetic switching paths
- Frustrated magnetism
- Spintronic devices

## Community and Support
- Open-source (MIT)
- Developed by Forschungszentrum Jülich (G. Bihlmayer group)
- Active GitHub repository

## Verification & Sources
**Primary sources**:
1. Homepage: https://spirit-code.github.io/
2. GitHub: https://github.com/spirit-code/spirit
3. Publication: G. P. Müller et al., Phys. Rev. B 99, 224414 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (FZ Jülich)
- Applications: Spin dynamics, GNEB, skyrmions, visualization
