# NEMO5

## Official Resources
- **Homepage**: https://nemo5.org/
- **Distribution**: https://nanohub.org/groups/nemo5_distribution
- **License**: Dual (Academic / Commercial via Silvaco)

## Overview
**NEMO5** (NanoElectronics MOdeling Tools 5) is a state-of-the-art **multiscale simulation framework** for nanoelectronics. It represents the culmination of decades of development (NEMO-1D, NEMO-3D, OMEN) and integrates atomistic tight-binding models with continuum approximations, strain engineering, phonon transport, and optical properties. It is designed to simulate realistic semiconductor devices (FinFETs, nanowires, quantum dots) with millions of atoms, leveraging petascale supercomputing resources.

**Scientific domain**: Nanoelectronics, TCAD, Semiconductor Physics
**Target user community**: Device physicists and advanced TCAD engineers

## Theoretical Methods
- **Transport**:
  - Non-Equilibrium Green's Function (NEGF).
  - Quantum Transmitting Boundary Method (QTBM) / Wave Function formalism.
  - Recursive Green's Function (RGF) for linear scaling in length.
- **Electronic Structure**:
  - Tight-Binding ($sp^3d^5s^*$, etc.).
  - Effective Mass and $k \cdot p$.
- **Lattice Dynamics**:
  - Valence Force Field (VFF) for strain relaxation.
  - Phonon band structure and transport.
- **Electrostatics**:
  - 1D/2D/3D Schrödinger-Poisson self-consistent solver.

## Capabilities
- **Device Simulation**:
  - CMOS scaling (FinFET, GAA-FET).
  - Tunneling FETs (TFETs).
  - Solar cells and LEDs.
  - Quantum dots for qubits.
- **Multiphysics**:
  - Impact of strain on mobility.
  - Phonon-limited mobility/scattering.
  - Optical absorption/emission spectra.
- **Scale**: From small atomistic clusters to realistic devices with >10 million atoms.

## Key Strengths
- **Massive Scalability**: optimized for hybrid MPI/OpenMP parallelism, scaling to hundreds of thousands of cores.
- **Versatility**: A unified tool that can switch between levels of theory (e.g., Effective Mass vs Tight Binding) within the same workflow.
- **Atomistic Strain**: Best-in-class handling of lattice relaxation in heterostructures.

## Inputs & Outputs
- **Inputs**: Tcl or Python scripts defining the geometry, material parameters, and physics solvers.
- **Outputs**:
  - VTK/Silvaco formats for 3D fields.
  - Band structures and I-V curves.

## Interfaces & Ecosystem
- **nanoHUB**: Available as a cloud tool.
- **Silvaco**: Commercialized as "Victory Atomistic".
- **Scripting**: Complex logic controlled via scripting languages.

## Performance Characteristics
- **Computational Cost**: High. Full NEGF with scattering is extremely expensive.
- **Efficiency**: Uses advanced RGF algorithms and multilevel parallelization.

## Comparison with Other Codes
- **vs. OMEN**: NEMO5 includes OMEN's transport kernels but adds broader multiphysics (optics, strain, continuum models).
- **vs. QuantumATK**: Similar scope; NEMO5 is academic/US-based, ATK is commercial (Synopsys).
- **vs. Kwant**: Kwant is a library for model Hamiltonians; NEMO5 is a full device simulator with material databases and Poisson solvers.

## Application Areas
- **End-of-Roadmap CMOS**: Predicting performance of 3nm/2nm node devices.
- **Quantum Computing**: Modeling spin qubits in Si/Ge heterostructures.
- **Alloy Disorder**: Statistical variability in random alloys.

## Community and Support
- **Development**: Purdue University (Gerhard Klimeck group) and NCN.
- **Source**: Available to academic groups; commercial support via Silvaco.

## Verification & Sources
- **Website**: [https://nemo5.org/](https://nemo5.org/)
- **Primary Publication**: S. Steiger et al., IEEE Trans. Nanotechnol. 10, 1464 (2011).
- **Verification status**: ✅ VERIFIED
  - Benchmark code in the field of atomistic TCAD.
