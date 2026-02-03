# Metadynamics

## Official Resources
- Homepage: Method - Implemented in PLUMED, CP2K, VASP, etc.
- Documentation: https://www.plumed.org/doc (PLUMED implementation)
- Source Repository: https://github.com/plumed/plumed2 (Reference implementation)
- License: Varies by implementation (LGPL for PLUMED)

## Overview
Metadynamics is an enhanced sampling method used to reconstruct free energy surfaces and accelerate rare events in molecular dynamics simulations. It works by adding a history-dependent bias potential (usually Gaussian hills) to selected collective variables (CVs), encouraging the system to explore new regions of phase space.

**Scientific domain**: Enhanced sampling, free energy calculations, rare events
**Target user community**: Computational chemists, materials scientists, physicists

## Theoretical Methods
- Standard Metadynamics
- Well-Tempered Metadynamics (WT-MetaD)
- Bias-Exchange Metadynamics
- Parallel Tempering Metadynamics
- Transition Path Sampling
- Reweighting techniques

## Capabilities (CRITICAL)
- Reconstructing free energy landscapes (FES)
- Accelerating transitions between metastable states
- Sampling high-dimensional reaction coordinates
- Estimating equilibrium properties from non-equilibrium runs
- Implemented in: PLUMED (plugin for many codes), CP2K, DESMO-J, and native implementations in some MD codes

**Sources**: Laio and Parrinello, PNAS 99, 12562 (2002)

## Inputs & Outputs
- **Input**: Definition of collective variables (CVs), Gaussian height/width, deposition rate
- **Output**: Time series of CVs, HILLS file (list of added biases), Reconstructed Free Energy Surface

## Interfaces & Ecosystem
- **PLUMED**: The de facto standard library for metadynamics
- **VASP**: Native or PLUMED interface
- **CP2K**: Strong internal implementation
- **GROMACS/LAMMPS/NAMD**: Via PLUMED

## Workflow and Usage
1. **Identify CVs**: Choose variables that describe the slow degrees of freedom.
2. **Setup**: Configure metadynamics parameters (bias factor, sigma, height).
3. **Run**: Execute MD simulation with metadynamics bias active.
4. **Analysis**: Sum the deposited Gaussians to obtain the negative free energy surface.
5. **Convergence**: Check for diffusive behavior in CV space.

## Performance Characteristics
- **Overhead**: Requires calculation of CVs and bias update (usually low overhead)
- **Convergence**: Depends critically on the choice of CVs
- **Scalability**: Parallelizable (multiple walkers)

## Application Areas
- Crystal nucleation
- Chemical reactions
- Protein folding/unfolding
- Ligand binding
- Phase transitions

## Community and Support
- Large community centered around PLUMED
- Active development of new variants (e.g., OPES)
- Theoretical chemistry community

## Verification & Sources
**Primary sources**:
1. PLUMED: https://www.plumed.org/
2. Publication: A. Laio and M. Parrinello, PNAS 99, 12562 (2002)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Method: STANDARD
- Documentation: COMPREHENSIVE (via PLUMED)
- Applications: Free energy, enhanced sampling
