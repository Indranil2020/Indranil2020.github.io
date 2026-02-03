# libAtoms/Quippy

## Official Resources
- Homepage: https://libatoms.github.io/
- Documentation: https://libatoms.github.io/QUIP/
- Source Repository: https://github.com/libAtoms/QUIP
- License: GNU General Public License v2.0 (QUIP), MIT (GAP)

## Overview
libAtoms/Quippy is a software package for molecular dynamics and atomistic simulations, primarily known for its implementation of the Gaussian Approximation Potential (GAP) machine learning interatomic potential. The package consists of the QUIP core (written in Fortran) and Quippy (Python bindings), providing a powerful and flexible environment for running simulations with advanced potentials.

**Scientific domain**: Machine learning potentials, molecular dynamics, Gaussian Approximation Potentials (GAP)  
**Target user community**: Materials scientists, researchers using ML potentials

## Theoretical Methods
- Gaussian Approximation Potentials (GAP)
- Gaussian Process Regression (GPR)
- Smooth Overlap of Atomic Positions (SOAP) descriptors
- Classical Molecular Dynamics
- Geometry Optimization
- QM/MM embedding (LOTF - Learn on the Fly)
- Tight-binding models (via descriptors)

## Capabilities (CRITICAL)
- Training and using GAP potentials
- Running MD with machine learning potentials
- Calculating SOAP descriptors
- Integration with ASE (Atomic Simulation Environment)
- LOTF (Learn on the Fly) hybrid molecular dynamics
- Interfaces to LAMMPS and CP2K
- Flexible Python scripting via Quippy

**Sources**: libAtoms website, Phys. Rev. Lett. 104, 136403 (2010)

## Inputs & Outputs
- **Input formats**: XYZ (extended XYZ), classical potential parameter files (.xml)
- **Output data types**: Trajectories (XYZ), forces, energies, virials

## Interfaces & Ecosystem
- **ASE**: Deep integration as an ASE calculator
- **LAMMPS**: Can use QUIP potentials via `pair_style quip`
- **CP2K**: Integration for QM/MM
- **Python**: Full scripting via Quippy

## Workflow and Usage
1. Generate training data (DFT calculations)
2. Train GAP potential: `gap_fit` command
3. Run simulation: Use ASE with `QUIP` calculator or LAMMPS with `pair_style quip`
4. Analysis: Use Quippy/ASE tools

## Performance Characteristics
- GAP potentials are computationally more expensive than simple classical potentials (EAM/LJ) but much cheaper than DFT
- Parallelization via MPI (in QUIP and LAMMPS)
- Scaling depends on descriptor cutoffs and basis set size

## Application Areas
- Silicon crystallization and defects
- Carbon materials (amorphous carbon, graphene)
- Phase diagrams of metals (Iron, Tungsten)
- Radiation damage simulations
- High-accuracy property prediction

## Community and Support
- Open-source (GPL/MIT)
- GitHub repository
- Active development (Csányi, Kermode groups)
- Tutorial materials available

## Verification & Sources
**Primary sources**:
1. Homepage: https://libatoms.github.io/
2. GitHub: https://github.com/libAtoms/QUIP
3. Publication: A. P. Bartók et al., Phys. Rev. Lett. 104, 136403 (2010)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Cambridge/Warwick)
- Applications: GAP potentials, SOAP descriptors, MD with ML potentials
