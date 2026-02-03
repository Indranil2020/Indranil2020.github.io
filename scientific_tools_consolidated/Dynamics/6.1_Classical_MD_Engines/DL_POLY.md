# DL_POLY

## Official Resources
- Homepage: https://www.scd.stfc.ac.uk/software/dl_poly/
- Documentation: https://www.scd.stfc.ac.uk/software/dl_poly/documentation
- Source Repository: https://gitlab.com/DL_POLY_Classic/dl_poly_classic (Classic version) / DL_POLY_4 (Licensed)
- License: Proprietary (DL_POLY_4) / BSD (DL_POLY_Classic)

## Overview
DL_POLY is a general-purpose classical molecular dynamics simulation package developed at Daresbury Laboratory. It is designed to run on a wide range of computers, from single processor workstations to massively parallel supercomputers. DL_POLY handles a very wide variety of molecular systems including macromolecules, polymers, ionic systems, solutions, and surfaces.

**Scientific domain**: Classical molecular dynamics, materials science, chemistry  
**Target user community**: Academic and industrial researchers in materials and chemistry

## Theoretical Methods
- Classical Molecular Dynamics
- Domain Decomposition (DL_POLY_4)
- Replicated Data (DL_POLY_Classic)
- Rigid Body Dynamics
- Shell Model for Polarization
- Multiple Timestep Algorithms
- Free Energy Methods (Thermodynamic Integration)
- Metadynamics

## Capabilities (CRITICAL)
- Simulation of huge systems (millions of atoms with DL_POLY_4)
- Extensive range of force fields and potentials
- Ionic materials (oxides, minerals)
- Biological systems (proteins, DNA)
- Polymers and macromolecules
- Metals and alloys (EAM, Sutton-Chen)
- Non-equilibrium MD (shear, thermal gradients)
- Parallel performance (MPI)

**Sources**: STFC website, Mol. Simulat. 28, 95 (2002)

## Key Strengths

### Materials Focus:
- Ionic materials (oxides, minerals)
- Metals (EAM, Sutton-Chen)
- Shell model polarization
- Radiation damage

### Parallelization:
- Domain decomposition (DL_POLY_4)
- Excellent HPC scaling
- Millions of atoms

### Versatility:
- Wide range of potentials
- Multiple ensembles
- Non-equilibrium MD

## Inputs & Outputs
- **Input formats**: CONTROL (simulation parameters), CONFIG (coordinates), FIELD (force field), TABLE (tabulated potentials)
- **Output data types**: HISTORY (trajectory), OUTPUT (log), REVCON (restart), STATIS (statistics)

## Interfaces & Ecosystem
- **GUI**: Java-based GUI available
- **Analysis**: DL_FIELD, DL_ANALYSER
- **Python**: Analysis scripts
- **VMD**: Visualization support

## Workflow and Usage
1. Prepare system: Generate CONFIG and FIELD files (using DL_FIELD or other tools)
2. Define control: Create CONTROL file
3. Run: `DLPOLY.X`
4. Analysis: Process HISTORY and STATIS files

## Performance Characteristics
- DL_POLY_4: Excellent scaling on massively parallel systems (Domain Decomposition)
- DL_POLY_Classic: Good for smaller systems (Replicated Data)
- Optimized for HPC environments

## Computational Cost
- Excellent parallel scaling (DL_POLY_4)
- Efficient for ionic systems
- Good for large systems
- Overall: HPC-optimized

## Best Practices
- Use DL_POLY_4 for large systems
- Validate potential parameters
- Use appropriate cutoffs
- Check energy conservation

## Limitations & Known Constraints
- DL_POLY_4 requires license
- Less biomolecular focus
- Older interface style
- Limited GPU support

## Application Areas
- Solid state materials (defects, diffusion)
- Ionic liquids and molten salts
- Biomolecular simulations
- Surface science and catalysis
- Radiation damage (cascades)

## Comparison with Other Codes
- **vs LAMMPS**: DL_POLY better ionic/shell model, LAMMPS more potentials
- **vs GROMACS**: DL_POLY materials focus, GROMACS biomolecular
- **Unique strength**: Shell model polarization, ionic materials, radiation damage

## Community and Support
- Developed by STFC Daresbury Laboratory
- User workshops and training
- Mailing list
- Classic version is open source

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.scd.stfc.ac.uk/software/dl_poly/
2. GitLab: https://gitlab.com/DL_POLY_Classic/dl_poly_classic
3. Publication: I.T. Todorov et al., J. Mater. Chem. 16, 1911 (2006)

**Secondary sources**:
1. DL_POLY documentation
2. STFC training materials
3. Published materials science applications

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: MIXED (Classic open, DL_POLY_4 licensed)
- Development: ACTIVE (STFC)
- Applications: MD, materials, ionic systems, parallel computing
