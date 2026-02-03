# LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator)

## Official Resources
- Homepage: https://www.lammps.org/
- Documentation: https://docs.lammps.org/
- Source Repository: https://github.com/lammps/lammps
- License: GNU General Public License v2.0

## Overview
LAMMPS is a classical molecular dynamics code with a focus on materials modeling. It's an acronym for Large-scale Atomic/Molecular Massively Parallel Simulator. LAMMPS has potentials for solid-state materials (metals, semiconductors) and soft matter (biomolecules, polymers) and coarse-grained or mesoscopic systems. It can be used to model atoms or, more generically, as a parallel particle simulator at the atomic, meso, or continuum scale.

**Scientific domain**: Molecular dynamics, materials simulation, soft matter, coarse-grained modeling  
**Target user community**: Materials scientists, physicists, chemists, engineers

## Theoretical Methods
- Classical molecular dynamics (Newton's equations)
- Brownian and Langevin dynamics
- Energy minimization
- Nonequilibrium molecular dynamics (NEMD)
- Grand Canonical Monte Carlo (GCMC)
- Dissipative Particle Dynamics (DPD)
- Peridynamics
- Smooth Particle Hydrodynamics (SPH)
- Time-Dependent Ginzburg-Landau (TDGL)
- Reaction dynamics (ReaxFF)
- Spin dynamics
- Electron force field (eFF)

## Capabilities (CRITICAL)
- Atomic, polymeric, biological, solid-state, granular, coarse-grained simulations
- Massive parallelization (MPI, OpenMP, GPU, Kokkos)
- Huge variety of interatomic potentials (Lennard-Jones, EAM, MEAM, Tersoff, ReaxFF, AIREBO, COMB, SW, etc.)
- Advanced constraints and boundary conditions
- On-the-fly analysis and post-processing
- Coupling with other codes (quantum, FE)
- Python interface
- User-extendable via C++ classes

**Sources**: LAMMPS documentation, Comp. Phys. Comm. 183, 1136 (2012)

## Key Strengths

### Versatility:
- Huge variety of potentials (100+)
- Materials and soft matter
- Coarse-grained to atomistic
- User-extensible via C++

### Parallelization:
- Excellent MPI scaling
- GPU acceleration (Kokkos, GPU package)
- Billions of atoms possible
- Load balancing

### Ecosystem:
- Python interface
- ASE integration
- PLUMED support
- Extensive community packages

## Inputs & Outputs
- **Input formats**: Text-based script files, data files for initial structure
- **Output data types**: Dump files (text/binary/custom), log files (thermodynamics), restart files, XTC, DCD

## Interfaces & Ecosystem
- **Python**: Full Python wrapper
- **ASE**: Integrated via ASE calculator
- **OVITO**: Visualization standard
- **Moltemplate/Topotools**: Input generation
- **PLUMED**: Enhanced sampling
- **Phonopy**: Phonon calculations
- **VMD**: Visualization

## Workflow and Usage
1. Prepare structure (data file)
2. Write input script (units, potential, fixes, run)
3. Run: `lmp_mpi -in input.in`
4. Visualize output (dump file)

## Performance Characteristics
- Highly scalable (millions to billions of particles)
- Optimized for MPI and accelerators (GPU/Kokkos)
- Load balancing for inhomogeneous systems

## Computational Cost
- Scales linearly with atoms for short-range
- Long-range (Ewald/PPPM) adds overhead
- GPU provides 10-100x speedup
- Overall: Highly efficient for materials

## Best Practices
- Use appropriate units for your system
- Choose neighbor list settings carefully
- Validate potential for your application
- Use restart files for long runs
- Check energy conservation in NVE

## Limitations & Known Constraints
- Less optimized for biomolecules than GROMACS
- Steeper learning curve than some codes
- Some packages require compilation
- Documentation can be overwhelming

## Application Areas
- Metals and alloys (defects, mechanics)
- Polymers and biomolecules
- Nanostructures and interfaces
- Granular materials
- Shock physics
- Thermal transport
- Chemical reactions (ReaxFF)

## Comparison with Other Codes
- **vs GROMACS**: LAMMPS more materials-focused, GROMACS faster for biomolecules
- **vs AMBER**: LAMMPS open-source with more potentials, AMBER better force fields for bio
- **vs OpenMM**: LAMMPS more general, OpenMM more flexible custom forces
- **Unique strength**: Unmatched potential variety, materials science focus, extensibility

## Community and Support
- Open-source (GPL v2)
- Very active mailing list
- Large user community
- Frequent stable releases
- Extensive workshops and tutorials

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.lammps.org/
2. Documentation: https://docs.lammps.org/
3. GitHub: https://github.com/lammps/lammps
4. Publication: S. Plimpton, J. Comp. Phys. 117, 1 (1995)

**Secondary sources**:
1. LAMMPS tutorials and workshops
2. Extensive published applications
3. Community packages documentation

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (Sandia National Labs, Temple University)
- Applications: MD, materials, soft matter, parallel computing, huge potential library
