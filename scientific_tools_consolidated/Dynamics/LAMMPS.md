# LAMMPS

## Official Resources
- Homepage: https://www.lammps.org/
- Documentation: https://docs.lammps.org/
- Source Repository: https://github.com/lammps/lammps
- License: GNU General Public License v2.0

## Overview
LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is a classical molecular dynamics code with a focus on materials modeling. Developed at Sandia National Laboratories, it is the most widely used open-source MD code, featuring excellent parallel scaling, extensive force field support, and comprehensive capabilities for simulating atomic, polymeric, biological, solid-state, granular, and coarse-grained systems.

**Scientific domain**: Classical molecular dynamics, atomistic simulations, materials modeling  
**Target user community**: Materials scientists, chemists, physicists studying dynamics of atomic systems

## Theoretical Methods
- Classical molecular dynamics (MD)
- Energy minimization
- Canonical ensemble (NVT)
- Isothermal-isobaric ensemble (NPT)
- Microcanonical ensemble (NVE)
- Grand canonical Monte Carlo (GCMC)
- Steered molecular dynamics (SMD)
- Replica exchange molecular dynamics (REMD)
- Non-equilibrium MD (NEMD)
- Rigid body dynamics
- Coarse-grained models
- Dissipative particle dynamics (DPD)
- Smoothed particle hydrodynamics (SPH)
- Peridynamics

## Capabilities (CRITICAL)
- Classical MD simulations
- Energy minimization
- Various ensembles (NVE, NVT, NPT, etc.)
- Multiple time integration schemes
- Extensive force field library (CHARMM, AMBER, OPLS, ReaxFF, etc.)
- Reactive force fields (ReaxFF, COMB)
- Machine learning potentials (ACE, NequIP, DeePMD, etc.)
- Many-body potentials (EAM, MEAM, Tersoff, Stillinger-Weber)
- Coarse-grained models
- Rigid body dynamics
- Constraint algorithms (SHAKE, RATTLE)
- Long-range electrostatics (Ewald, PPPM, MSM)
- Enhanced sampling (metadynamics, umbrella sampling)
- Non-equilibrium simulations
- Shock physics
- Crack propagation
- Grain boundaries and interfaces
- Polymers and soft matter
- Biological systems
- Massively parallel (MPI+OpenMP)
- GPU acceleration (CUDA, OpenCL, Kokkos)
- Excellent scaling to millions of cores
- Python interface (PyLAMMPS)
- On-the-fly visualization

**Sources**: Official LAMMPS documentation (https://www.lammps.org/), confirmed in 7/7 source lists

## Key Strengths

### Performance:
- Excellent parallel scaling
- Multi-GPU support
- Millions of cores demonstrated
- Highly optimized
- Leadership-class supercomputer ready

### Force Fields:
- 50+ pair styles
- Reactive (ReaxFF, COMB)
- ML potentials integrated
- User-extensible
- Comprehensive library

### Versatility:
- Atoms to continuum
- Multiple material types
- Multiscale methods
- Diverse applications

### Community:
- Large user base
- Active development
- Extensive documentation
- Regular releases
- Open-source

### Extensibility:
- Package system
- User-contributed packages
- Python interface
- Library mode
- Easy to modify

## Inputs & Outputs
- **Input formats**:
  - LAMMPS input script (in.*)
  - Data file (structure, topology)
  - Various coordinate formats
  - Force field parameter files
  
- **Output data types**:
  - Dump files (trajectories)
  - Thermodynamic output (log file)
  - Restart files
  - Custom output via compute/fix
  - XYZ, DCD, NetCDF formats

## Interfaces & Ecosystem
- **Pre/post-processing**:
  - Pizza.py tools
  - OVITO for visualization
  - VMD, ParaView
  - ASE integration
  
- **Python interface**:
  - PyLAMMPS
  - Library mode
  - Jupyter notebooks
  - Automated workflows
  
- **Machine learning**:
  - DeePMD interface
  - PLUMED for enhanced sampling
  - ACE, NequIP, MACE potentials
  
- **Multiscale**:
  - QM/MM coupling
  - FEM coupling
  - Coarse-graining tools
  
- **Parallelization**:
  - MPI standard
  - OpenMP threads
  - GPU acceleration (CUDA, Kokkos)
  - Hybrid MPI+OpenMP+GPU

## Workflow and Usage

### Example Input Script:
```lammps
# Lennard-Jones fluid simulation

units           lj
atom_style      atomic

lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box

mass            1 1.0

pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5

velocity        all create 1.0 87287

fix             1 all nvt temp 1.0 1.0 0.1

thermo          100
dump            1 all custom 1000 dump.lammpstrj id type x y z

run             10000
```

### Common Workflow:
1. Prepare structure (data file)
2. Write input script
3. Run: `lmp_serial < in.file`
4. Analyze trajectory
5. Compute properties

## Advanced Features

### Reactive Force Fields:
- ReaxFF for bond breaking/formation
- COMB for charge transfer
- Dynamic bonding
- Chemical reactions

### Machine Learning Potentials:
- ACE (Atomic Cluster Expansion)
- NequIP (neural equivariant interatomic potentials)
- DeePMD integration
- SNAP (spectral neighbor analysis)
- GAP interface

### Enhanced Sampling:
- Metadynamics (via PLUMED)
- Umbrella sampling
- Temperature replica exchange
- Steered MD
- Bias potentials

### Multiscale Methods:
- Coarse-graining
- Atomistic to continuum
- Coupled simulations
- Hierarchical modeling

### Non-Equilibrium:
- NEMD for transport properties
- Shock simulations
- Deformation
- Flow

### Specialized Physics:
- Granular materials
- Peridynamics (fracture)
- SPH (fluids)
- DPD (soft matter)
- Rigid body dynamics

## Performance Characteristics
- **Scaling**: Excellent weak scaling to millions of cores
- **GPU**: 5-50x speedup depending on system
- **Efficiency**: >80% parallel efficiency at scale
- **Typical systems**: 100-10,000,000 atoms
- **Timestep**: 0.5-2 fs typical

## Computational Cost
- **Force evaluation**: Dominant cost
- **Long-range**: Ewald/PPPM expensive
- **ReaxFF**: Much more expensive than simple potentials
- **ML potentials**: Moderate overhead
- **Scaling**: Near-linear for large systems

## Limitations & Known Constraints
- **Classical mechanics**: No electronic structure
- **Force fields**: Quality depends on parameterization
- **Learning curve**: Steep; many options
- **Input scripts**: Can be complex
- **Documentation**: Extensive but overwhelming
- **GPU**: Not all features supported
- **Platform**: Linux primarily, Windows/macOS possible

## Comparison with Other MD Codes
- **vs GROMACS**: LAMMPS more versatile, GROMACS faster for biomolecules
- **vs NAMD**: LAMMPS more general, NAMD specialized for biomolecules
- **vs AMBER**: LAMMPS open-source, more force fields
- **vs DL_POLY**: LAMMPS better scaling and features
- **Unique strength**: Versatility, performance, extensibility, community

## Application Areas

### Materials Science:
- Mechanical properties
- Fracture and failure
- Phase transitions
- Grain boundaries
- Defects

### Soft Matter:
- Polymers
- Colloids
- Emulsions
- Liquid crystals

### Biomolecules:
- Proteins and DNA
- Membranes
- Drug delivery
- Solvation

### Solid-State:
- Crystals
- Interfaces
- Thermal transport
- Mechanical behavior

### Geophysics:
- Minerals
- Friction
- Seismology applications

## Best Practices

### System Preparation:
- Use appropriate force field
- Check initial structure
- Minimize energy first
- Equilibrate carefully

### Simulation Setup:
- Choose appropriate ensemble
- Set timestep conservatively (1 fs typical)
- Use constraints if needed
- Monitor conserved quantities

### Equilibration:
- Gradual heating/cooling
- Check energy drift
- Sufficient equilibration time
- Monitor target properties

### Production:
- Long enough for statistics
- Frequent trajectory dumps
- Compute properties on-the-fly
- Save restart files

### Performance:
- Domain decomposition for parallelization
- Load balancing
- GPU acceleration when available
- Optimize communication

## Community and Support
- Open-source on GitHub
- Large user community
- Active mailing list
- Forum for questions
- Workshops and tutorials
- Regular releases
- Extensive manual

## Educational Resources
- Comprehensive manual
- Tutorial examples
- Workshop materials
- Video tutorials
- Published books
- Active community Q&A

## High-Performance Computing
- Optimized for DOE supercomputers
- Exascale-ready
- GPU acceleration mature
- Kokkos performance portability
- Production use on leadership facilities

## Verification & Sources
**Primary sources**:
1. Official website: https://www.lammps.org/
2. Documentation: https://docs.lammps.org/
3. GitHub repository: https://github.com/lammps/lammps
4. S. Plimpton, J. Comp. Phys. 117, 1 (1995) - Original LAMMPS paper
5. A. P. Thompson et al., Comp. Phys. Comm. 271, 108171 (2022) - LAMMPS update

**Secondary sources**:
1. LAMMPS manual and tutorials
2. Published MD studies using LAMMPS (>20,000 citations)
3. Workshop materials
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in ALL 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v2)
- Community support: Very active mailing list, forum
- Academic citations: >25,000
- Active development: Weekly updates, regular releases
- Benchmark validation: Extensively validated, industry standard
- HPC optimization: Production use on exascale systems
- Industry standard: Most widely used open-source MD code
