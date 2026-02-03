# GROMACS (GROningen MAchine for Chemical Simulations)

## Official Resources
- Homepage: https://www.gromacs.org/
- Documentation: https://manual.gromacs.org/
- Source Repository: https://gitlab.com/gromacs/gromacs
- License: GNU Lesser General Public License v2.1

## Overview
GROMACS is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for systems with hundreds to millions of particles. It is primarily designed for biochemical molecules like proteins, lipids, and nucleic acids that have a lot of complicated bonded interactions, but since GROMACS is extremely fast at calculating the nonbonded interactions (that usually dominate simulations), many groups are also using it for research on non-biological systems, e.g. polymers.

**Scientific domain**: Biomolecular simulation, molecular dynamics, soft matter  
**Target user community**: Biochemists, biophysicists, polymer physicists

## Theoretical Methods
- Classical molecular dynamics (Newton's/Langevin equations)
- Energy minimization (Steepest descent, Conjugate gradient)
- Free energy calculations (Thermodynamic Integration, BAR, Bennett)
- Replica exchange (Temperature, Hamiltonian)
- Normal Mode Analysis
- Essential Dynamics
- Non-equilibrium MD (pulling, flow)
- Coarse-grained models (Martini)

## Capabilities (CRITICAL)
- Extremely high performance (SIMD, GPU acceleration)
- Simulation of proteins, lipids, nucleic acids, polymers
- Advanced free energy methods
- Flexible force field support (AMBER, CHARMM, GROMOS, OPLS, Martini)
- Analysis tools (RMSD, RDF, H-bonds, clustering, etc.)
- Parallelization (MPI + OpenMP + GPU offload)
- Checkpointing and restarts
- Virtual sites for removing degrees of freedom

**Sources**: GROMACS documentation, SoftwareX 1-2, 19 (2015)

## Key Strengths

### Performance:
- World-leading single-node speed
- Excellent GPU utilization
- SIMD optimizations
- Dynamic load balancing

### Biomolecular Focus:
- Optimized for proteins/lipids/DNA
- Extensive analysis tools
- Free energy methods
- Coarse-grained (Martini)

### Ecosystem:
- PLUMED integration
- BioExcel workflows
- Extensive tutorials
- Professional support available

## Inputs & Outputs
- **Input formats**: pdb/gro (structure), top (topology), mdp (run parameters), tpr (compiled run input)
- **Output data types**: xtc/trr (trajectory), edr (energy), log (log file), gro (final structure)

## Interfaces & Ecosystem
- **PLUMED**: Native patching supported
- **VMD/PyMOL**: Visualization
- **Python**: gmxapi for Python control
- **BioExcel**: Integration in workflows
- **CP2K/QMCPACK**: QM/MM interfaces

## Workflow and Usage
1. Prepare topology: `gmx pdb2gmx -f protein.pdb -o protein.gro -p topol.top`
2. Define box and solvate: `gmx editconf`, `gmx solvate`
3. Add ions: `gmx genion`
4. Energy minimization: `gmx grompp`, `gmx mdrun`
5. Equilibration (NVT/NPT): `gmx grompp`, `gmx mdrun`
6. Production run: `gmx mdrun -v -deffnm prod`

## Performance Characteristics
- World-leading single-node performance
- Excellent GPU utilization (CUDA/OpenCL/SYCL)
- Efficient domain decomposition for MPI
- Dynamic load balancing

## Computational Cost
- Fastest single-node performance available
- GPU provides 10-50x speedup
- Efficient for 10K-10M atoms
- Overall: Industry-leading efficiency

## Best Practices
- Use latest version for GPU support
- Choose appropriate integrator (md, sd)
- Validate force field for your system
- Use checkpointing for long runs
- Check for LINCS/SETTLE warnings

## Limitations & Known Constraints
- Less flexible than LAMMPS for custom potentials
- Biomolecular focus (less materials)
- Complex topology format
- Steep learning curve for advanced features

## Application Areas
- Protein folding and dynamics
- Drug discovery (binding affinity)
- Membrane biophysics
- Polymer melts and solutions
- Micelles and self-assembly

## Comparison with Other Codes
- **vs LAMMPS**: GROMACS faster for biomolecules, LAMMPS more versatile potentials
- **vs AMBER**: GROMACS open-source, AMBER better GPU for some systems
- **vs NAMD**: GROMACS faster single-node, NAMD better multi-node scaling
- **Unique strength**: Fastest biomolecular MD, excellent free energy methods

## Community and Support
- Open-source (LGPL v2.1)
- Huge user community
- Active forum (gmx-users)
- Annual workshops and conferences
- Professional support available

## Verification & Sources
**Primary sources**:
1. Homepage: https://www.gromacs.org/
2. GitLab: https://gitlab.com/gromacs/gromacs
3. Publication: Abraham et al., SoftwareX 1-2, 19 (2015)

**Secondary sources**:
1. GROMACS tutorials
2. BioExcel documentation
3. Extensive published applications (>50,000 citations)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitLab)
- Development: ACTIVE (KTH, Max Planck, etc.)
- Applications: MD, biochemistry, high performance, GPU acceleration
