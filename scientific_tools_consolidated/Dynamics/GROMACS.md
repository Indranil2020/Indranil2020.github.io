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

## Application Areas
- Protein folding and dynamics
- Drug discovery (binding affinity)
- Membrane biophysics
- Polymer melts and solutions
- Micelles and self-assembly

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

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitLab)
- Development: ACTIVE (KTH, Max Planck, etc.)
- Applications: MD, biochemistry, high performance, GPU acceleration
