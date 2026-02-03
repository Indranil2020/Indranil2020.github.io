# NEB (Nudged Elastic Band)

## Official Resources
- Homepage: Method - Implemented in many codes (VASP, ASE, LAMMPS, VTST)
- Documentation: https://theory.cm.utexas.edu/vtsttools/neb.html (VTST Tools)
- Source Repository: https://github.com/henkelmanlab/vtstscripts
- License: Varies by implementation (MIT for VTST scripts)

## Overview
Nudged Elastic Band (NEB) is a method for finding the minimum energy path (MEP) between two stable states of a system, typically used to determine transition states and activation barriers for chemical reactions, diffusion processes, and phase transitions. It is not a single software code but a method implemented in virtually all major electronic structure and molecular dynamics packages (VASP, ASE, LAMMPS, DFTB+, etc.).

**Scientific domain**: Reaction pathways, transition states, activation energy, saddle point search  
**Target user community**: Computational chemists, materials scientists, catalysis researchers

## Theoretical Methods
- Nudged Elastic Band (NEB)
- Climbing Image NEB (CI-NEB) for precise saddle points
- Doubly Nudged Elastic Band (DNEB)
- Free Energy NEB
- Tangent estimation
- Optimization algorithms (L-BFGS, Quick-Min, FIRE)

## Capabilities (CRITICAL)
- Finding Minimum Energy Paths (MEP)
- Locating transition states (saddle points)
- Calculating activation energy barriers
- Investigating reaction mechanisms
- Atomic diffusion pathways
- Solid-state phase transitions
- Method implemented in: VASP (via VTST), ASE, LAMMPS, CP2K, Quantum ESPRESSO, GPAW, etc.

**Sources**: G. Henkelman et al., J. Chem. Phys. 113, 9901 (2000)

## Inputs & Outputs
- **Input formats**: Initial and Final structures (POSCAR/xyz), NEB parameters (spring constant, optimizer)
- **Output data types**: Image structures along path, energy barrier profile, tangent forces

## Interfaces & Ecosystem
- **VTST Tools**: Standard implementation for VASP
- **ASE**: Python-based NEB implementation compatible with many calculators
- **LAMMPS**: NEB for classical potentials
- **Transition Path Theory**: Related methods

## Workflow and Usage
1. Optimize Initial (IS) and Final (FS) states
2. Generate initial path: Interpolate images (linear or IDPP)
3. Run NEB: Relax images perpendicular to the path
4. Climbing Image: Turn on CI-NEB to converge to saddle point
5. Analysis: Plot energy vs reaction coordinate

## Performance Characteristics
- Computationally expensive: Requires N images × optimization steps
- Parallelization: Images can be run in parallel (one image per group)
- Convergence: Can be slow on flat potential energy surfaces

## Application Areas
- Catalysis (reaction barriers)
- Diffusion in batteries (ion hopping)
- Surface diffusion
- Chemical reaction mechanisms
- Conformational changes

## Community and Support
- Henkelman Group (UT Austin) for VTST
- ASE community
- VASP forum
- Widely used and cited method

## Verification & Sources
**Primary sources**:
1. VTST Tools: https://theory.cm.utexas.edu/vtsttools/
2. ASE NEB: https://wiki.fysik.dtu.dk/ase/ase/neb.html
3. Publication: G. Henkelman, B.P. Uberuaga, and H. Jonsson, J. Chem. Phys. 113, 9901 (2000)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (VTST/ASE)
- Documentation: COMPREHENSIVE
- Method: STANDARD in field
- Applications: Transition states, reaction barriers, diffusion
