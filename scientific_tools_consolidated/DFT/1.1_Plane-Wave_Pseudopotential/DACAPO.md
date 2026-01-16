# DACAPO

## Official Resources
- Homepage: https://wiki.fysik.dtu.dk/dacapo/
- Documentation: https://wiki.fysik.dtu.dk/dacapo/
- Source Repository: Historic repositories at DTU; integrated into ASE history
- License: Open Source (GPL)

## Overview
DACAPO is a pioneering plane-wave pseudopotential Density Functional Theory (DFT) code. Developed within the CAMPOS project at the Technical University of Denmark (DTU), it was the original electronic structure backend for the Atomic Simulation Environment (ASE). It utilizes ultrasoft pseudopotentials and was heavily used for surface science and catalysis research, laying the groundwork for many modern computational workflows.

**Scientific domain**: Surface science, catalysis, materials physics
**Target user community**: Historic user base in catalysis; largely superseded by GPAW and Quantum ESPRESSO but remains a reference code.

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Ultrasoft Pseudopotentials (USPP)
- Vanderbilt construction
- LDA and GGA functionals (PW91, RPBE)
- Iterative diagonalization (Davidson, Pulay)

## Capabilities
- Ground-state energy and forces
- Geometry optimization (via ASE)
- Molecular Dynamics (MD)
- Transition state search (Nudged Elastic Band via ASE)
- Surface adsorption energy calculations
- Parallelization (MPI)

## Key Strengths

### Historic Significance:
- Pioneered the python-driven workflow approach via ASE.
- Developed specifically for surface science problems (slabs, adsorption).

### Methods:
- Robust implementation of Ultrasoft Pseudopotentials allowing lower cutoffs.
- Excellent convergence for transition metals.

## Inputs & Outputs
- **Input formats**:
  - Python scripts (via ASE)
  - NetCDF files for densities and potentials
  
- **Output data types**:
  - NetCDF output files
  - Text logs
  - Trajectory files (ASE .traj)

## Interfaces & Ecosystem
- **ASE**: The Atomic Simulation Environment was originally built *around* DACAPO.
- **NetCDF**: Uses NetCDF for robust binary I/O of large fields.

## Performance Characteristics
- **Speed**: competitive in its era; highly optimized Fortran.
- **Parallelization**: MPI parallelization over k-points and domain decomposition.

## Computational Cost
- **Efficiency**: Highly efficient for its time due to Ultrasoft Pseudopotentials (lower energy cutoffs required).
- **Scaling**: Good MPI scaling on legacy clusters; not optimized for modern GPU acceleration.
- **Memory**: Standard plane-wave memory requirements ($O(N^2)$ storage for wavefunctions).

## Best Practices

### Legacy Workflow:
- **Use via ASE**: Always drive DACAPO calculations using the ASE Python interface, which preserves the input/output logic.
- **Pseudopotentials**: Ensure you have the compatible USPP files (often hard to find now); sticking to the `cmp` database included with old distributions is recommended.

## Community and Support
- **Hosting**: Archived on [DTU Physics](https://wiki.fysik.dtu.dk/dacapo/).
- **Status**: **Legacy/Archived**. Active development ceased ~2010.
- **Support**: Limited. Queries are best directed to the [ASE Users mailing list](https://listserv.fysik.dtu.dk/mailman/listinfo/ase-users) for historic knowledge.

## Verification & Sources
**Primary sources**:
1. DTU Wiki: https://wiki.fysik.dtu.dk/dacapo/
2. "The CAMPOS Project" history
3. User experience with ASE < 3.0

**Confidence**: VERIFIED - Code is a known historic pillar of the field.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: Historic Surface Science
- Key Feature: ASE Backend

## Limitations & Known Constraints
- **Maintenance**: Effectively legacy software. Most development has moved to GPAW or standardizing interfaces to QE/VASP.
- **Installation**: Can be difficult to compile on modern systems due to older library dependencies.

## Comparison with Other Codes
- **vs GPAW**: GPAW is the modern successor from the same DTU group, offering PAW/Grid/LCAO modes and modern Python integration.
- **vs VASP**: Similar capabilities for surfaces; DACAPO was the open-source alternative in the early 2000s.

