# fcDMFT

## Official Resources
- Homepage: https://github.com/ZhuGroup-Yale/fcdmft
- Documentation: https://github.com/ZhuGroup-Yale/fcdmft/tree/master/examples
- Source Repository: https://github.com/ZhuGroup-Yale/fcdmft
- License: Open Source (check repo for specific license, often GPL or MIT)

## Overview
fcDMFT is a Python-based software package designed for *ab initio* Full Cell Dynamical Mean-Field Theory (DMFT) calculations. Built upon the PySCF quantum chemistry framework, it extends standard embedding theories to treat solid-state systems with high accuracy. It focuses on GW+DMFT and HF+DMFT methodologies, enabling the study of electronic correlations in periodic crystals using quantum chemical solvers without downfolding to a small subspace.

**Scientific domain**: Condensed matter physics, Quantum chemistry of solids, GW+DMFT
**Target user community**: Researchers in electronic structure theory, quantum chemistry, and correlated materials

## Theoretical Methods
- Full Cell DMFT (no downfolding to minimal basis)
- GW+DMFT (G0W0+DMFT)
- Hartree-Fock + DMFT (HF+DMFT)
- Quantum Chemical Impurity Solvers (CCSD-GF, DMRG, FCI)
- Periodic RPA and GW
- CAS-CI impurity treatment

## Capabilities (CRITICAL)
- **Full Cell Embedding**: Treats the entire unit cell in the embedding scheme, avoiding some supercell approximations and basis set truncation errors typical of downfolding.
- **GW Integration**: Combines GW quasiparticle energies with DMFT correlations for superior spectral properties.
- **Quantum Chemistry Solvers**: Utilizes solvers like Coupled Cluster Green's Function (CCGF), DMRG (via `Block2`), and FCI (via `CheMPS2`).
- **PySCF Integration**: leveraged PySCF for integrals, mean-field methods, and periodic boundary conditions.
- **Parallelization**: Mixed MPI and OpenMP parallelism for performance.

## Key Features

### Quantum Chemistry Solvers:
- Interfaces with `Block2` for DMRG calculations.
- Interfaces with `CheMPS2` for Full Configuration Interaction (FCI).
- Implements Coupled-Cluster Green's Function (CCSD-GF) solvers naturally.

### Full Cell Approach:
- Performs DMFT calculations in the full unit cell basis, capturing non-local correlations more effectively than single-site approximations.

### GW+DMFT:
- Advanced integration of GW screening with DMFT local correlations for accurate spectral properties.

## Inputs & Outputs
- **Input formats**:
  - Python scripts using PySCF objects (e.g., `run_dmft.py`).
  - Geometry and basis set definitions via PySCF `gto.Cell`.
  - Scripts like `si_gw.py` for precursor GW steps.
- **Output data types**:
  - Green's functions (HDF5 or numpy arrays)
  - Self-energies (frequency dependent)
  - Quasiparticle energies
  - Spectral functions

## Interfaces & Ecosystem
- **PySCF**: Core dependency for integral generation, mean-field (DFT/HF), and PBC tools.
- **Solvers**: Interfaces with `Block2` (DMRG) and `CheMPS2` (FCI).
- **Python**: Fully scriptable Python API for custom workflows.

## Workflow and Usage
Users write Python scripts (examples in `/fcdmft/examples`):
1.  Define cell and basis set using PySCF (`gto.Cell`).
2.  Perform a mean-field (DFT or HF) or GW calculation (`si_gw.py`).
3.  Derive the impurity Hamiltonian and GW double-counting terms (`si_set_ham.py`).
4.  Construct the impurity Hamiltonian.
5.  Invoke the fcDMFT solver loop (`run_dmft.py`).

## Performance Characteristics
- **Scalability**: Uses MPI/OpenMP to scale on clusters.
- **Cost**: Computationally intensive due to full cell treatment and quantum chemical solvers, but offering high accuracy.


## Comparison with Other Codes
- **vs EDMFTF**: EDMFTF (Wien2k) assumes locality in real space and focuses on forces; fcDMFT (PySCF) focuses on quantum chemistry solvers and full unit cell embedding.
- **vs Standard DMFT**: fcDMFT avoids downfolding to a minimal basis, treating the full cell, which is more expensive but potentially more accurate.
- **Unique strength**: Seamless integration with quantum chemistry solvers (DMRG, CCSD) via PySCF.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/ZhuGroup-Yale/fcdmft
2. Example scripts: `/fcdmft/examples/Si` in the repository
3. Associated Publications: Research from the Zhu Group (Yale) on Full Cell DMFT methods.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Integration: Deeply integrated with PySCF
- Active development: Research code from Yale group
