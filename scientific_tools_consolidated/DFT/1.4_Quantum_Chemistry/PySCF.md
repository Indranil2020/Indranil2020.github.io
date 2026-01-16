# PySCF

## Official Resources
- Homepage: https://pyscf.org/
- Documentation: https://pyscf.org/user.html
- Source Repository: https://github.com/pyscf/pyscf
- License: Apache License 2.0

## Overview
PySCF (Python-based Simulations of Chemistry Framework) is an open-source quantum chemistry package with emphasis on ab initio methods for molecules and crystals. Built entirely in Python with performance-critical sections in C, it provides a simple, lightweight, and efficient platform for developing and testing quantum chemistry methods with excellent scriptability and integration with the Python scientific ecosystem.

**Scientific domain**: Quantum chemistry, method development, molecular and periodic systems  
**Target user community**: Quantum chemists, method developers, researchers needing Python-scriptable quantum chemistry

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF, GHF)
- Density Functional Theory (DFT)
  - LDA, GGA, meta-GGA, hybrid functionals
  - Range-separated hybrids
  - Double-hybrid functionals
  - Extensive Libxc integration
- Møller-Plesset (MP2, MP3)
- Coupled Cluster (CCSD, CCSD(T), EOM-CCSD, IP/EA-EOM-CC)
- Multireference methods (CASCI, CASSCF, NEVPT2, ic-MRCI)
- Full Configuration Interaction (FCI)
- Density Matrix Renormalization Group (DMRG via Block/CheMPS2)
- Selected CI (SHCI, heat-bath CI)
- Algebraic Diagrammatic Construction (ADC)
- GW approximation (G0W0, evGW)
- Random Phase Approximation (RPA)
- Time-Dependent DFT (TDDFT)
- Relativistic methods (X2C, DKH, Breit-Pauli SO)
- Periodic boundary conditions (PBC) with k-point sampling
- Density Matrix Embedding Theory (DMET)

## Capabilities (CRITICAL)
- Ground-state electronic structure (molecules and solids)
- Geometry optimization
- Vibrational frequencies
- Excited states (TDDFT, EOM-CC, ADC)
- GW quasiparticle energies
- Periodic systems (Gamma-point and full k-point sampling)
- Gaussian and mixed Gaussian/plane-wave basis
- Molecular properties (dipole, quadrupole, polarizability)
- NMR shielding tensors
- Spin-orbit coupling
- Embedding methods (DMET, QM/MM)
- Multi-reference calculations
- FCI for benchmarking
- Custom method development via Python API
- Interface to external codes (DMRG, SHCI)
- GPU acceleration (GPU4PySCF for DFT)
- Integrals via libcint (efficient)

**Sources**: Official PySCF documentation, cited in 7/7 source lists

## Key Strengths

### Python Integration:
- Native Python interface
- NumPy/SciPy integration
- Jupyter notebook friendly
- Easy scripting and automation
- Rapid prototyping

### Method Development:
- Clean, readable codebase
- Modular architecture
- Easy to extend
- Large contributor community
- Active research platform

### Periodic Systems:
- Full k-point sampling
- Gamma-point approximation
- Periodic DFT, HF, MP2, CCSD
- Crystalline solids
- Band structure calculations

### Open-Source Ecosystem:
- Apache 2.0 license
- Active GitHub development
- pyscf-forge extension modules
- PySCF-NAO, PySCF-properties
- GPU4PySCF for acceleration

## Inputs & Outputs
- **Input formats**:
  - Python scripts (native interface)
  - XYZ coordinate files
  - Molecular input via ASE, pymatgen
  - Checkpoint files for restart
  
- **Output data types**:
  - Python objects (energies, wavefunctions)
  - Molecular orbitals
  - Density matrices
  - Property calculations
  - Checkpoint files (.chk)
  - FCIDUMP format for external codes

## Interfaces & Ecosystem
- **Framework integrations**:
  - ASE calculator interface
  - pymatgen structure handling
  - pyscf-forge extension modules
  - PySCF-NAO - numerical atomic orbitals
  
- **External codes**:
  - Block/Block2 for DMRG
  - CheMPS2 for DMRG
  - Dice for SHCI
  - FCIDUMP for external CI codes
  - Molden format for visualization
  
- **Python ecosystem**:
  - NumPy, SciPy for linear algebra
  - h5py for HDF5 I/O
  - Jupyter notebooks
  - mpi4py for parallel execution


## Workflow and Usage

### Python Scripting (Standard):
PySCF is designed to be used as a Python library.

```python
from pyscf import gto, scf, cc

# Define molecule
mol = gto.M(
    atom = 'O 0 0 0; H 0 1 0; H 0 0 1',
    basis = 'ccpvdz',
    symmetry = True
)

# Mean-field calculation
mf = scf.RHF(mol).run()

# Post-Hartree-Fock
mycc = cc.CCSD(mf).run()
et = mycc.ccsd_t()
```

### Common Objects:
- `gto.M`: Molecule/Cell object
- `scf.RHF`/`scf.RKS`: Mean-field objects
- `cc.CCSD`: Coupled cluster object
- `pbc.gto.Cell`: Periodic cell object

## Advanced Features

### Periodic Boundary Conditions (PBC):
- `pyscf.pbc` module
- Full k-point sampling
- Gamma-point optimization
- Gaussian function density fitting (GDF)
- Pseudopotentials (GTH)

### Custom Hamiltonians:
- Ideal for model systems (Hubbard, Heisenberg)
- Easy modification of core Hamiltonian
- Embedding environments
- Quantum computing interfaces

### Embedding Methods:
- Density Matrix Embedding Theory (DMET)
- Projection-based embedding
- QM/MM interfaces
- Local correlation approaches

### Relativistic Effects:
- `mol.x2c()` for X2C Hamiltonian
- `scf.DHF` for 4-component Dirac-Hartree-Fock
- Spin-orbit coupling handling
- Model potentials

### Molecular Dynamics:
- Born-Oppenheimer MD (via `pyscf.md`)
- Geometric integration
- Gradients for most methods
- Ab initio MD capability

## Performance Characteristics
- **Speed**: Optimized C backend (libcint, libxc) ensures competitive speed for integrals and SCF
- **Scalability**: MPI support via `mpi4py` (good, but less scalable than NWChem for massive jobs)
- **Memory**: Can be memory hungry for large CC calculations
- **Tensor Operations**: Leverages NumPy/block-sparse operations
- **GPU**: GPU4PySCF extension offers massive speedups for DFT

## Computational Cost
- **DFT**: Efficient, O(N^3) or better with density fitting
- **CCSD(T)**: O(N^7), standard scaling
- **Periodic DFT**: Expensive with large k-point grids
- **FCI**: Exponential cost, limited to ~16 electrons/orbitals
- **DMRG**: Polynomial cost, feasible for larger active spaces

## Comparison with Other Codes
- **vs PSI4**: PySCF is more "toolkit" oriented, better for custom method development and periodic systems; PSI4 has better SAPT.
- **vs VASP**: PySCF uses Gaussian basis for solids (all-electron or PP), VASP uses Plane-Waves; VASP is industry standard for solids, PySCF offers quantum chemistry methods (CC) for solids.
- **vs Gaussian**: PySCF is free/open-source, scriptable; Gaussian is a black-box commercial tool.
- **Unique strength**: Unmatched scriptability (Python), excellent handling of both molecular and periodic systems in one framework.

## Best Practices

### Memory Usage:
- Set `mol.max_memory` (in MB)
- Use density fitting (`mf.density_fit()`) to save memory/time
- Use `pyscf.lib.num_threads()` to control parallelism

### Periodic Calculations:
- Check k-point convergence
- Use `GDF` (Gaussian Density Fitting) for stability
- Verify pseudopotentials (GTH)

### Method Development:
- Use `pyscf.tools` for debugging
- Leverage `numpy` integration
- Check `examples/` directory (very comprehensive)

## Community and Support
- **GitHub**: Very active development and issue tracking
- **Documentation**: Extensive API reference and user guide
- **Examples**: Hundreds of example scripts in the repository
- **Discord/Slack**: Active developer community
- **License**: Apache 2.0 (permissive open source)

## Application Areas

### Method Development:
- New algorithm implementation
- Rapid prototyping
- Testing and benchmarking
- Educational purposes
- Research extensions

### Periodic Systems:
- Solid-state DFT
- Band structure calculations
- Crystalline materials
- Periodic coupled cluster

### Embedding:
- Density Matrix Embedding Theory
- QM/MM calculations
- Multiscale modeling
- Local correlation

## Limitations & Known Constraints
- **Python overhead**: Slower than pure Fortran/C++ for routine HPC calculations
- **Memory**: Python memory management overhead for very large calculations
- **Parallelization**: Threading and MPI available but not as optimized as commercial codes
- **Basis sets**: Primarily Gaussian-type; via basis_set_exchange
- **System size**: Molecular focus; periodic features maturing
- **Documentation**: Good but assumes Python and quantum chemistry knowledge
- **Learning curve**: Requires Python programming skills
- **Platform**: Primarily Linux/macOS; Windows requires care

## Verification & Sources
**Primary sources**:
1. Official website: https://pyscf.org/
2. Documentation: https://pyscf.org/user.html
3. GitHub repository: https://github.com/pyscf/pyscf
4. Q. Sun et al., WIREs Comput. Mol. Sci. 8, e1340 (2018) - PySCF paper
5. Q. Sun et al., J. Chem. Phys. 153, 024109 (2020) - Recent developments

**Secondary sources**:
1. PySCF tutorials and examples
2. Published method development using PySCF
3. Jupyter notebook examples
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, Apache 2.0)
- Community support: Very active (GitHub, Google group)
- Academic citations: >2,000 (main paper)
- Active development: Continuous updates, large contributor base
- Specialized strength: Python integration, method development, periodic systems, embedding, open-source
