# PSI4

## Official Resources
- Homepage: https://psicode.org/
- Documentation: https://psicode.org/psi4manual/master/
- Source Repository: https://github.com/psi4/psi4
- License: GNU Lesser General Public License v3.0

## Overview
PSI4 is an open-source suite of ab initio quantum chemistry programs designed for efficient, high-accuracy simulations of molecular properties. It emphasizes modern software engineering practices, native Python integration, and provides state-of-the-art coupled cluster, density functional, and symmetry-adapted perturbation theory methods. The latest version is PSI4 1.9.1 (February 2024).

**Scientific domain**: Quantum chemistry, molecular properties, method development, non-covalent interactions  
**Target user community**: Quantum chemists, researchers studying intermolecular interactions

## Theoretical Methods
- Hartree-Fock (RHF, UHF, ROHF)
- Density Functional Theory (DFT)
  - LDA, GGA, meta-GGA, hybrid functionals
  - Double-hybrid functionals
  - Extensive Libxc integration
- Møller-Plesset (MP2, MP3, MP4)
- Density-fitted MP2 (DF-MP2) and DLPNO-MP2
- Coupled Cluster (CCSD, CCSD(T), FNOCC)
- Orbital-optimized methods (OO-MP2, OO-CCSD)
- Symmetry-Adapted Perturbation Theory (SAPT0 to SAPT2+3)
  - F/I-SAPT (functional-group/intramolecular)
  - High-spin open-shell SAPT0 (v1.9+)
- Algebraic Diagrammatic Construction (ADC)
- Equation-of-motion coupled cluster (EOM-CCSD)
- Time-Dependent DFT/HF (TDSCF)
- Multi-configurational SCF (MCSCF)
- Full Configuration Interaction (FCI)
- Density Cumulant Theory (DCT)
- Dispersion corrections (DFT-D3, DFT-D4)
- Solvation models (PCM via PCMSolver)
- Scalar relativistic Hamiltonians

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Geometry optimization and transition states
- Vibrational frequencies and thermochemistry
- Excited states (TDDFT, EOM-CC, ADC)
- Intermolecular interaction analysis (SAPT)
- Energy decomposition analysis (EDA)
- Molecular properties (dipole, quadrupole, polarizability)
- NMR chemical shifts
- Response properties
- Orbital analysis
- Natural Bond Orbital (NBO) analysis via interface
- Basis set extrapolation (automatic)
- Composite methods
- Density fitting for efficiency
- GPU acceleration (BrianQC interface)
- Native Python API
- QCSchema standardized I/O
- Integration with MolSSI QCArchive infrastructure

**Sources**: Official PSI4 documentation, cited in 7/7 source lists

## Key Strengths

### SAPT Excellence:
- State-of-the-art SAPT implementations
- F/I-SAPT for functional group analysis
- Detailed energy decomposition
- Non-covalent interaction specialization
- Benchmark-quality interaction energies

### Python Integration:
- Native Python API (not wrapper)
- Psi4NumPy for educational purposes
- Jupyter notebook support
- Seamless NumPy/SciPy integration
- Scriptable workflows

### Open-Source Ecosystem:
- LGPL v3 license
- Active GitHub development
- MolSSI QCARCHIVE integration
- Plugin architecture
- Community contributions

### Modern Software:
- Clean C++ and Python codebase
- Plugin system for extensions
- Continuous integration
- Extensive testing suite
- Active maintenance

## Inputs & Outputs
- **Input formats**:
  - Python scripts (native interface)
  - Simple input files (.dat)
  - XYZ coordinate files
  - Z-matrix input
  - QCSchema JSON format
  
- **Output data types**:
  - Detailed output files
  - Energies, gradients, Hessians
  - Molecular orbitals (Molden format)
  - Checkpoint files
  - Wavefunction objects in Python
  - JSON output via QCSchema

## Interfaces & Ecosystem
- **Python integration**:
  - Native Python API
  - Psi4NumPy educational modules
  - Integration with NumPy, SciPy
  - Jupyter notebook support
  
- **External programs**:
  - CFOUR interface for high-level CC
  - MRCC interface for arbitrary-order CC
  - CheMPS2 for DMRG
  - Libxc for DFT functionals
  - BrianQC for GPU acceleration
  - PCMSolver for solvation
  
- **Workflow tools**:
  - OptKing for geometry optimization
  - QCEngine for standardized I/O
  - QCArchive for distributed computing
  - ASE interface


## Workflow and Usage

### Python API (Recommended):
PSI4 is most powerful when used as a Python library.

```python
import psi4

# Define molecule
psi4.set_memory('4 GB')
mol = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

# Run calculation
en, wfn0 = psi4.energy('pbe0/def2-svp', return_wfn=True)
psi4.optimize('pbe0/def2-svp')
```

### Psithon Input Format:
A Python-like input file format for standalone execution.

```python
molecule {
  O
  H 1 0.96
  H 1 0.96 2 104.5
}

set basis def2-svp
energy('pbe0')
```

### Running PSI4:
```bash
psi4 input.dat output.dat
```

## Advanced Features

### Symmetry-Adapted Perturbation Theory (SAPT):
- Decomposes interaction energy into electrostatic, exchange, induction, and dispersion
- SAPT0, SAPT2, SAPT2+, SAPT2+3 methods
- F/I-SAPT for functional group analysis
- Intramolecular non-covalent interactions

### CBS Extrapolation:
- Automated Complete Basis Set extrapolation
- `energy('ccsd(t)/cbs')` syntax
- Configurable schemes (Helgaker, Feller, etc.)
- Essential for high-accuracy thermochemistry

### Density Cumulant Theory (DCT):
- ODC-12 method
- Description of static correlation
- Alternative to multireference methods
- Efficient implementation

### Orbital Optimization:
- OO-MP2 and OO-CCSD
- Addresses spin contamination in open-shell systems
- Improved results for radicals and bond breaking
- Variational optimization of orbitals

### QCArchive Integration:
- Native support for MolSSI QCArchive
- Systematic data generation
- Distributed computing support
- Machine learning dataset creation

## Performance Characteristics
- **Speed**: Efficient density fitting (DF) makes MP2/CCSD faster than many competitors
- **Scalability**: Shared memory parallelism (OpenMP) is good; distributed memory (MPI) is limited compared to NWChem
- **Memory**: Can handle large basis sets with DF approximations
- **Python Overhead**: Minimal, core is C++
- **Bottlenecks**: Disk I/O for large coupled cluster calculations

## Computational Cost
- **DFT**: Efficient, especially with density fitting
- **SAPT0**: O(N^5), tractable for ~50-100 atoms
- **SAPT2+3**: O(N^7), very expensive
- **DLPNO-MP2**: Linear scaling (very efficient)
- **Canonical CCSD(T)**: O(N^7), expensive
- **Education**: Very low cost for small pedagogical examples

## Comparison with Other Codes
- **vs Gaussian**: PSI4 is open-source, better for non-covalent interactions (SAPT); Gaussian has more functionals and solvation models.
- **vs ORCA**: Both strong interactions focus; PSI4 has native Python API, ORCA has DLPNO-CCSD(T) (PSI4 has DLPNO-MP2).
- **vs PySCF**: Both Python-based; PSI4 higher level "black box" quantum chemistry; PySCF more "toolkit" style for tensors/solids.
- **vs Molpro**: Both strong in coupled cluster; PSI4 is free/open-source.
- **Unique strength**: Best-in-class SAPT, native Python API for complex workflows, QCArchive support.

## Best Practices

### Memory Management:
- Always set memory explicitly (`psi4.set_memory`)
- Default is often too low for CC methods
- Use density fitting (`set scf_type df`) for speed

### SAPT Calculations:
- Use `sapt0` for large systems
- Use `sapt2+3` only for small benchmarks
- Use mixed basis sets (e.g., jun-cc-pVDZ) for balanace

### Python Scripting:
- Loop over geometries for PES scans
- Use Python for post-processing results
- Save wavefunctions for restart

## Community and Support
- **Forum**: Active discourse forum
- **GitHub**: Transparent development, easy contribution
- **Tutorials**: Extensive Psi4NumPy collection
- **Education**: Widely used in computational chemistry courses
- **MolSSI**: Supported by Molecular Sciences Software Institute

## Application Areas

### Non-Covalent Interactions:
- Hydrogen bonding analysis
- π-stacking interactions
- Van der Waals complexes
- Host-guest chemistry
- Drug-receptor binding

### Method Development:
- New method implementation
- Plugin architecture
- Benchmarking studies
- Algorithm testing

### Education:
- Psi4NumPy tutorials
- Transparent code
- Interactive learning
- Web-based UI (Psi4-WebUI)

## Limitations & Known Constraints
- **Molecular focus**: Not designed for periodic systems
- **System size**: Limited by CC scaling; ~50-100 atoms for DFT
- **Basis sets**: Gaussian-type only
- **Parallelization**: Threading and limited MPI; varies by method
- **Memory**: High-level methods memory-intensive
- **Learning curve**: Moderate; Python knowledge helpful
- **Documentation**: Excellent but assumes quantum chemistry background
- **Platform**: Linux, macOS; Windows via WSL

## Verification & Sources
**Primary sources**:
1. Official website: https://psicode.org/
2. Documentation: https://psicode.org/psi4manual/master/
3. GitHub repository: https://github.com/psi4/psi4
4. D. G. A. Smith et al., J. Chem. Phys. 152, 184108 (2020) - PSI4 1.4
5. R. M. Parrish et al., J. Chem. Theory Comput. 13, 3185 (2017) - PSI4 1.1
6. Latest release: PSI4 1.9.1 (February 2024)

**Secondary sources**:
1. PSI4 tutorials and workshops
2. Psi4NumPy educational modules
3. QCArchive integration documentation
4. Confirmed in 7/7 source lists (claude, g, gr, k, m, q, z)

**Confidence**: CONFIRMED - Appears in all 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, LGPL v3)
- Community support: Very active (forum, GitHub)
- Academic citations: >2,500 (various versions)
- Active development: Regular releases, modern codebase
- Specialized strength: SAPT, Python integration, non-covalent interactions, open-source
