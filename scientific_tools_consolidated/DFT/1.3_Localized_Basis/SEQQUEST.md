# SEQQUEST

## Official Resources
- Homepage: https://dft.sandia.gov/quest/
- Documentation: https://dft.sandia.gov/quest/SeqsDocumentationPage.html
- Developer: Sandia National Laboratories
- License: Available to friendly users (contact required)

## Overview
SEQQUEST (Sequential Quantum Electronic Structure Tool) is a general-purpose electronic structure code developed at Sandia National Laboratories for computing energies and forces of molecules, periodic surfaces (slabs), and bulk solids. It uses the LCAO approach with contracted-Gaussian basis sets and norm-conserving pseudopotentials, with linear-scaling algorithms for Hamiltonian generation.

**Scientific domain**: Surfaces, solids, molecules, defects, interfaces  
**Target user community**: Materials scientists studying surface chemistry, catalysis, and solid-state materials

## Theoretical Methods
- Density Functional Theory (DFT)
- Linear Combination of Atomic Orbitals (LCAO)
- Contracted-Gaussian basis sets
- Norm-conserving pseudopotentials
- LDA and GGA exchange-correlation functionals
- Spin-polarized calculations
- Linear-scaling Hamiltonian construction
- Self-consistent field (SCF)

## Capabilities (CRITICAL)
- Ground-state electronic structure
- Total energy calculations
- Force calculations
- Geometry optimization
- Periodic systems (1D, 2D, 3D)
- Surface slabs
- Point defects
- Molecular calculations
- Spin-polarized magnetism
- Band structure and DOS
- Charge density analysis
- Linear-scaling Hamiltonian generation

**Sources**: Sandia National Laboratories, DOE Office of Science

## Key Strengths

### Linear-Scaling Hamiltonian:
- O(N) Hamiltonian matrix construction
- Efficient for large systems
- Exploits sparsity
- Parallel implementation

### LCAO Efficiency:
- Compact Gaussian basis sets
- Localized description
- Efficient matrix elements
- Reduced computational cost

### Surface Science Focus:
- Optimized for slabs
- Adsorption studies
- Surface reconstructions
- Adsorbate-surface interactions

### Sandia Quality:
- Extensively tested
- Production-ready
- DOE-funded development
- Reliable results

## Inputs & Outputs
- **Input formats**:
  - Native SEQQUEST input format
  - Atomic coordinates
  - Basis set files
  - Pseudopotential files
  
- **Output data types**:
  - Total energies
  - Forces
  - Optimized structures
  - Band structure
  - DOS
  - Charge densities

## Interfaces & Ecosystem
- **Web portals**:
  - NSF nanoHUB access
  - memsHUB integration
  - Web-based job submission
  
- **Preprocessing**:
  - Structure generation tools
  - Slab builders
  - Defect constructors

## Advanced Features

### Surface Calculations:
- Slab geometry optimization
- Adsorption energies
- Surface energy calculations
- Work function determination
- Surface reconstructions

### Defect Studies:
- Point defects in solids
- Vacancy formation energies
- Dopant calculations
- Defect levels

### Magnetic Systems:
- Spin-polarized DFT
- Magnetic moment calculations
- Ferromagnetic ordering
- Antiferromagnetic states

## Performance Characteristics
- **Speed**: Efficient LCAO implementation
- **Accuracy**: Reliable DFT results
- **System size**: Hundreds to thousand+ atoms
- **Memory**: Moderate requirements
- **Parallelization**: MPI parallel

## Computational Cost
- **Linear scaling**: Hamiltonian generation O(N)
- **SCF**: Cubic diagonalization
- **Typical**: Workstation to cluster
- **Large systems**: Efficient on modest resources

## Limitations & Known Constraints
- **Availability**: Not freely distributed
- **Hybrid functionals**: Limited support
- **All-electron**: Pseudopotential only
- **Documentation**: Internal focus
- **Community**: Sandia-centered
- **Learning curve**: Requires training

## Comparison with Other Codes
- **vs SIESTA**: Similar LCAO approach
- **vs VASP**: SEQQUEST LCAO vs VASP plane-wave
- **vs Gaussian**: SEQQUEST materials focus
- **Unique strength**: Sandia-quality, surface science, linear-scaling Hamiltonian

## Application Areas

### Catalysis:
- Surface reactions
- Adsorbate binding
- Reaction barriers
- Catalyst design

### Semiconductor Defects:
- Dopant energetics
- Native defects
- Defect complexes
- Carrier concentrations

### Surface Chemistry:
- Oxidation
- Corrosion initiation
- Surface passivation
- Interface formation

### Energy Materials:
- Battery materials
- Solar cell interfaces
- Fuel cell catalysts
- Hydrogen storage

## Best Practices

### Basis Set Optimization:
- Use appropriate contraction
- Test convergence
- Balance accuracy and cost

### Surface Models:
- Sufficient slab thickness
- Vacuum gap size
- Symmetric vs asymmetric slabs
- k-point convergence

## Community and Support
- Sandia National Laboratories
- DOE funding
- nanoHUB/memsHUB access
- Internal documentation

## Verification & Sources
**Primary sources**:
1. Sandia DFT: https://dft.sandia.gov/quest/
2. DOE Office of Science: Sandia computational materials

**Confidence**: VERIFIED - Sandia National Laboratories official code

**Verification status**: ✅ VERIFIED
- Source code: Restricted (Sandia)
- Academic use: Collaborator access
- Documentation: Internal
- Production use: DOE projects
- Specialty: Surface science, linear-scaling Hamiltonian
