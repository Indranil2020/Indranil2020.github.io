# OpenMolcas

## Official Resources
- Homepage: https://www.molcas.org/
- Documentation: https://molcas.gitlab.io/OpenMolcas/sphinx/
- Source Repository: https://gitlab.com/Molcas/OpenMolcas
- License: GNU Lesser General Public License v2.1

## Overview
OpenMolcas is an open-source quantum chemistry software package with special emphasis on multiconfigurational methods. It is the successor to Molcas and excels at treating systems with strong electron correlation, excited states, and complex electronic structures requiring multi-reference approaches. It is particularly strong in photochemistry and spectroscopy.

**Scientific domain**: Quantum chemistry, multiconfigurational methods, photochemistry, spectroscopy  
**Target user community**: Researchers studying multi-reference systems, excited states, and photochemical processes

## Theoretical Methods
- Complete Active Space SCF (CASSCF)
- Restricted Active Space SCF (RASSCF)
- Multi-Configurational Perturbation Theory (CASPT2, RASPT2)
- Multi-State CASPT2 (MS-CASPT2)
- Extended Multi-State CASPT2 (XMS-CASPT2)
- Multi-Reference CI (MRCI)
- Density Matrix Renormalization Group (DMRG)
- Hartree-Fock and DFT
- Coupled Cluster (CCSD, CCSD(T))
- Møller-Plesset perturbation theory (MP2)
- Time-Dependent DFT
- Quantum Monte Carlo (QMC)
- Relativistic methods (Douglas-Kroll, AMFI)
- Cholesky decomposition techniques
- Atomic Natural Orbital (ANO) basis sets

## Capabilities (CRITICAL)
- Multiconfigurational ground states (CASSCF, RASSCF)
- Excited states (MS-CASPT2, MRCI)
- Potential energy surfaces
- Conical intersection optimization
- Surface hopping dynamics (SHARC interface)
- Non-adiabatic coupling elements
- Spin-orbit coupling
- Geometry optimization (ground and excited states)
- Vibrational frequencies
- Transition state searches
- Reaction path calculations (IRC, MEP)
- Spectroscopic properties
- Transition dipole moments and oscillator strengths
- Absorption and emission spectra
- Phosphorescence and fluorescence
- Photochemical reaction mechanisms
- Chromophore modeling
- Relativistic effects
- Solvent effects (PCM)
- Large active spaces via DMRG

**Sources**: Official OpenMolcas documentation (https://gitlab.com/Molcas/OpenMolcas), confirmed in 6/7 source lists

## Key Strengths

### Multiconfigurational Methods:
- State-of-the-art CASSCF/CASPT2 implementation
- RASSCF for extended active spaces
- DMRG for very large active spaces
- Efficient algorithms and approximations
- Cholesky decomposition for integrals

### Excited States:
- Multi-state treatments (MS-CASPT2, XMS-CASPT2)
- Correct state mixing and interactions
- Conical intersections
- Non-adiabatic couplings
- Photochemistry and photophysics

### Photochemistry:
- Potential energy surface exploration
- Excited state optimization
- Surface hopping dynamics
- Radiationless transitions
- Photoisomerization mechanisms

### Specialized Features:
- ANO-type basis sets optimized for correlations
- Cholesky decomposition for efficiency
- Relativistic methods (SO-CASPT2)
- Protein chromophore embedding (QM/MM)

## Inputs & Outputs
- **Input formats**:
  - OpenMolcas input file (.input)
  - Molecule specification (XYZ, internal)
  - Keyword-based commands
  - Module-based structure
  
- **Output data types**:
  - Energies (ground and excited states)
  - Wavefunctions and CI coefficients
  - Molecular orbitals
  - Spectroscopic properties
  - Gradients and Hessians
  - Non-adiabatic couplings
  - Formatted output files

## Interfaces & Ecosystem
- **Dynamics interfaces**:
  - SHARC for surface hopping
  - COBRAMM for QM/MM dynamics
  - Newton-X interface
  
- **QM/MM**:
  - COBRAMM for proteins
  - Electrostatic embedding
  - Polarizable embedding
  
- **Visualization**:
  - Luscus (dedicated GUI)
  - Molden format support
  - Export to standard formats
  
- **Analysis tools**:
  - LOPROP for properties
  - Natural orbitals analysis
  - Density analysis

## Workflow and Usage

### Module-Based Structure:
OpenMolcas uses a modular approach:

```
&GATEWAY
 coord=water.xyz
 basis=ANO-L-VDZP

&SEWARD

&SCF

&RASSCF
 Spin=1
 nActEl=8 0 0
 Inactive=3
 Ras2=6

&CASPT2
 Multistate=3 1 2 3
```

### Common Modules:
- **GATEWAY**: Molecular geometry and basis sets
- **SEWARD**: Integral evaluation
- **SCF**: Hartree-Fock calculations
- **RASSCF/CASSCF**: Multiconfigurational SCF
- **CASPT2**: Perturbation theory
- **ALASKA**: Gradients
- **SLAPAF**: Geometry optimization
- **MCLR**: Response properties

## Advanced Features

### Active Space Selection:
- Flexible definition of active space
- Orbital selection strategies
- Natural orbital analysis
- Automated active space selection tools

### Multi-State Treatments:
- MS-CASPT2 for state averaging
- XMS-CASPT2 for extended multi-state
- State interaction via effective Hamiltonian
- Correct treatment of near-degeneracies

### DMRG:
- Large active spaces (>20 orbitals)
- Tensor network methods
- Reduced computational cost
- Maintained accuracy

### Conical Intersections:
- Optimization of CI geometries
- Branching plane analysis
- Minimum energy crossing points
- Sloped and peaked topologies

### Cholesky Decomposition:
- Reduced storage of two-electron integrals
- Faster integral evaluation
- Controllable accuracy
- Enables larger systems

## Performance Characteristics
- **Efficiency**: Optimized for multiconfigurational calculations
- **Scaling**: Active space size critical
- **Parallelization**: MPI and OpenMP support
- **Memory**: Active space dominates memory
- **Typical systems**: 10-100 atoms depending on active space

## Limitations & Known Constraints
- **Learning curve**: Very steep; requires understanding of CASSCF/CASPT2
- **Active space**: Manual selection challenging
- **Computational cost**: Multiconfigurational methods expensive
- **Documentation**: Comprehensive but technical
- **Input format**: Module-based requires understanding
- **Basis sets**: ANO basis sets recommended but not mandatory
- **Platform**: Linux, macOS
- **Not for routine DFT**: Specialized for multi-reference

## Comparison with Other Codes
- **vs Gaussian**: OpenMolcas superior for multi-reference
- **vs ORCA**: OpenMolcas more specialized CASPT2
- **vs Q-Chem**: OpenMolcas stronger in photochemistry
- **vs Molpro**: Similar capabilities, different implementations
- **Unique strength**: CASPT2 variants, photochemistry, open-source

## Application Areas

### Photochemistry:
- Photoisomerization reactions
- Photoexcitation dynamics
- Fluorescence and phosphorescence
- Vision chemistry (retinal)

### Biophysics:
- Protein chromophores
- Photoreceptor proteins
- Light-harvesting complexes
- Enzyme photocatalysis

### Transition Metal Chemistry:
- d-d transitions
- Charge-transfer states
- Magnetic properties
- Catalytic mechanisms

### Spectroscopy:
- UV-Vis absorption spectra
- Emission spectra
- Circular dichroism
- Two-photon absorption

### Method Development:
- New multiconfigurational methods
- Embedding techniques
- Dynamics algorithms

## Best Practices

### Active Space Selection:
- Include orbitals with significant changes
- Start small and expand systematically
- Check natural orbital occupations
- Validate with larger active spaces

### CASPT2 Calculations:
- Use MS-CASPT2 or XMS-CASPT2 for excited states
- Include sufficient states in averaging
- Check IPEA shift dependence
- Use imaginary level shift if needed

### Basis Sets:
- ANO-type basis sets recommended
- At least double-zeta quality
- Include diffuse functions for Rydberg states
- Correlation-consistent also acceptable

### Convergence:
- Tight SCF convergence before RASSCF
- Multiple initial orbital guesses
- Check CI convergence
- Monitor state characters

## Community and Development
- Open-source on GitLab
- International development consortium
- Regular releases
- Active user community
- Workshops and schools

## Verification & Sources
**Primary sources**:
1. Official website: https://www.molcas.org/
2. Documentation: https://molcas.gitlab.io/OpenMolcas/sphinx/
3. GitLab repository: https://gitlab.com/Molcas/OpenMolcas
4. I. Fdez. Galván et al., J. Chem. Theory Comput. 15, 5925 (2019) - OpenMolcas
5. F. Aquilante et al., J. Chem. Phys. 152, 214117 (2020) - Modern methods

**Secondary sources**:
1. OpenMolcas manual and tutorials
2. Published photochemistry applications
3. CASPT2 methodology papers
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab, LGPL v2.1)
- Community support: Mailing list, workshops
- Academic citations: >2,000 (Molcas/OpenMolcas combined)
- Active development: Regular updates
- Benchmark validation: Gold standard for multiconfigurational methods
- Specialized strength: Photochemistry, excited states, multi-reference systems
