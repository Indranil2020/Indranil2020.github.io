# Dalton

## Official Resources
- Homepage: https://www.daltonprogram.org/
- Documentation: https://daltonprogram.org/documentation/
- Source Repository: https://gitlab.com/dalton/dalton
- License: GNU Lesser General Public License v2.1

## Overview
Dalton is a powerful quantum chemistry program with particular emphasis on molecular properties, response theory, and environment effects. Originally developed in Scandinavia, it excels at computing spectroscopic properties, electric and magnetic molecular properties, and various response properties using advanced wave function and DFT methods.

**Scientific domain**: Quantum chemistry, molecular properties, spectroscopy, response theory  
**Target user community**: Researchers studying molecular properties, spectroscopy, and environment effects

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT)
- Multi-Configurational Self-Consistent Field (MCSCF)
- Coupled Cluster (CC2, CCSD, CC3)
- Multi-reference CI (MRCI)
- Møller-Plesset perturbation theory (MP2)
- Time-Dependent DFT (TDDFT)
- Polarizable Embedding (PE)
- QM/MM with polarizable force fields
- Response theory (linear, quadratic, cubic)
- Relativistic methods (DKH, X2C)
- Frozen Density Embedding (FDE)

## Capabilities (CRITICAL)
- Ground state energies (HF, DFT, CC)
- Geometry optimization
- Vibrational frequencies (analytical)
- Excited states (TDDFT, CC2, MCSCF)
- Linear response properties
- Quadratic and cubic response properties
- NMR shielding and spin-spin coupling
- EPR parameters (g-tensors, hyperfine coupling)
- Optical rotation and circular dichroism
- Magnetic circular dichroism (MCD)
- Two-photon absorption
- Raman and resonance Raman spectra
- Vibrational circular dichroism (VCD)
- Electric and magnetic properties
- Polarizabilities and hyperpolarizabilities
- Solvent effects (PCM, PE model)
- QM/MM calculations
- Environment effects on properties
- Relativistic property calculations

**Sources**: Official Dalton documentation (https://www.daltonprogram.org/), verified in multiple source lists

## Key Strengths

### Molecular Properties:
- Comprehensive property calculations
- Response theory implementation
- High accuracy for spectroscopic properties
- NMR parameters (chemical shifts, couplings)
- Magnetic properties (EPR, magnetizabilities)
- Optical properties (OR, CD, MCD)

### Environment Effects:
- Polarizable Embedding (PE) model
- QM/MM with advanced force fields
- Frozen Density Embedding
- PCM solvation
- Accurate environment effects on properties

### Response Theory:
- Linear response for excited states
- Quadratic response for two-photon absorption
- Cubic response for nonlinear optics
- Analytical derivatives
- Time-dependent properties

### Wave Function Methods:
- MCSCF for multi-reference systems
- Coupled Cluster for high accuracy
- CC response theory
- Multi-reference CI

## Inputs & Outputs
- **Input formats**:
  - Dalton input file (.dal)
  - Molecule file (.mol)
  - Text-based input
  - Structured keyword format
  
- **Output data types**:
  - Energies and properties
  - Molecular orbitals
  - Response properties
  - Spectroscopic parameters
  - Population analyses
  - Property tensors
  - Formatted output

## Interfaces & Ecosystem
- **QM/MM interfaces**:
  - Interface to AMBER
  - Interface to GROMACS
  - Generic MM interfaces
  - Polarizable Embedding
  
- **Visualization**:
  - Molden format support
  - Export to standard formats
  - Property visualization tools
  
- **Scripting**:
  - Python scripting possible
  - Batch processing utilities
  - Workflow integration

## Workflow and Usage

### Input Structure:
Dalton uses two input files:
1. **.dal file**: Computational details
2. **.mol file**: Molecular structure and basis

### Example Input (.dal):
```
**DALTON INPUT
.RUN WAVE FUNCTION
**WAVE FUNCTIONS
.DFT
 B3LYP
**RESPONSE
*LINEAR
.DIPLEN
**END OF DALTON INPUT
```

### Example Molecule (.mol):
```
BASIS
6-31G*
Water molecule

Atomtypes=2
Charge=8.0 Atoms=1
O   0.000000  0.000000  0.117176
Charge=1.0 Atoms=2
H   0.000000  0.757200 -0.469706
H   0.000000 -0.757200 -0.469706
```

## Advanced Features

### Polarizable Embedding:
- Include environment polarization explicitly
- QM region embedded in polarizable MM
- Accurate spectroscopic properties in solution
- Protein environment effects
- Multi-layer embedding schemes

### Response Property Calculations:
- Systematic hierarchy of response
- Linear: Excitation energies, transition moments
- Quadratic: Two-photon absorption, Raman
- Cubic: Third-order susceptibilities
- Damped response for frequency-dependent properties

### Magnetic Properties:
- NMR chemical shifts (GIAO method)
- Spin-spin coupling constants
- EPR g-tensors and hyperfine couplings
- Magnetizabilities
- Magnetic circular dichroism

### Chiroptical Properties:
- Optical rotation
- Electronic circular dichroism (ECD)
- Vibrational circular dichroism (VCD)
- Raman optical activity (ROA)
- Magnetic circular dichroism (MCD)

## Performance Characteristics
- **Efficiency**: Good for property calculations
- **Scaling**: Standard quantum chemistry scaling
- **Parallelization**: MPI support
- **Memory**: Moderate requirements
- **Typical systems**: 10-200 atoms depending on method

## Limitations & Known Constraints
- **Learning curve**: Steep; input format requires understanding
- **Documentation**: Comprehensive but technical
- **GUI**: Limited graphical interface options
- **Very large systems**: Not optimized for huge molecules
- **Some methods**: Specialized; not all QC methods available
- **Platform**: Primarily Linux/Unix
- **Input format**: Two-file system requires care

## Comparison with Other Codes
- **vs Gaussian**: Dalton stronger in response properties
- **vs ORCA**: Dalton more specialized for properties
- **vs Q-Chem**: Dalton unique PE and FDE methods
- **vs ADF**: Similar emphasis on properties
- **Unique strength**: Response theory, environment embedding, magnetic properties

## Application Areas

### Spectroscopy:
- NMR spectroscopy predictions
- UV-Vis absorption spectra
- Circular dichroism studies
- Raman and ROA spectroscopy
- EPR parameter calculations

### Chirality:
- Optical rotation
- Absolute configuration determination
- Chiroptical spectroscopy
- Chiral recognition

### Biochemistry:
- Protein chromophores
- Environment effects on spectra
- QM/MM studies of enzymes
- Protein-ligand interactions

### Materials:
- Nonlinear optical properties
- Molecular response properties
- Polarizabilities
- Hyperpolarizabilities

## Best Practices

### Input Preparation:
- Understand keyword hierarchy
- Use appropriate basis sets
- Choose method for property
- Consider symmetry

### Property Calculations:
- Converge ground state first
- Check basis set requirements
- Validate with experiments
- Consider environment effects

### Environment Modeling:
- Use PE for biomolecules
- QM/MM for extended systems
- Validate embedding region
- Test convergence with embedding

### Spectroscopy:
- Include solvent effects
- Use appropriate functionals
- Check gauge origin for magnetic
- Consider relativistic effects

## Community and Support
- Open-source on GitLab
- Academic development
- User mailing list
- Regular workshops
- Active Nordic community

## Verification & Sources
**Primary sources**:
1. Official website: https://www.daltonprogram.org/
2. Documentation: https://daltonprogram.org/documentation/
3. GitLab repository: https://gitlab.com/dalton/dalton
4. K. Aidas et al., WIREs Comput. Mol. Sci. 4, 269 (2014) - Dalton review
5. T. Helgaker et al., Chem. Rev. 112, 543 (2012) - Response theory

**Secondary sources**:
1. Dalton manual and tutorials
2. Published applications in spectroscopy
3. PE and QM/MM method papers
4. Verified in multiple source lists

**Confidence**: VERIFIED - Well-established in spectroscopy community

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitLab, LGPL v2.1)
- Community support: Mailing list, workshops
- Academic citations: >1,500
- Active development: Regular updates
- Benchmark validation: Extensive spectroscopy benchmarks
- Specialized strength: Molecular properties and environment effects
