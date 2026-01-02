# Turbomole

## Official Resources
- Homepage: https://www.turbomole.org/
- Documentation: https://www.turbomole.org/turbomole/documentation/
- Source Repository: Commercial (source available to licensees)
- License: Commercial with academic licenses available

## Overview
Turbomole is a highly efficient quantum chemistry program package for ab initio electronic structure calculations. Developed at the University of Karlsruhe, it is renowned for exceptional computational efficiency, particularly for large molecules, and features advanced algorithms including RI (Resolution of Identity) approximations for accelerated calculations.

**Scientific domain**: Quantum chemistry, electronic structure, molecular properties, spectroscopy  
**Target user community**: Computational chemists studying molecular systems, reactions, and properties

## Theoretical Methods
- Hartree-Fock (HF)
- Density Functional Theory (DFT) with many functionals
- Møller-Plesset perturbation theory (MP2, RI-MP2)
- Coupled Cluster (CC2, CCSD, CCSD(T), ADC(2))
- Time-Dependent DFT (TDDFT)
- Multi-reference methods (MRCI, CASSCF)
- Semiempirical methods (DFTB)
- Resolution of Identity (RI) approximations
- Multipole Accelerated RI (marij)
- Conductor-like Screening Model (COSMO) solvation
- Excited state calculations
- Ground and excited state gradients
- Numerical frequencies and vibrational analysis

## Capabilities (CRITICAL)
- Ground state energies (HF, DFT, post-HF)
- Geometry optimization (efficient algorithms)
- Transition state searches
- Reaction path following (IRC)
- Vibrational frequencies and IR/Raman spectra
- Thermochemistry
- Excited states (TDDFT, CC2, ADC(2))
- Absorption and emission spectra
- Circular dichroism (CD) spectra
- Magnetic properties (NMR, EPR)
- Electric properties (polarizabilities, hyperpolarizabilities)
- Solvation effects (COSMO, COSMO-RS)
- Periodic boundary conditions (DFT only)
- Dispersion corrections (DFT-D3, DFT-D4)
- Relativistic effects (DKH, X2C)
- Systems up to several hundred atoms
- Exceptional efficiency for large molecules

**Sources**: Official Turbomole documentation (https://www.turbomole.org/), confirmed in multiple source lists

## Key Strengths

### Computational Efficiency:
- RI approximation dramatically speeds DFT and MP2
- marij (multipole accelerated RI) for exchange integrals
- Efficient algorithms for large molecules (100-1000 atoms)
- Low memory requirements
- Linear scaling for selected methods
- Efficient parallelization

### Resolution of Identity (RI):
- RI-DFT: Faster DFT calculations
- RI-MP2: MP2 with near-linear scaling
- RI-CC2: Efficient coupled cluster
- Minimal loss of accuracy
- Standard auxiliary basis sets available

### User Interface:
- Command-line tools for flexible scripting
- define utility for input generation
- Module-based architecture
- Well-designed workflow
- Extensive default settings

## Inputs & Outputs
- **Input formats**:
  - coord (Turbomole coordinate format)
  - control (main input file)
  - basis (basis set specifications)
  - Generated via define utility or manually
  
- **Output data types**:
  - Energies and gradients
  - Optimized geometries
  - Vibrational frequencies
  - Molecular orbitals
  - Spectroscopic properties
  - Population analyses
  - Thermodynamic data

## Interfaces & Ecosystem
- **GUIs available**:
  - TmoleX (official graphical interface)
  - Integration with other GUIs (ParaTools, etc.)
  
- **Workflow tools**:
  - Python scripting interfaces
  - Integration with workflow managers
  - Batch job submission tools
  
- **Visualization**:
  - Export to standard formats
  - Compatible with VMD, Molden, etc.
  - TmoleX integrated visualization

## Module Structure

Turbomole uses a modular design with specialized programs:

### Ground State:
- **dscf**: DFT and HF SCF calculations
- **ridft**: RI-DFT (faster DFT)
- **ricc2**: RI-CC2 and related methods
- **mpgrad**: MP2 gradients

### Geometry Optimization:
- **jobex**: Job execution and optimization driver
- **statpt**: Stationary point optimization
- **relax**: Geometry relaxation

### Excited States:
- **escf**: TDDFT excitation energies
- **egrad**: Excited state gradients
- **ricc2**: CC2, ADC(2) excited states

### Properties:
- **aoforce**: Analytical second derivatives
- **NumForce**: Numerical frequencies
- **mpshift**: NMR shielding tensors
- **escf**: Optical properties

## Workflow and Usage

### Typical Workflow:

1. **Setup**:
   ```bash
   # Generate input with define
   define
   # Or use TmoleX GUI
   ```

2. **Single Point**:
   ```bash
   # DFT calculation
   ridft > ridft.out
   ```

3. **Optimization**:
   ```bash
   # Automatic optimization
   jobex -ri -c 200
   ```

4. **Frequencies**:
   ```bash
   # Analytical frequencies
   aoforce > aoforce.out
   ```

5. **Excited States**:
   ```bash
   # TDDFT excitations
   escf > escf.out
   ```

### Define Utility:
- Interactive input generation
- Molecule definition
- Basis set selection
- Method specification
- Reasonable defaults
- Saves to control file

## Advanced Features

### RI Approximations:
- RI-J: Coulomb integrals
- RI-K (marij): Exchange integrals
- RI-MP2: Correlated calculations
- Auxiliary basis sets optimized
- Near-linear scaling achieved

### Excited State Dynamics:
- Excited state gradients
- Conical intersection optimizations
- Non-adiabatic coupling elements
- Surface hopping capabilities (with extensions)

### Solvation Models:
- COSMO: Implicit solvation
- COSMO-RS: Realistic solvation free energies
- Wide range of solvents
- Parameterized for various properties

### Dispersion Corrections:
- DFT-D3: Grimme's D3 correction
- DFT-D4: Latest dispersion model
- Essential for non-covalent interactions
- Automatic parameter selection

## Performance Characteristics
- **Efficiency**: Among the fastest quantum chemistry codes
- **Scaling**: Near-linear for RI methods
- **Parallelization**: SMP and MPI
- **Memory**: Low memory footprint
- **Typical systems**:
  - Small molecules: seconds to minutes
  - 100 atoms: minutes to hours (DFT)
  - 500 atoms: hours to days (RI-DFT)

## Limitations & Known Constraints
- **Commercial software**: License required (academic discounts available)
- **Learning curve**: Moderate; command-line interface
- **Documentation**: Comprehensive but technical
- **GUIs**: TmoleX recommended for ease of use
- **Periodic systems**: Limited to DFT
- **Very large systems**: Not as specialized as plane-wave codes
- **Some methods**: Not all quantum chemistry methods available
- **Platform**: Linux, macOS, Windows

## Comparison with Other Codes
- **vs Gaussian**: Turbomole often faster, especially RI methods
- **vs ORCA**: Similar capabilities, different strengths
- **vs Q-Chem**: Comparable, Turbomole very efficient
- **vs PSI4**: Turbomole commercial but faster
- **Unique strength**: RI efficiency, large molecule handling

## Application Areas

### Organic Chemistry:
- Reaction mechanisms
- Conformational analysis
- Thermochemistry
- Spectroscopy predictions

### Photochemistry:
- Excited state properties
- Absorption and emission spectra
- Photochemical pathways
- Fluorescence/phosphorescence

### Materials Chemistry:
- Molecular materials
- Supramolecular systems
- Host-guest complexes
- Organic semiconductors

### Drug Design:
- Molecular properties
- Solvation effects
- Binding energies
- QSAR studies

## Best Practices

### Input Preparation:
- Use define for consistent setups
- Choose appropriate basis sets
- Consider RI methods for speed
- Include dispersion for non-covalent

### Convergence:
- Default settings usually good
- Monitor SCF convergence
- Check geometry optimization
- Validate frequencies (no imaginary for minima)

### Efficiency:
- Use RI methods when possible
- Appropriate basis set size
- Parallelize effectively
- Start with smaller basis for testing

## Licensing and Access
- **Academic licenses**: Available at reduced cost
- **Commercial licenses**: Full price for industry
- **Support**: Professional support included
- **Updates**: Regular releases and bug fixes
- **Training**: Workshops and tutorials available

## Verification & Sources
**Primary sources**:
1. Official website: https://www.turbomole.org/
2. Documentation: https://www.turbomole.org/turbomole/documentation/
3. R. Ahlrichs et al., Chem. Phys. Lett. 162, 165 (1989) - Turbomole program
4. F. Weigend and M. Häser, Theor. Chem. Acc. 97, 331 (1997) - RI-MP2
5. C. Hättig and F. Weigend, J. Chem. Phys. 113, 5154 (2000) - RI-CC2

**Secondary sources**:
1. Turbomole user manual
2. Published applications across chemistry
3. Benchmark studies and comparisons
4. Confirmed in multiple source lists

**Confidence**: CONFIRMED - Well-established commercial code

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: Available to licensees
- Community support: Professional support, user forums
- Academic citations: >5,000
- Active development: Regular releases
- Benchmark validation: Extensive validation published
- Wide adoption: Standard tool in computational chemistry
