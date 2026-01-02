# BAGEL

## Official Resources
- Homepage: https://nubakery.org/
- Documentation: https://nubakery.org/manual.html
- Source Repository: https://github.com/nubakery/bagel
- License: GNU General Public License v3.0

## Overview
BAGEL (Brilliantly Advanced General Electronic-structure Library) is a modern quantum chemistry package specializing in relativistic quantum chemistry and methods for excited states. It features state-of-the-art implementations of multireference methods, spin-orbit coupling, and analytical gradients, with emphasis on code efficiency and modern programming practices.

**Scientific domain**: Relativistic quantum chemistry, excited states, multireference methods  
**Target user community**: Researchers studying heavy elements, excited states, and systems requiring relativistic or multireference treatments

## Theoretical Methods
- Complete Active Space SCF (CASSCF)
- Multi-State CASPT2 (MS-CASPT2, XMS-CASPT2)
- CASPT2 with DMRG reference
- Multi-Reference CI (MRCI)
- Dirac-Hartree-Fock (relativistic)
- Dirac-Coulomb CASSCF
- Relativistic coupled cluster (Dirac-CCSD(T))
- Equation-of-Motion Coupled Cluster (EOM-CC)
- Hartree-Fock and DFT
- Time-Dependent DFT
- MP2 and RI-MP2
- Spin-orbit coupling (CASSCF, CASPT2)
- Relativistic gradients (analytical)
- Nuclear Energy Gradients (NEG)

## Capabilities (CRITICAL)
- Relativistic CASSCF (Dirac-Coulomb)
- Spin-orbit coupled states
- Relativistic excited states
- Analytical gradients (HF, CASSCF, relativistic)
- Geometry optimization (ground and excited states)
- Conical intersections
- Transition state searches
- Multi-state treatments
- Heavy element chemistry
- Actinide and lanthanide systems
- Transition metal complexes
- Spectroscopic properties
- Potential energy surfaces
- Non-adiabatic coupling elements
- Large active spaces via DMRG
- Efficient parallel implementation
- Modern C++ codebase

**Sources**: Official BAGEL documentation (https://nubakery.org/), confirmed in 6/7 source lists

## Key Strengths

### Relativistic Methods:
- Full four-component Dirac-Coulomb
- Dirac-CASSCF for heavy elements
- Relativistic CASPT2
- Spin-orbit CASSCF and CASPT2
- Analytical relativistic gradients

### Analytical Gradients:
- Efficient gradient implementations
- Gradients for CASSCF
- Gradients for relativistic methods
- Enables geometry optimization
- Transition state searches

### Modern Code Design:
- C++14 implementation
- Object-oriented architecture
- Template metaprogramming
- Efficient memory management
- Good parallel scaling

### Multireference Methods:
- State-of-the-art CASPT2
- DMRG-CASPT2 for large systems
- Multi-state treatments
- Spin-orbit coupling

## Inputs & Outputs
- **Input formats**:
  - JSON-based input file
  - Modern, readable syntax
  - Structured format
  - Clear hierarchical organization
  
- **Output data types**:
  - Energies (ground and excited)
  - Gradients and Hessians
  - Molecular orbitals (spinors)
  - Spectroscopic properties
  - Spin-orbit coupling elements
  - Formatted output

## Interfaces & Ecosystem
- **Parallelization**:
  - MPI for distributed memory
  - OpenMP for shared memory
  - Hybrid MPI+OpenMP
  - Efficient task distribution
  
- **Libraries**:
  - Uses modern linear algebra libraries
  - Optimized BLAS/LAPACK
  - Efficient integral evaluation
  
- **Visualization**:
  - Export to standard formats
  - Compatible with visualization tools

## Workflow and Usage

### JSON Input Format:
BAGEL uses a modern JSON-based input:

```json
{
  "bagel" : [
    {
      "title" : "molecule",
      "basis" : "svp",
      "df_basis" : "svp-jkfit",
      "angstrom" : true,
      "geometry" : [
        {"atom" : "O", "xyz" : [0.0, 0.0, 0.0]},
        {"atom" : "H", "xyz" : [0.0, 0.0, 1.0]},
        {"atom" : "H", "xyz" : [0.0, 1.0, 0.0]}
      ]
    },
    {
      "title" : "hf"
    },
    {
      "title" : "casscf",
      "nstate" : 3,
      "nclosed" : 3,
      "nact" : 6,
      "state" : [1, 1, 1]
    }
  ]
}
```

### Common Calculation Types:
- **HF**: Hartree-Fock
- **CASSCF**: Multiconfigurational SCF
- **CASPT2**: Perturbation theory
- **DHF**: Dirac-Hartree-Fock
- **FORCE**: Gradient calculation

## Advanced Features

### Spin-Orbit CASSCF:
- Variational spin-orbit coupling
- State interaction approach
- Accurate for heavy elements
- Analytical gradients available

### DMRG Integration:
- Large active spaces
- DMRG-CASSCF
- DMRG-CASPT2
- Efficient tensor network methods

### Relativistic Gradients:
- Analytical gradients for Dirac-Coulomb
- Geometry optimization of heavy systems
- Transition states with relativity
- Efficient implementations

### Multi-State Methods:
- MS-CASPT2 and XMS-CASPT2
- State averaging
- Correct state interactions
- Spin-orbit multi-state

## Performance Characteristics
- **Efficiency**: Modern implementation, very efficient
- **Scaling**: Good parallel scaling
- **Parallelization**: Excellent MPI and OpenMP
- **Memory**: Efficient memory management
- **Typical systems**: 10-100 atoms depending on method

## Computational Comparison
- **Gradients**: BAGEL often faster than competitors
- **Relativistic**: Comparable to DIRAC
- **CASPT2**: Competitive with OpenMolcas
- **Parallel**: Excellent scaling

## Limitations & Known Constraints
- **Learning curve**: Steep; requires quantum chemistry expertise
- **JSON format**: Different from traditional codes
- **Documentation**: Good but assumes background
- **Active space**: Manual selection needed
- **Compilation**: Requires modern C++ compiler
- **Platform**: Linux, macOS (requires compilation)
- **Specialized**: Not for routine DFT

## Comparison with Other Codes
- **vs DIRAC**: BAGEL has analytical gradients
- **vs OpenMolcas**: BAGEL more modern code, similar methods
- **vs ORCA**: BAGEL specialized relativity and gradients
- **vs Molpro**: Similar capabilities, BAGEL open-source
- **Unique strength**: Relativistic gradients, modern code, efficiency

## Application Areas

### Heavy Element Chemistry:
- Actinide complexes
- Lanthanide spectroscopy
- Heavy transition metals
- Superheavy elements

### Photochemistry:
- Excited state dynamics
- Conical intersections with relativity
- Spin-orbit effects on photochemistry
- Heavy atom photocatalysis

### Transition Metal Catalysis:
- Mechanism elucidation
- Spin-orbit effects
- Geometry optimization
- Reaction barriers

### Spectroscopy:
- Spin-orbit split spectra
- Relativistic effects on spectra
- Heavy element NMR/EPR
- Absorption and emission

## Best Practices

### Input Preparation:
- Use JSON validators
- Check input syntax carefully
- Start with small test cases
- Understand hierarchical structure

### Method Selection:
- Dirac-Coulomb for heavy elements
- CASSCF for multireference
- CASPT2 for dynamic correlation
- Gradients for optimization

### Active Space:
- Select active space carefully
- Use natural orbitals
- Start small and expand
- Validate results

### Convergence:
- Tight convergence for gradients
- Check SCF and CASSCF convergence
- Monitor state characters
- Test different initial guesses

## Community and Development
- Open-source on GitHub
- Active development by Shiozaki group
- Regular updates
- Modern software engineering practices
- Issue tracking on GitHub

## Educational Resources
- Official manual
- Example inputs
- Publication list
- Tutorial materials

## Verification & Sources
**Primary sources**:
1. Official website: https://nubakery.org/
2. Documentation: https://nubakery.org/manual.html
3. GitHub repository: https://github.com/nubakery/bagel
4. T. Shiozaki, WIREs Comput. Mol. Sci. 8, e1331 (2018) - BAGEL overview
5. M. K. MacLeod and T. Shiozaki, J. Chem. Phys. 142, 051103 (2015) - Gradients

**Secondary sources**:
1. BAGEL manual and examples
2. Published applications
3. Methodology papers
4. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: CONFIRMED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Documentation: COMPREHENSIVE and ACCESSIBLE
- Source code: OPEN (GitHub, GPL v3)
- Community support: GitHub issues, mailing list
- Academic citations: >150
- Active development: Regular updates
- Modern codebase: C++14, efficient implementation
- Specialized strength: Relativistic analytical gradients, modern design
