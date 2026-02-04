# IrRep

## Official Resources
- Homepage: https://irrep.readthedocs.io/
- Documentation: https://irrep.readthedocs.io/en/latest/
- Source Repository: https://github.com/stepan-tsirkin/irrep
- License: GNU General Public License v3.0

## Overview
IrRep is a Python code for calculating the irreducible representations of Bloch states in ab-initio calculations. It interfaces with Quantum ESPRESSO, VASP, and Abinit to determine the symmetry properties of electronic bands, which is essential for identifying topological invariants, enforcing selection rules, and understanding band connectivity.

**Scientific domain**: Symmetry analysis, irreducible representations, topological materials  
**Target user community**: Topological physics researchers, DFT users

## Theoretical Methods
- Group theory
- Irreducible representations of space groups
- Trace of symmetry operators
- Character tables
- Topological quantum chemistry (TQC)
- Elementary Band Representations (EBR)

## Capabilities (CRITICAL)
- Calculation of irreps for all space groups (with/without SOC)
- Identification of topological invariants
- Analysis of band connectivity
- Interface with Quantum ESPRESSO, VASP, Abinit
- Calculation of symmetry eigenvalues
- Trace of symmetry operators

**Sources**: IrRep documentation, Comp. Phys. Comm. 262, 107836 (2021)

## Inputs & Outputs
- **Input formats**: DFT wavefunctions (WAVECAR/tmp.pp), structural data
- **Output data types**: Irrep labels, traces, character tables

## Interfaces & Ecosystem
- **Quantum ESPRESSO**: Native interface
- **VASP**: Supported via WAVECAR
- **Abinit**: Supported
- **Python**: Core implementation language

## Workflow and Usage
1. Perform DFT calculation.
2. Run IrRep: `python -m irrep <input_flags>`
3. Analyze output for symmetry labels.

## Performance Characteristics
- Fast analysis of wavefunctions
- Efficient handling of large basis sets

## Application Areas
- Topological insulators
- Weyl semimetals
- Symmetry-protected topological phases
- Band structure analysis

## Community and Support
- Open-source (GPL v3)
- GitHub repository
- Developed by Stepan Tsirkin (University of Zurich)

## Verification & Sources
**Primary sources**:
1. Homepage: https://irrep.readthedocs.io/
2. GitHub: https://github.com/stepan-tsirkin/irrep
3. Publication: M. Iraola, J. L. Mañes, B. Bradlyn, T. Neupert, M. G. Vergniory, S. S. Tsirkin, Comp. Phys. Comm. 262, 107836 (2021)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE (Tsirkin)
- Applications: Irreps, topological phases, DFT interface
