# NBO

## Official Resources
- Homepage: https://nbo7.chem.wisc.edu/
- Alternate Homepage: https://nbo.chem.wisc.edu/
- Publication: F. Weinhold et al., J. Comput. Chem. 40, 1497 (2019)
- Related Utility: NBOPro@Jmol via official NBO pages

## Overview
NBO is the Natural Bond Orbital program for analyzing localized and delocalized chemical bonding in wavefunctions. It provides a broad set of orbital, population, donor-acceptor, and bond-character descriptors that are widely used across computational chemistry for chemically intuitive interpretation of molecular electronic structure.

**Scientific domain**: Orbital-based chemical bonding analysis, population analysis, donor-acceptor interactions  
**Target user community**: Quantum chemists and electronic-structure users seeking chemically intuitive orbital and bonding descriptors

## Theoretical Methods
- Natural Atomic Orbitals (NAOs)
- Natural Hybrid Orbitals (NHOs)
- Natural Bond Orbitals (NBOs)
- Natural Localized Molecular Orbitals (NLMOs)
- Natural Population Analysis (NPA)
- Second-order donor-acceptor perturbation analysis

## Capabilities (CRITICAL)
- Localized orbital description of chemical bonding
- Natural population analysis and atomic charges
- Donor-acceptor interaction analysis via second-order perturbation theory
- Bonding, antibonding, lone-pair, and Rydberg orbital characterization
- Broad integration with quantum chemistry workflows and educational use
- Widely cited standard tool for orbital-based bonding analysis

**Sources**: Official NBO pages and JCC publication on NBO 7.0

## Key Strengths

### Chemical Interpretability:
- Intuitive localized bonding picture
- Strong orbital language for chemists
- Donor-acceptor analysis is widely adopted
- Useful for both routine and advanced bonding studies

### Method Breadth:
- Charges and populations
- Bond orbital analysis
- Hyperconjugation and delocalization interpretation
- Localized molecular orbital representations

### Ecosystem Maturity:
- Long-established software family
- Formal NBO 7 publication
- Additional visualization utilities such as NBOPro@Jmol

## Inputs & Outputs
- **Input formats**:
  - Quantum-chemistry wavefunction information through supported interfaces and NBO-style input workflows

- **Output data types**:
  - Natural populations and charges
  - Bonding and antibonding orbitals
  - Donor-acceptor interaction energies
  - Localized orbital descriptors and reports

## Workflow and Usage
1. Run a compatible electronic-structure calculation with NBO analysis enabled.
2. Generate NBO outputs from the wavefunction.
3. Inspect populations, orbitals, and donor-acceptor interactions.
4. Use the results to interpret bonding, delocalization, and reactivity.

## Performance Characteristics
- Efficient post-processing of wavefunction information
- Deep analysis output rather than minimal summary metrics
- Standard tool in molecular electronic-structure interpretation

## Limitations & Known Constraints
- **Orbital framework**: Provides a localized-orbital view rather than density-topology analysis
- **Licensing/distribution**: Not a simple open-source package workflow
- **Scope**: Best complemented by QTAIM or density-based tools when topological information is needed

## Comparison with Other Tools
- **vs JANPA**: NBO is broader and more established; JANPA offers an open-source NPA-centered workflow
- **vs EDDB**: NBO emphasizes localized orbitals and donor-acceptor interactions, whereas EDDB emphasizes delocalization measures
- **vs QTAIM tools**: NBO is orbital-based, not critical-point/basin-based
- **Unique strength**: The canonical localized-orbital framework for chemically intuitive bonding analysis

## Application Areas
- Bonding and antibonding analysis
- Hyperconjugation and donor-acceptor studies
- Charge and population analysis
- Interpretation of molecular electronic structure and reactivity

## Community and Support
- Long-established official project pages
- Strong literature presence
- Widely taught and used in computational chemistry

## Verification & Sources
**Primary sources**:
1. Homepage: https://nbo7.chem.wisc.edu/
2. Alternate homepage: https://nbo.chem.wisc.edu/
3. F. Weinhold et al., J. Comput. Chem. 40, 1497 (2019)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Official homepage: ACCESSIBLE
- Publication: AVAILABLE
- Community adoption: VERY STRONG
- Primary use case: Orbital-based chemical bonding and population analysis
