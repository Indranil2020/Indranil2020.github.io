# MagneticKP

## Official Resources
- Homepage: https://github.com/zhangzeyingvv/MagneticKP
- Documentation: https://github.com/zhangzeyingvv/MagneticKP (README)
- Source Repository: https://github.com/zhangzeyingvv/MagneticKP
- License: GNU General Public License v3.0

## Overview
MagneticKP is a software package (available in Mathematica and Python) for efficiently constructing k·p effective Hamiltonians for magnetic and non-magnetic crystals. Developed by Zhang et al., it implements a novel "Iterative Simplification Algorithm" (ISA) to solve the symmetry constraints imposed by Magnetic Space Groups (MSGs), allowing for the rapid generation of k·p models expanded to arbitrary order in k.

**Scientific domain**: k·p perturbation theory, Magnetic Space Groups, Effective Models
**Target user community**: Condensed matter theorists, topological materials researchers

## Theoretical Methods
- k·p Perturbation Theory
- Magnetic Space Group (MSG) Symmetry
- Method of Invariants
- Iterative Simplification Algorithm (ISA)
- Direct-Product Decomposition Algorithm (DDA)

## Capabilities (CRITICAL)
**Category**: Magnetic k·p model generator
- **Model Construction**: Generates symmetry-allowed k·p Hamiltonians
- **Symmetry**: Supports all 1651 Magnetic Space Groups
- **Expansion**: Arbitrary order in wavenumber k
- **Basis**: Spinless (single-valued) and Spinful (double-valued) representations
- **Performance**: High efficiency using ISA algorithm
- **Versions**:
  - **Mathematica**: Full feature set, symbolic manipulation
  - **Python**: Core functionality, integration with Python stack

**Sources**: GitHub repository, Publication (Comput. Phys. Comm. 283, 108575 (2023))

## Key Strengths

### Algorithmic Efficiency:
- ISA algorithm is significantly faster than traditional projection methods
- Handles high-order expansions and high-dimensional representations

### Magnetic Symmetry:
- Native support for magnetic symmetries
- Essential for magnetic topological insulators and semimetals

### Dual Implementation:
- Mathematica for symbolic derivation
- Python for numerical integration/workflows

## Status
- **Type**: Mathematica/Python Package
- **Development**: Active
- **Maintainer**: Zeying Zhang (zhangzeyingvv)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/zhangzeyingvv/MagneticKP
2. Publication: Zhang, Z., et al. "MagneticKP: A package for quickly constructing k·p models..." Comput. Phys. Comm. 283, 108575 (2023).

**Confidence**: VERIFIED - Research Tool

**Verification status**: ✅ CONFIRMED
- GitHub: ACCESSIBLE
- **Note**: Replaces incorrect reference to `andrewfeng12`. Companion tool to MagneticTB.
