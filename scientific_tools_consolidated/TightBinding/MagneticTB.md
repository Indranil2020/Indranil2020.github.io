# MagneticTB

## Official Resources
- Homepage: https://github.com/zhangzeyingvv/MagneticTB
- Documentation: https://github.com/zhangzeyingvv/MagneticTB (README & Examples)
- Source Repository: https://github.com/zhangzeyingvv/MagneticTB
- License: GNU General Public License v3.0

## Overview
MagneticTB is a Mathematica package designed for the automated construction of tight-binding models for magnetic and non-magnetic materials involving arbitrary magnetic space groups (MSGs). Developed by Zhang et al., it allows users to generate symmetry-allowed tight-binding Hamiltonians by simply inputting the magnetic space group number and orbital information. It significantly simplifies the modeling of magnetic topological materials.

**Scientific domain**: Magnetic Space Groups, Tight-Binding Models, Topological Magnetism
**Target user community**: Researchers in topological magnetism and spintronics

## Theoretical Methods
- Magnetic Space Group (MSG) Theory
- Symmetry-Adapted Tight-Binding Models
- Corepresentations of MSGs
- Group Theory (Wyckoff positions, Irreps)

## Capabilities (CRITICAL)
**Category**: Mathematica magnetic TB package
- **Model Construction**: Automatically generates symmetry-allowed Hamiltonian matrices
- **Magnetic Symmetry**: Supports all 1651 Magnetic Space Groups
- **Input**:
  - MSG Number (BNS or OG setting)
  - Wyckoff positions of atoms
  - Orbital basis (s, p, d, f)
- **Output**: Analytical or numerical Hamiltonian parameters
- **Analysis**:
  - Symmetry operators in matrix form
  - Band structure calculation
  - Interface with other codes (e.g., Z2Pack, WannierTools)

**Sources**: GitHub repository, Publication (Comput. Phys. Comm. 270, 108153 (2022))

## Key Strengths

### Automated Symmetry Analysis:
- Eliminates manual derivation of symmetry constraints
- Handles complex magnetic symmetries automatically
- Ensures Hamiltonians respect all MSG symmetries

### Versatility:
- Applicable to both magnetic and non-magnetic crystals
- Handles spin-orbit coupling (SOC) naturally

### User-Friendly:
- High-level Mathematica interface
- Minimal input requirements

## Status
- **Type**: Mathematica Package
- **Development**: Active
- **Maintainer**: Zeying Zhang (zhangzeyingvv)
- **Language**: Mathematica

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/zhangzeyingvv/MagneticTB
2. Publication: Zhang, Z., et al. "MagneticTB: A package for tight-binding model of magnetic and non-magnetic materials." Comput. Phys. Comm. 270, 108153 (2022).

**Confidence**: VERIFIED - Specialized Research Tool

**Verification status**: ✅ CONFIRMED
- GitHub: ACCESSIBLE
- **Note**: Replaces incorrect reference to `andrewfeng12` or `xfzhang-phys`. The correct developer is `zhangzeyingvv`. Requires Mathematica.
