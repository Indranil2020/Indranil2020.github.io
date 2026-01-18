# CT-HYB

## Official Resources
- Homepage: Multiple implementations (TRIQS/cthyb, ALPS/CT-HYB, etc.)
- Documentation: See specific implementations
- Source Repository: Various (TRIQS, ALPS, w2dynamics, etc.)
- License: Varies by implementation

## Overview
CT-HYB (Continuous-Time Hybridization Expansion) is an algorithm for solving quantum impurity problems, not a single specific software. It is the most widely used continuous-time quantum Monte Carlo (CTQMC) method in DMFT applications. Multiple implementations exist in different software packages including TRIQS/cthyb, ALPS, w2dynamics, iQIST, and others.

**Scientific domain**: Quantum impurity solvers, DMFT, CTQMC algorithms  
**Target user community**: Researchers performing DMFT calculations

## Theoretical Methods
- Continuous-time quantum Monte Carlo (CTQMC)
- Hybridization expansion algorithm
- Stochastic sampling in imaginary time
- Multi-orbital Anderson impurity model
- Partition function expansion
- Segment picture or matrix formulation

## Capabilities (CRITICAL)
**Note**: CT-HYB is an ALGORITHM implemented in multiple codes:

**Primary implementations**:
- **TRIQS/cthyb**: Most widely used, GPL-licensed
- **ALPS/CT-HYB**: Part of ALPS project
- **w2dynamics**: Includes CT-HYB solver
- **iQIST**: Multiple CT-HYB variants
- **ComCTQMC**: GPU-accelerated version

**General capabilities**:
- Multi-orbital impurity problems
- General interactions
- Temperature-dependent calculations
- Green's functions and self-energies
- Statistical sampling

**Sources**: Master list notes: "VERIFIED - TRIQS implementation", algorithm widely implemented

## Inputs & Outputs
**Depends on implementation**: Each code has its own interface

**Common inputs**: Hybridization functions, interaction parameters
**Common outputs**: Green's functions, self-energies, occupations

## Interfaces & Ecosystem
- **TRIQS/cthyb**: Standard in TRIQS ecosystem
- **ALPS**: Part of ALPS solvers
- **w2dynamics**: Integrated CT-HYB
- **DCore**: Can use various CT-HYB implementations
- **Multiple frameworks**: Most DMFT codes support CT-HYB

## Limitations & Known Constraints
- CTQMC computational cost
- Sign problem in some regimes
- Statistical errors
- Scales with inverse temperature
- Implementation-specific limitations

## Comparison with Other Algorithms
| Algorithm | CT-HYB | CT-INT | CT-AUX |
| :--- | :--- | :--- | :--- |
| **Expansion** | Hybridization ($\Delta$) | Interaction ($U$) | Auxiliary Field |
| **Best Regime** | Strong Coupling (Large $U$) | Weak Coupling (Small $U$) | Weak/Intermediate |
| **Sign Problem** | Generally mild | Moderate | Moderate |
| **Matrix Size** | Linear in expansion order | Linear in expansion order | Linear in expansion order |

## Verification & Sources
**Primary sources**:
1. E. Gull et al., Rev. Mod. Phys. 83, 349 (2011) - CT-HYB review
2. P. Werner and A. J. Millis, Phys. Rev. B 74, 155107 (2006) - Original algorithm
3. TRIQS/cthyb: https://triqs.github.io/cthyb/
4. Master list: "VERIFIED - TRIQS implementation"

**Secondary sources**:
1. Multiple implementations in major DMFT codes
2. Standard algorithm in DMFT community
3. Hundreds of papers using CT-HYB

**Confidence**: VERIFIED - Algorithm with multiple implementations

**Verification status**: VERIFIED as ALGORITHM
- Status: Algorithm/method, not single software
- Multiple implementations: CONFIRMED
- Most common: TRIQS/cthyb
- Widely used in DMFT: CONFIRMED
- Users should specify which implementation
