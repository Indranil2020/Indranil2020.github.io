# CT-QMC

## Official Resources
- Homepage: Multiple implementations (w2dynamics primary reference)
- Documentation: See specific implementations
- Source Repository: https://github.com/w2dynamics/w2dynamics (w2dynamics)
- License: Varies by implementation

## Overview
CT-QMC (Continuous-Time Quantum Monte Carlo) is a general term for the family of continuous-time quantum Monte Carlo algorithms used as impurity solvers in DMFT. This includes CT-HYB, CT-INT, CT-AUX, and other variants. The term CT-QMC often refers generically to these methods, with w2dynamics being a major implementation providing multiple CT-QMC algorithms.

**Scientific domain**: Quantum impurity solvers, DMFT, CTQMC algorithms  
**Target user community**: Researchers performing DMFT calculations

## Theoretical Methods
- Continuous-time quantum Monte Carlo (general family)
- Hybridization expansion (CT-HYB)
- Interaction expansion (CT-INT)
- Auxiliary field expansion (CT-AUX)
- Segment representation
- Worm sampling algorithms
- Stochastic methods for impurity problems

## Capabilities (CRITICAL)
**Note**: CT-QMC is a FAMILY OF ALGORITHMS, not single software

**Major implementations**:
- **w2dynamics**: Comprehensive CT-QMC package (CT-HYB, CT-INT)
- **TRIQS/cthyb**: CT-HYB variant
- **iQIST**: Multiple CT-QMC algorithms
- **ALPS**: CT-HYB implementation
- **ComCTQMC**: GPU-accelerated CT-QMC

**General capabilities** (algorithm-dependent):
- Multi-orbital impurity problems
- Temperature-dependent calculations
- Various interaction types
- Green's functions and observables
- Statistical sampling methods

**Sources**: Master list notes: "VERIFIED - w2dynamics solver", CT-QMC family widely used

## Inputs & Outputs
**Implementation-dependent**: Each variant has specific I/O

**Common elements**:
- Input: Hybridization/Weiss functions, interactions, temperature
- Output: Green's functions, self-energies, occupations

## Interfaces & Ecosystem
- **w2dynamics**: Complete CT-QMC solver package
- **TRIQS**: CT-HYB via TRIQS/cthyb
- **DCore**: Interfaces to multiple CT-QMC solvers
- **solid_dmft**: Can use various CT-QMC implementations
- **DMFT frameworks**: Standard impurity solver approach

## Limitations & Known Constraints
- Computational cost scales with inverse temperature
- Sign problem in certain regimes
- Statistical errors from Monte Carlo
- Algorithm-specific limitations
- Memory intensive for large problems
- Requires substantial compute resources

## Verification & Sources
**Primary sources**:
1. w2dynamics: https://github.com/w2dynamics/w2dynamics
2. E. Gull et al., Rev. Mod. Phys. 83, 349 (2011) - CT-QMC review
3. A. N. Rubtsov et al., Phys. Rev. B 72, 035122 (2005) - CT-INT
4. P. Werner et al., Phys. Rev. B 74, 155107 (2006) - CT-HYB
5. Master list: "VERIFIED - w2dynamics solver"

**Secondary sources**:
1. Multiple CT-QMC implementations
2. Standard in DMFT community
3. Extensive literature on CT-QMC methods

**Confidence**: VERIFIED - Family of algorithms with many implementations

**Verification status**: ✅ VERIFIED as ALGORITHM FAMILY
- Status: Family of algorithms, not single software
- Multiple implementations: CONFIRMED
- w2dynamics: Major reference implementation
- Widely used in DMFT: CONFIRMED
- Users should specify which CT-QMC variant
