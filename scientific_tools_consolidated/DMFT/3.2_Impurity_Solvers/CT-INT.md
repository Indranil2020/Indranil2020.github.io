# CT-INT

## Official Resources
- Homepage: Multiple implementations (ALPSCore/CT-INT, w2dynamics, ComCTQMC)
- Documentation: See specific implementations
- Source Repository: https://github.com/ALPSCore/CT-INT (ALPSCore version)
- License: Varies by implementation (ALPSCore: GPL v2)

## Overview
CT-INT (Continuous-Time Interaction Expansion) is a quantum Monte Carlo algorithm for solving impurity problems, implemented in several DMFT software packages. Unlike CT-HYB which expands in the hybridization, CT-INT expands in the interaction term. Multiple implementations exist including ALPSCore/CT-INT, w2dynamics, iQIST, and ComCTQMC.

**Scientific domain**: Quantum impurity solvers, DMFT, CTQMC algorithms  
**Target user community**: Researchers performing DMFT calculations, especially at weak-to-intermediate coupling

## Theoretical Methods
- Continuous-time quantum Monte Carlo (CTQMC)
- Interaction expansion algorithm
- Stochastic sampling
- Anderson impurity model
- Weak-to-intermediate coupling regime
- Determinant Monte Carlo techniques

## Capabilities (CRITICAL)
**Note**: CT-INT is an ALGORITHM implemented in multiple codes:

**Primary implementations**:
- **ALPSCore/CT-INT**: Open-source, ALPSCore-based
- **w2dynamics**: Includes CT-INT solver
- **iQIST**: Multiple CT-INT variants (begonia, lavender)
- **ComCTQMC**: GPU-accelerated, includes CT-INT

**General capabilities**:
- Multi-orbital impurity problems
- Weak-to-intermediate coupling strength
- Onsite Coulomb interactions
- Complex Weiss functions
- Temperature-dependent calculations
- MPI parallelization

**Sources**: Master list notes: "VERIFIED - CT-INT solver (ComCTQMC)", multiple implementations

## Inputs & Outputs
**Depends on implementation**: Each code has its own interface

**Common inputs**: Weiss function, interaction parameters, temperature
**Common outputs**: Green's functions, self-energies, observables

## Interfaces & Ecosystem
- **ALPSCore/CT-INT**: Standalone using ALPSCore libraries
- **w2dynamics**: Integrated CT-INT option
- **iQIST**: Multiple optimized CT-INT solvers
- **ComCTQMC**: GPU-accelerated version
- **DMFT frameworks**: Used in various DFT+DMFT workflows

## Limitations & Known Constraints
- Sign problem more severe than CT-HYB
- Best for weak-to-intermediate coupling
- Strong coupling can have convergence issues
- Statistical noise from Monte Carlo
- Implementation-specific constraints
- Typically slower than CT-HYB for strong coupling

## Comparison with Other Algorithms
| Algorithm | CT-INT | CT-HYB |
| :--- | :--- | :--- |
| **Expansion Parameter** | Interaction ($U$) | Hybridization ($\Delta$) |
| **Best Regime** | Weak Coupling | Strong Coupling |
| **Complexity** | $O(k^3)$ in perturbation order | $O(k^3)$ in perturbation order |
| **Sign Problem** | Can be severe at large $U$ | Generally better at large $U$ |

## Verification & Sources
**Primary sources**:
1. ALPSCore/CT-INT: https://github.com/ALPSCore/CT-INT
2. A. N. Rubtsov et al., Phys. Rev. B 72, 035122 (2005) - Original CT-INT
3. E. Gull et al., Europhys. Lett. 82, 57003 (2008) - Worm sampling
4. Master list: "VERIFIED - CT-INT solver (ComCTQMC)"

**Secondary sources**:
1. w2dynamics documentation
2. iQIST implementations
3. CT-INT algorithm papers

**Confidence**: VERIFIED - Algorithm with multiple implementations

**Verification status**: ✅ VERIFIED as ALGORITHM
- Status: Algorithm/method, not single software
- Multiple implementations: CONFIRMED
- ALPSCore version: Open source
- Complementary to CT-HYB: Different coupling regimes
- Users should specify which implementation
