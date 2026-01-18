# COMSUITE

## Official Resources
- Homepage: https://github.com/rutgersphysics/COMSUITE
- Documentation: https://github.com/rutgersphysics/COMSUITE/blob/master/README.md
- Source Repository: https://github.com/rutgersphysics/COMSUITE
- License: See repository copyright notice

## Overview
COMSUITE (Combination Suite) is a computational materials physics code for simulating correlated quantum materials using Dynamic Mean Field Theory (DMFT) and its extensions. It provides an integrated framework for DFT+DMFT calculations with sophisticated impurity solvers and multiple methodological approaches including DMFT, cluster DMFT, and Gutzwiller approximations.

**Scientific domain**: Strongly correlated materials, DFT+DMFT, many-body physics  
**Target user community**: Researchers performing advanced DMFT calculations on correlated materials

## Theoretical Methods
- DMFT (single-site and cluster)
- DFT+DMFT with charge self-consistency
- Gutzwiller approximation (rotationally invariant slave boson)
- Continuous-time quantum Monte Carlo solvers
- LDA+DMFT, GGA+DMFT
- GW+DMFT extensions
- Dual fermion approaches

## Capabilities (CRITICAL)
- Comprehensive DFT+DMFT framework
- Multiple DMFT methodologies
- Charge self-consistent calculations
- Cluster DMFT calculations
- Gutzwiller variational approach (ComRISB)
- Advanced impurity solvers integration
- Realistic materials calculations
- Spectral function calculations
- Integration with DFT codes
- Wannier function downfolding

**Sources**: Official COMSUITE repository (https://github.com/rutgersphysics/COMSUITE), confirmed in 6/7 source lists

## Inputs & Outputs
**Input formats**:
- DFT outputs
- Wannier90 projections
- Configuration files
- Interaction parameters

**Output data types**:
- Self-energies
- Green's functions
- Spectral functions
- Observables
- Convergence data

## Interfaces & Ecosystem
- **Part of Comscope suite**: Related to ComDMFT, ComCTQMC, ComRISB
- **DFT integration**: Multiple DFT code support
- **Impurity solvers**: CTQMC and other solvers
- **Rutgers Physics**: Developed by Rutgers research group

## Limitations & Known Constraints
- Complex installation and setup
- Requires deep DMFT expertise
- Documentation limited
- Computational cost high for advanced calculations
- Platform: Linux/Unix
- Active development status varies


## Performance Characteristics
- **Solvers**: Efficient CTQMC and RISB solvers
- **Parallelization**: MPI support for large-scale calculations
- **Gutzwiller**: ComRISB provides fast, low-cost approximate solutions
- **ComLowH**: optimized downfolding to low-energy Hamiltonians

## Comparison with Other Codes
- **vs TRIQS**: COMSUITE is a specialized package for LDA+DMFT calculations with predefined workflows, whereas TRIQS is a general library for building DMFT tools
- **vs EDMFTF**: Both provide DMFT capabilities, but COMSUITE emphasizes the RISB (Gutzwiller) approach alongside DMFT
- **vs DFT+DMFT**: COMSUITE offers a complete suite including Wannier construction and multiple solvers
- **Unique strength**: Integrated Gutzwiller (RISB) and DMFT with charge self-consistency

## Best Practices
- **Workflow**: Use ComWann for basis construction, then ComRISB/ComCTQMC for solving
- **Solver Selection**: Use RISB for quick estimates or wider bandwidths, CTQMC for exact results
- **Self-consistency**: Ensure charge self-consistency for strongly correlated oxides

## Verification & Sources
**Primary sources**:
1. GitHub repository: https://github.com/rutgersphysics/COMSUITE
2. Related Comscope packages and documentation
3. Rutgers University research group

**Secondary sources**:
1. Published papers using COMSUITE
2. Rutgers DMFT research
3. Confirmed in 6/7 source lists (claude, g, gr, k, m, q)

**Confidence**: VERIFIED - Appears in 6 of 7 independent source lists

**Verification status**: ✅ VERIFIED
- Official repository: ACCESSIBLE (GitHub)
- Documentation: LIMITED (README)
- Source code: OPEN
- Development: Rutgers Physics group
- Comprehensive DMFT suite
- Part of Comscope project
