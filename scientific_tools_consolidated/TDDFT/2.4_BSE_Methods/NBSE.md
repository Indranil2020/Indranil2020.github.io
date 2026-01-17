# NBSE (NIST Core-Level BSE Solver)

## Official Resources
- **Homepage**: [OCEAN Website](http://feff.phys.washington.edu/OCEAN/) (Distributed as part of OCEAN)
- **Documentation**: Referenced in OCEAN documentation
- **License**: Ask OCEAN developers

## Overview
**NBSE** (NIST Bethe-Salpeter Equation solver) is the core solver engine used within the **OCEAN** (Obtaining Core Excitations using ABINIT and NBSE) package. Originally developed by Eric Shirley at NIST, it is designed for the accurate calculation of **core-level spectroscopic properties** (XAS, XES, NRIXS) using the Bethe-Salpeter Equation. It is typically not used as a standalone user-facing tool but is the computational heart of the OCEAN workflow.

## Scientific Domain
- **Core-Level Spectroscopy**: XAS, XES
- **Bethe-Salpeter Equation**: Core-hole implementation
- **X-ray Physics**: Synchrotron theoretical support

## Relationship with OCEAN
- **OCEAN**: The user-facing package that interfaces DFT (ABINIT/QE) with NBSE.
- **NBSE**: The numeric solver that handles the BSE matrix construction and diagonalization for core excitations.

## Recommendation
> **Use OCEAN**: For researchers wanting to use NBSE, the correct approach is to install and use the **OCEAN** package (#100), which wraps NBSE with necessary pre/post-processing and DFT interfaces.

## Verification
- **Status**: ✅ VERIFIED (as component)
- **Identity**: Confirmed as "Obtaining Core Excitations using ABINIT and **NBSE**".
- **Primary Use**: Backend for OCEAN.
