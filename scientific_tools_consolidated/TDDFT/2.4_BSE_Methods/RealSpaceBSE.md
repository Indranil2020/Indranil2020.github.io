# Real-Space-BSE

## Official Resources
- **Source Code**: https://github.com/AlexBuccheri/Bethe-Salpeter
- **Developer**: Alex Buccheri
- **License**: See Repository

## Overview
**Real-Space-BSE** is a specialized Fortran 2003 library implementing the **Bethe-Salpeter Equation (BSE)** in a **real-space** formalism. Unlike many plane-wave or k-space codes, this tool is specifically designed for large **molecular systems**, utilizing a localized basis set to efficiently compute optical absorption spectra and optical band gaps for systems with thousands of atoms.

## Theoretical Methods
- **Real-Space Formulation**: Solves the BSE directly in real-space coordinates.
- **Minimal Tight-Binding Basis**: Uses a simplified basis for efficiency.
- **Two-Center Approximation**: Approximates integrals to speed up computation.
- **Full BSE**: Solves the full equation, not restricted to the Tamm-Dancoff Approximation (TDA).
- **Monopole Approximation**: Used for Coulomb terms beyond nearest neighbors.

## Capabilities
- **Optical Spectra**: Calculation of optical absorption cross-sections.
- **Optical Band Gaps**: Direct determination of optical gaps including excitonic binding.
- **Large Systems**: Tested on molecular systems up to **6000 atoms** (as of 2017).
- **Molecule-Only**: Specifically for finite molecular systems, not periodic solids.

## Implementation & Tech Stack
- **Language**: Fortran 2003 (F2003 library).
- **Status**: Research code / Library.
- **Parallelization**: Real-space methods are naturally parallelizable, though specific MPI/OpenMP details depend on the driver.

## Application Areas
- **Large Molecules**: Organic semiconductors, polymers, large biomolecules.
- **Nanostructures**: Finite quantum dots or clusters where periodic boundaries are artificial.
- **Excitonics**: Study of spatial exciton extent in real space.
