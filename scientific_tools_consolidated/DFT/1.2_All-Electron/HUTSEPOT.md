# HUTSEPOT

## Official Resources
- Homepage: https://www.jku.at/institut-fuer-theoretische-physik/forschung/abteilung-fuer-vielteilchensysteme/research/hutsepot
- Access: Contact developers at JKU Linz (Registration Required)
- License: Academic License

## Overview
HUTSEPOT is a versatile, all-electron Density Functional Theory (DFT) code based on the Korringa-Kohn-Rostoker (KKR) Green's function method. Developed at the Institute for Theoretical Physics at Johannes Kepler University (JKU) Linz, it is specifically designed to handle complex electronic structure problems involving disorder, surfaces, and nanostructures where periodic boundary conditions might be broken or require special treatment.

**Scientific domain**: Disordered alloys, surface physics, magnetism
**Target user community**: Solid-state physicists, researchers in magnetism and spintronics

## Theoretical Methods
- **KKR Green's Function**: Solves the Kohn-Sham equations using the multiple scattering theory approach.
- **Coherent Potential Approximation (CPA)**: Rigorous treatment of substitutional disorder in alloys without using supercells.
- **Relativistic Effects**: Fully relativistic implementation for accurate description of spin-orbit coupling and heavy elements.
- **Exchange-Correlation**: LDA, GGA, and self-interaction correction (SIC) schemes.

## Capabilities
- **Disordered Systems**: Efficient calculation of electronic structure in random alloys.
- **Low-Dimensional Systems**: Surfaces, interfaces, and clusters embedded in bulk.
- **Magnetism**: Non-collinear magnetism, magnetic anisotropy, and spin excitations.
- **Spectroscopy**: Calculation of Bloch spectral functions to compare with ARPES.

## Key Strengths
### Handling Disorder
- The KKR-CPA implementation is a gold standard for calculating properties of substitutionally disordered alloys, offering advantages over supercell approaches in terms of efficiency and averaging.

### Green's Function Embedding
- Allows for the precise embedding of impurities or clusters into semi-infinite hosts, making it ideal for studying defects and surface adsorbates.

## Inputs & Outputs
- **Input**: Structured text files defining geometry, potentials, and KKR parameters.
- **Output**: DOS, spectral functions, magnetic moments, total energies.

## Interfaces & Ecosystem
- **Development**: Maintained by the Vielteilchensysteme group at JKU Linz.
- **Language**: Fortran.
- **Parallelization**: MPI/OpenMP hybrid parallelization.

## Computational Cost
- **Scaling**: $O(N^3)$ generally, but prefactors can be large due to multiple scattering.
- **Efficiency**: Very efficient for calculating properties of impurities and disordered systems where supercells would be prohibitively large.

## Verification & Sources
**Primary sources**:
1. JKU Linz Website: https://www.jku.at/institut-fuer-theoretische-physik/forschung/abteilung-fuer-vielteilchensysteme/research/hutsepot
2. "HUTSEPOT: A multi-purpose multiple scattering code" (Author: Arthur Ernst et al.)

## Community and Support
- **Support**: Direct contact with developers (arthur.ernst@jku.at).
- **Availability**: Code is available upon request/registration via GitLab.

**Confidence**: VERIFIED
**Status**: Active, Academic
