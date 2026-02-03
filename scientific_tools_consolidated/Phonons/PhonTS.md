# PhonTS (Phonon Transport Simulator)

## Official Resources
- **Repository**: [CPC Program Library](https://data.mendeley.com/datasets/vn5wph88hy/1) (Mendeley Data / CPC)
- **Paper**: [Computer Physics Communications 192, 196-204 (2015)](https://doi.org/10.1016/j.cpc.2015.01.017)
- **License**: Standard CPC License (typically free for academic use)

## Overview
PhonTS (Phonon Transport Simulator) is a software package for calculating lattice thermal conductivity from first principles. Developed by Aleksandr V. Chernatynskiy and colleagues, it solves the phonon Boltzmann transport equation (BTE) using input from molecular dynamics or lattice dynamics calculations. It is distributed through the Computer Physics Communications (CPC) Program Library.

**Scientific domain**: Lattice thermal conductivity, Phonon transport
**Target user community**: Materials scientists, Thermal transport researchers

## Theoretical Methods
- Phonon Boltzmann Transport Equation (BTE)
- Iterative solution for BTE
- Relaxation Time Approximation (RTA)
- Integration with molecular dynamics (for lifetimes/scattering)
- Isotope scattering effects

## Capabilities (CRITICAL)
- Calculation of lattice thermal conductivity
- Phonon lifetime analysis
- Spectral thermal conductivity
- Accumulation of thermal conductivity with mean free path
- Anisotropic thermal conductivity tensors

## Key Strengths
- **BTE solver**: Full iterative and RTA solutions
- **Spectral analysis**: Mode-resolved thermal conductivity
- **Accumulation functions**: Mean free path analysis
- **Published methodology**: Well-documented in CPC publication

## Inputs & Outputs
- **Inputs**: Phonon frequencies, group velocities, lifetimes (scattering rates)
- **Outputs**: Thermal conductivity tensor, spectral contributions, accumulation functions

## Interfaces & Ecosystem
- **Force constants**: From phonopy, phono3py, or MD simulations
- **MD codes**: Interfaces for lifetime extraction
- **Standalone**: Self-contained BTE solver

## Performance Characteristics
- **Computational cost**: Moderate; depends on q-point grid
- **Scalability**: Handles standard phonon transport calculations
- **Purpose**: Research-grade thermal transport

## Computational Cost
- Force constant/lifetime generation: External (DFT or MD)
- PhonTS BTE solution: Minutes to hours
- Iterative BTE more expensive than RTA

## Limitations & Known Constraints
- **Requires external input**: Phonon properties from other codes
- **Documentation**: Primarily via publication
- **Community**: Academic user base
- **Development**: CPC library distribution

## Comparison with Other Codes
- **vs ShengBTE/phono3py**: PhonTS alternative BTE solver
- **vs Phoebe**: PhonTS more traditional approach
- **Use case**: Academic research, methodology comparison

## Application Areas
- Lattice thermal conductivity research
- Phonon transport analysis
- Spectral thermal conductivity studies
- Mean free path accumulation analysis
- Academic thermal transport research

## Best Practices
- Converge phonon properties before BTE solution
- Test RTA vs iterative convergence
- Validate against experimental thermal conductivity
- Appropriate q-point grid for convergence

## Community and Support
- CPC Program Library distribution
- Academic support via authors
- Publication-based documentation
- Free for academic use

## Development
- Aleksandr V. Chernatynskiy and colleagues
- CPC library maintenance
- Academic research code

## Research Impact
PhonTS provides a well-documented BTE solver for lattice thermal conductivity calculations, enabling detailed spectral and accumulation analysis of phonon transport.

## Verification & Sources
**Primary sources**:
1.  A. V. Chernatynskiy et al., "Phonon Transport Simulator (PhonTS)", *Comput. Phys. Commun.* **192**, 196 (2015).
2.  Mendeley Data: https://data.mendeley.com/datasets/vn5wph88hy/1

**Confidence**: VERIFIED
**Verification status**: ✅ VERIFIED
- **Status**: Academic research code (CPC Library).
- **Documentation**: Described in the associated publication.
