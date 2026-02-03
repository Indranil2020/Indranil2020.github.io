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

## Inputs & Outputs
- **Inputs**: Phonon frequencies, group velocities, lifetimes (scattering rates)
- **Outputs**: Thermal conductivity tensor, spectral contributions, accumulation functions

## Verification & Sources
**Primary sources**:
1.  A. V. Chernatynskiy et al., "Phonon Transport Simulator (PhonTS)", *Comput. Phys. Commun.* **192**, 196 (2015).
2.  Mendeley Data: https://data.mendeley.com/datasets/vn5wph88hy/1

**Confidence**: VERIFIED
**Verification status**: ✅ VERIFIED
- **Status**: Academic research code (CPC Library).
- **Documentation**: Described in the associated publication.
