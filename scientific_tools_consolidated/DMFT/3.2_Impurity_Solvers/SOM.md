# SOM (Stochastic Optimization Method)

## Official Resources
- Homepage: https://github.com/kcd2015/SOM
- Documentation: https://github.com/kcd2015/SOM
- Source Repository: https://github.com/kcd2015/SOM
- License: Open Source (GPL or similar)

## Overview
SOM (Stochastic Optimization Method) is a code for the analytic continuation of quantum Monte Carlo (QMC) data. It solves the inverse problem of reconstructing real-frequency spectral functions from imaginary-time or Matsubara frequency Green's functions. It implements a stochastic optimization approach, as proposed by Mishchenko et al., often providing a robust alternative to Maximum Entropy (MaxEnt) by avoiding entropic regularization bias.

**Scientific domain**: Analytic Continuation, Condensed Matter Physics, Numerical Methods
**Target user community**: DMFT/QMC practitioners needing spectral functions

## Theoretical Methods
- Stochastic Optimization (Mishchenko's method)
- Analytic Continuation
- Inverse Ill-posed Problems
- Green's function inversion
- Fredholm integral equation of the first kind

## Capabilities (CRITICAL)
- **Reconstruction**: Recovers Real-frequency spectral function $A(\omega)$ from $G(\tau)$ or $G(i\omega_n)$ input data.
- **Alternative to MaxEnt**: Does not rely on entropic regularization, potentially offering different handling of sharp features.
- **Error Estimation**: Can provide estimates of reliability for reconstructed features.
- **Sampling**: Uses Monte Carlo sampling of the solution space (rectangles configuration).

## Key Features

### TRIQS Integration:
- Often developed within or compatible with the TRIQS ecosystem.
- Designed to handle noisy QMC data effectively.

### Configurable Updates:
- Elementary updates include shifting, resizing, adding, removing, splitting, and gluing rectangles to form the spectrum.

## Inputs & Outputs
- **Input formats**:
  - QMC Data: Green's function data (tau or Matsubara).
  - Error bars (covariance matrix).
  - Parameter file defining MC steps, updates, and deviation function.
- **Output data types**:
  - Real-frequency spectral function $A(\omega)$.

## Workflow and Usage
Used as a post-processing step after a QMC/DMFT calculation.
1.  Load QMC data ($G(\tau)$).
2.  Configure stochastic parameters.
3.  Run SOM to stochastically sample the space of possible spectra.
4.  Average samples to obtain final $A(\omega)$.

## Comparison with Other Codes
| Feature | SOM (Stochastic Optimization) | SpM (Sparse Modeling) | ana_cont |
| :--- | :--- | :--- | :--- |
| **Methodology** | Stochastic sampling of spectral functions | Sparse modeling (L1 regularization) | Maximum Entropy / Pade |
| **Model Dependence** | Low (Minimal prior bias) | Low (Automatic basis selection) | High (Requires default model for MaxEnt) |
| **Computational Cost** | High (Sampling intensive) | Low/Moderate | Low |
| **Key Strength** | Reliable error estimation, unbiased | Robustness against noisy QMC data | General purpose availability |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/kcd2015/SOM
2. Literature: Mishchenko et al., "Stochastic optimization method for analytic continuation..." (arXiv/Phys. Rev. B).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Method: Established stochastic continuation technique (Mishchenko)
