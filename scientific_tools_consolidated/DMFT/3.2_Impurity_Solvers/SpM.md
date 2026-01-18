# SpM (Sparse Modeling)

## Official Resources
- Homepage: https://github.com/SpM-lab/SpM
- Documentation: https://github.com/SpM-lab/SpM (README/Examples)
- Source Repository: https://github.com/SpM-lab/SpM
- License: Open Source (GPL/MIT)

## Overview
SpM is a software tool for the analytic continuation of imaginary-time Green's functions using **Sparse Modeling**. This approach utilizes the scarcity of information in the Matsubara Green's function to construct a compact representation (Intermediate Representation) and reconstruct the spectral function. It offers an alternative to Maximum Entropy that is often less sensitive to certain types of noise (overfitting) and requires fewer ad-hoc parameters.

**Scientific domain**: Analytic Continuation, Data Science in Physics, Machine Learning
**Target user community**: DMFT/QMC users, Data scientists in physics

## Theoretical Methods
- **Sparse Modeling**: L1-regularization (LASSO) techniques to find sparse solutions.
- **Intermediate Representation (IR)**: Uses the `irbasis` library (or `sparse-ir`) to define a compact, optimal basis for Green's functions tailored to the kernel.
- **Analytic Continuation**: Inversion of the spectral representation using the sparse coefficients.

## Capabilities (CRITICAL)
- **Spectral Reconstruction**: Recovers $A(\omega)$ from $G(\tau)$ or $G(i\omega_n)$.
- **Noise Tolerance**: Designed to handle statistical noise effectively via sparse regularization.
- **Basis Optimization**: Uses singular value expansion (SVE) of the kernel to define the most relevant basis functions (IR basis).
- **Automatic Parameter Selection**: Often determines regularization parameters automatically.

## Key Features

### Intermediate Representation:
- Leverages `irbasis` / `sparse-ir` libraries.
- Compactly represents Green's functions (both imaginary and real domains) with a minimal number of coefficients.

### Algorithmic Robustness:
- Sparse modeling provides a mathematically grounded way to select the most relevant features of the spectrum, avoiding "fitting the noise".

## Inputs & Outputs
- **Inputs**:
  - $G(\tau)$ or $G(i\omega_n)$.
  - Temperature $\beta$.
  - Error estimates.
- **Outputs**:
  - Spectral function $A(\omega)$.
  - Sparse coefficients vector.

## Workflow and Usage
Post-processing of QMC data.
1.  Compute the Intermediate Representation basis for the given $\beta$ and cutoff.
2.  Project input data $G(\tau)$ onto this basis.
3.  Solve the L1-regularized inversion problem to find sparse coefficients.
4.  Reconstruct $A(\omega)$ from coefficients.

## Comparison with Other Codes
| Feature | SpM (Sparse Modeling) | ana_cont | SOM |
| :--- | :--- | :--- | :--- |
| **Methodology** | L1 regularization, singular value decomposition | MaxEnt, Pade (standard) | Stochastic sampling of solutions |
| **Noise Handling** | Excellent (ignores noise components) | Moderate (depends on error interaction) | Good (statistically handled) |
| **Input** | Imaginary-time Green's functions | Matsubara data | Imaginary-time/Matsubara data |
| **Key Strength** | Stability and efficiency | Standard, widely used implementation | Unbiased spectral reconstruction |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/SpM-lab/SpM
2. Literature: J. Otsuki et al., "Sparse modeling approach to analytical continuation...", Phys. Rev. E.
3. Dependencies: `irbasis` or `sparse-ir`, `scikit-learn` (often used for LASSO).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Method: Sparse modeling is a recognized, modern, data-driven technique in the field.
