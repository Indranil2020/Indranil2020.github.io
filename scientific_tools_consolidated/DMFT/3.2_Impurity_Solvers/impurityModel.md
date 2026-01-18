# impurityModel

## Official Resources
- Homepage: https://github.com/JohanSchott/impurityModel
- Documentation: https://github.com/JohanSchott/impurityModel
- Source Repository: https://github.com/JohanSchott/impurityModel
- License: MIT License

## Overview
impurityModel is a Julia/Python-based exact diagonalization (ED) solver for the Anderson impurity model. It is specialized for simulating core-level spectroscopies (like XPS, XAS, RIXS) where finite-size effects of ED are less critical or where multiplet effects are dominant. It allows for the accurate simulation of local many-electron physics in core-level spectroscopy, particularly for correlated materials like transition metal oxides.

**Scientific domain**: Quantum impurity models, Spectroscopy (XPS, XAS, RIXS)
**Target user community**: Spectroscopists, DMFT researchers needing ED solvers

## Theoretical Methods
- Anderson Impurity Model (AIM)
- Exact Diagonalization (Lanczos algorithm)
- Core-level spectroscopy simulation
- Charge-transfer multiplet theory
- Spin-orbit coupling

## Capabilities (CRITICAL)
- **Spectral Functions**: Calculates static and dynamic correlation functions (Green's functions).
- **Spectroscopy**: Specialized for:
    - **XPS**: X-ray Photoemission Spectroscopy.
    - **XAS**: X-ray Absorption Spectroscopy (L-edges, K-edges).
    - **RIXS**: Resonant Inelastic X-ray Scattering (photon-in, photon-out).
- **Hamiltonian**: Flexible definition of impurity Hamiltonians including full Coulomb interaction matrix and bath hybridization.
- **Finite T**: Finite temperature calculations (via Lanczos).

## Key Features

### Lightweight ED:
- An efficient, lightweight alternative to heavy CT-QMC for problems suitable for ED (small bath sizes).
- Ideal for atomic-like limits and core-level physics where multiplets dominate.

### Julia/Python Interface:
- Easy to define models and parameters.
- Scriptable analysis workflows.

## Inputs & Outputs
- **Input formats**:
  - Model parameters (hopping `t`, interaction `U`, `J`, `slater`, hybridization `V`) defined in script.
- **Output data types**:
  - Spectra (Intensity vs Energy).
  - Ground state wavefunctions.

## Interfaces & Ecosystem
- **Standalone**: Can be used standalone for spectroscopy simulation.
- **Integration**: Potentially interfaced as a solver for DMFT loops when bath sizes are small.

## Performance Characteristics
- **Cost**: Exponential scaling with number of bath sites (typical ED limitation).
- **Speed**: Fast for small clusters standard in multiplet calculations.

## Comparison with Other Codes
| Feature | impurityModel (Johan Schott) | Quanty |
| :--- | :--- | :--- |
| **Type** | Python/Julia Library | Scripting Language (Lua based) |
| **Core Algorithm** | Exact Diagonalization (Lanczos) | Exact Diagonalization (Lanczos/Full) |
| **Specialization** | Core-level Spectroscopy (XAS, RIXS) for Impurity Models | General Quantum Many-Body Problems & Spectroscopy |
| **Flexibility** | Focused on standard impurity Hamiltonians | Highly flexible custom Hamiltonian definitions |
| **User Interface** | Python/Julia classes | Lua scripts |
| **Performance** | Optimized for small clusters | High performance, supports large basis sets |

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/JohanSchott/impurityModel
2. Literature: Standard application of Anderson Impurity Model to X-ray spectroscopy.

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub)
- Focus: ED solver for spectroscopy
