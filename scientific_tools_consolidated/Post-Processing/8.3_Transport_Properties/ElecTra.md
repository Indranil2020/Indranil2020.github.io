# ElecTra (ELECTRA)

## Official Resources
- Source Repository: https://github.com/PatrizioGraziosi/ELECTRA
- Preprint: https://arxiv.org/abs/2208.00745
- License: See repository

## Overview
ElecTra (often referenced as ELECTRA) is an open-source code for computing electronic and thermoelectric transport coefficients from a full electronic band structure by solving the linearized Boltzmann transport equation (BTE) under momentum-relaxation-time-type approximations. It targets semiconductor transport where scattering rates can depend on carrier energy, momentum, and band index.

**Scientific domain**: Electronic transport, thermoelectrics, mobility modeling  
**Target user community**: Researchers computing conductivity, mobility, Seebeck coefficient and related coefficients from first-principles band structures

## Theoretical Methods
- Linearized Boltzmann transport equation (electrons/holes)
- Momentum relaxation time approximation (energy- and k-dependent scattering)
- Transport integrals for conductivity and thermoelectric tensors

## Capabilities (CRITICAL)
- Electrical conductivity
- Seebeck coefficient
- Electronic thermal conductivity
- Carrier mobility
- Unipolar and bipolar transport regimes (depending on model setup)
- Scattering models (as supported by the code) for phonons/impurities via relaxation-time type inputs

## Inputs & Outputs
- **Input formats**: Band structure / k-grid dependent quantities as required by the code workflow (see repository documentation)
- **Output data types**: Temperature and chemical-potential dependent transport coefficients; tabulated outputs for analysis/plotting

## Interfaces & Ecosystem
- Designed for workflows using DFT-derived electronic structures (user-prepared inputs)

## Performance Characteristics
- Cost dominated by k-grid sampling and evaluation of transport integrals

## Limitations & Known Constraints
- Accuracy depends on the chosen scattering model and the quality/density of the underlying electronic structure.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/PatrizioGraziosi/ELECTRA
2. arXiv:2208.00745

**Confidence**: VERIFIED (public repository + citable preprint)

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
- Documentation: In repository
