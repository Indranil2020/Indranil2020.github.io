# TransOpt

## Official Resources
- Source Repository: https://github.com/yangjio4849/TransOpt
- License: See repository

## Overview
TransOpt is a transport post-processing package intended for VASP users to compute electrical transport coefficients, including Seebeck coefficients, electrical conductivities, and electronic thermal conductivities. The repository documentation describes two approaches: a momentum-matrix based method and a derivative method similar in spirit to the approach used by BoltzTraP-type workflows.

**Scientific domain**: Electronic transport, thermoelectrics (electronic part)  
**Target user community**: VASP users needing transport coefficients from band-structure data

## Theoretical Methods
- Semi-classical transport coefficient evaluation from electronic structure
- Derivative-based evaluation of band velocities (per code workflow)
- Momentum matrix method (per code workflow)

## Capabilities (CRITICAL)
- Seebeck coefficient
- Electrical conductivity (often reported as σ/τ depending on workflow)
- Electronic thermal conductivity (often reported as κe/τ depending on workflow)
- Transport vs temperature and chemical potential (workflow-dependent)

## Inputs & Outputs
- **Input formats**: VASP outputs required by TransOpt (see repository README)
- **Output data types**: Tabulated transport coefficients suitable for plotting and thermoelectric analysis

## Interfaces & Ecosystem
- VASP-focused post-processing

## Limitations & Known Constraints
- Requires VASP outputs at sufficient k-point density for converged derivatives/integrals.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/yangjio4849/TransOpt

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
