# topological-insulator-spin-hall

## Official Resources
- Source Repository: https://github.com/smfarzaneh/topological-insulator-spin-hall
- License: See repository

## Overview
This repository provides an ab initio workflow to compute spin Hall conductivity in topological insulators using Quantum ESPRESSO and Wannier90-derived interpolations (as described by the repository). It is a practical implementation-oriented code base for spin Hall transport calculations.

**Scientific domain**: Spin Hall conductivity, topological materials, electronic transport  
**Target user community**: Researchers performing spin Hall conductivity calculations from first-principles

## Theoretical Methods
- Wannier interpolation based evaluation of Berry-curvature-related response functions (workflow dependent)
- Linear-response transport for spin Hall conductivity (as implemented)

## Capabilities (CRITICAL)
- Spin Hall conductivity calculation workflow
- Integration of Quantum ESPRESSO + Wannier90 outputs

## Inputs & Outputs
- **Input formats**: Quantum ESPRESSO and Wannier90 workflow files as required by the repository
- **Output data types**: Spin Hall conductivity values and intermediate data products

## Interfaces & Ecosystem
- Quantum ESPRESSO
- Wannier90

## Limitations & Known Constraints
- Requires careful convergence of k-point meshes and Wannierization.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/smfarzaneh/topological-insulator-spin-hall

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
