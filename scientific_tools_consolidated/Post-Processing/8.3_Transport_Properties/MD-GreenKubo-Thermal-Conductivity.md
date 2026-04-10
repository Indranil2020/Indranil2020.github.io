# MD-GreenKubo-Thermal-Conductivity

## Official Resources
- Source Repository: https://github.com/erny123/MD-GreenKubo-Thermal-Conductivity
- License: See repository

## Overview
MD-GreenKubo-Thermal-Conductivity is a set of scripts and examples for computing thermal conductivity from molecular dynamics simulations using the Green–Kubo formalism. It targets the practical workflow of taking heat-flux outputs (e.g., from LAMMPS) and performing correlation/integration analysis.

**Scientific domain**: Thermal transport, molecular dynamics post-processing  
**Target user community**: MD users computing thermal conductivity via Green–Kubo relations

## Theoretical Methods
- Green–Kubo relations for thermal conductivity
- Equilibrium molecular dynamics (EMD) post-processing

## Capabilities (CRITICAL)
- Thermal conductivity evaluation from heat-flux autocorrelation
- Example workflows and scripts for data analysis and fitting (as provided)

## Inputs & Outputs
- **Input formats**: Heat-flux time series from MD (e.g., LAMMPS output) and associated metadata
- **Output data types**: Heat-flux autocorrelation functions; integrated thermal conductivity estimates

## Interfaces & Ecosystem
- Often used with LAMMPS-produced heat flux outputs

## Limitations & Known Constraints
- Requires careful statistical convergence (trajectory length, block averaging) and validation.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/erny123/MD-GreenKubo-Thermal-Conductivity

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
