# TEprop2D

## Official Resources
- Source Repository: https://github.com/artnugraha/TEprop2D
- License: See repository

## Overview
TEprop2D is a lightweight Fortran program to calculate thermoelectric transport properties of 2D materials using outputs produced by Quantum ESPRESSO and EPW. It is designed to compute standard thermoelectric quantities as a function of Fermi level and temperature for 2D systems.

**Scientific domain**: Thermoelectrics, 2D materials, electronic transport  
**Target user community**: Quantum ESPRESSO / EPW users studying thermoelectric properties of 2D materials

## Theoretical Methods
- Semi-classical transport coefficient evaluation from EPW-derived quantities (workflow-defined)

## Capabilities (CRITICAL)
- Seebeck coefficient
- Electrical conductivity (workflow-dependent scaling)
- Electronic thermal conductivity (workflow-dependent scaling)
- Thermoelectric power factor and ZT-related quantities (workflow-dependent)
- Directional transport outputs for 2D materials (as implemented)

## Inputs & Outputs
- **Input formats**: Quantum ESPRESSO and EPW outputs required by TEprop2D (see repository README)
- **Output data types**: Transport coefficients vs Fermi level/temperature; tabulated output files

## Interfaces & Ecosystem
- Designed to use Quantum ESPRESSO + EPW as upstream providers of electronic/EP quantities

## Limitations & Known Constraints
- Intended for 2D materials workflows; accuracy depends on EPW inputs and convergence.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/artnugraha/TEprop2D

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
