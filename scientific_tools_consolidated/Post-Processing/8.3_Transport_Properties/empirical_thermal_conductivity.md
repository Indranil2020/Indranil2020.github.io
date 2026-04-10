# empirical_thermal_conductivity

## Official Resources
- Source Repository: https://github.com/houzf/empirical_thermal_conductivity
- License: See repository

## Overview
empirical_thermal_conductivity is a tool to estimate lattice thermal conductivity using empirical models, including Clarke’s, Cahill–Pohl’s, and Slack’s models. It is designed for quick estimates and comparisons when full phonon-BTE calculations are not required.

**Scientific domain**: Thermal conductivity estimation, empirical models  
**Target user community**: Materials scientists needing fast baseline thermal conductivity estimates

## Theoretical Methods
- Clarke model
- Cahill–Pohl model
- Slack model

## Capabilities (CRITICAL)
- Compute thermal conductivity estimates from empirical relations
- Print supporting quantities such as Debye temperature and related parameters (as implemented)

## Inputs & Outputs
- **Input formats**: Material parameters as defined by the tool
- **Output data types**: Estimated thermal conductivity values and reported intermediate quantities

## Interfaces & Ecosystem
- Standalone estimator; can be integrated into screening workflows

## Limitations & Known Constraints
- Empirical models provide approximate estimates; may not capture strong anharmonicity/complex scattering.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/houzf/empirical_thermal_conductivity

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
