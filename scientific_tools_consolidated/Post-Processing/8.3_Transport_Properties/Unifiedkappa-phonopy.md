# Unifiedkappa-phonopy

## Official Resources
- Source Repository: https://github.com/yimavxia/Unifiedkappa-phonopy
- License: See repository

## Overview
Unifiedkappa-phonopy is a set of scripts and a modified phonopy-related workflow for analyzing lattice thermal conductivity, including separating diagonal and off-diagonal contributions. It is typically used in conjunction with phonon-BTE solvers (e.g., ShengBTE) and phonon workflow data.

**Scientific domain**: Phonon transport analysis, lattice thermal conductivity  
**Target user community**: Researchers analyzing phonon thermal conductivity beyond standard diagonal approximations

## Theoretical Methods
- Post-processing/analysis of phonon thermal transport quantities
- Decomposition of thermal conductivity contributions (as implemented)

## Capabilities (CRITICAL)
- Analysis workflows for thermal conductivity contributions
- Utilities supporting data preparation and interpretation for phonon transport studies

## Inputs & Outputs
- **Input formats**: Data products from phonon workflows (see repository documentation)
- **Output data types**: Thermal conductivity decompositions and tabulated results

## Interfaces & Ecosystem
- Typically used with phonopy/phono3py and phonon-BTE workflows

## Limitations & Known Constraints
- Intended as an analysis layer; requires upstream phonon calculations.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/yimavxia/Unifiedkappa-phonopy

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
