# kg4vasp

## Official Resources
- Source Repository: https://github.com/conodipaola/kg4vasp
- License: See repository

## Overview
kg4vasp is a post-processing workflow for computing transport properties using the generalized Kubo-Greenwood approach from first-principles molecular dynamics with VASP. It is intended for calculating electrical conductivity and related quantities from VASP matrix elements and/or compatible inputs.

**Scientific domain**: Electronic transport, Kubo-Greenwood conductivity, first-principles MD post-processing  
**Target user community**: VASP users computing conductivity/transport coefficients from MD trajectories

## Theoretical Methods
- Generalized Kubo-Greenwood formalism

## Capabilities (CRITICAL)
- Electrical conductivity (frequency-dependent and/or DC depending on workflow)
- Derived transport quantities supported by the workflow (see repository)

## Inputs & Outputs
- **Input formats**: VASP outputs and/or intermediate files as required by kg4vasp (see repository README)
- **Output data types**: Transport coefficients in tabulated form; plots/scripts (repository dependent)

## Interfaces & Ecosystem
- VASP focused; includes patching instructions for specific versions (see repository)

## Limitations & Known Constraints
- Depends on VASP version compatibility and available matrix element outputs.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/conodipaola/kg4vasp

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
