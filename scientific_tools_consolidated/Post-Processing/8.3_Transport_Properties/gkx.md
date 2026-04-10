# gkx

## Official Resources
- Source Repository: https://github.com/sirmarcel/gkx
- License: See repository

## Overview
gkx is a JAX-based workflow for computing thermal conductivity using Green–Kubo relations, designed to integrate with modern machine-learned interatomic potentials and to support efficient anharmonic thermal transport calculations from molecular dynamics trajectories.

**Scientific domain**: Thermal transport, Green–Kubo, anharmonicity, ML potentials  
**Target user community**: Researchers computing thermal conductivity from MD using modern JAX/ML potential stacks

## Theoretical Methods
- Green–Kubo (fluctuation–dissipation) relations for thermal conductivity
- Equilibrium molecular dynamics (EMD) workflow components

## Capabilities (CRITICAL)
- Thermal conductivity from EMD trajectories
- Workflow tooling for Green–Kubo correlation functions (as implemented)
- JAX-based pipelines enabling integration with differentiable/ML workflows

## Inputs & Outputs
- **Input formats**: Trajectory/heat-flux or on-the-fly heat flux workflow inputs (see repository)
- **Output data types**: Correlation functions and thermal conductivity estimates; tabulated outputs

## Interfaces & Ecosystem
- JAX ecosystem; intended to be combined with ML potential tooling (see repository)

## Limitations & Known Constraints
- Requires careful statistical convergence and validation against reference calculations.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/sirmarcel/gkx

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
