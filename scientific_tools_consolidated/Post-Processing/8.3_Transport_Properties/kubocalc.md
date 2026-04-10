# kubocalc

## Official Resources
- Source Repository: https://github.com/janbbeck/kubocalc
- License: See repository

## Overview
kubocalc is a Kubo-Greenwood based plugin/workflow for Quantum ESPRESSO intended to compute electronic transport coefficients such as electrical conductivity, electronic thermal conductivity, and the Seebeck coefficient from first-principles molecular dynamics (or sampled configurations) and electronic structure data.

**Scientific domain**: Electronic transport, Kubo-Greenwood conductivity, warm dense matter / liquids / disordered systems (typical use cases)  
**Target user community**: Quantum ESPRESSO users computing transport coefficients with Kubo-Greenwood approaches

## Theoretical Methods
- Kubo-Greenwood formalism for linear-response transport
- Evaluation of frequency-dependent and DC transport coefficients (workflow-dependent)

## Capabilities (CRITICAL)
- Electrical conductivity
- Electronic thermal conductivity
- Seebeck coefficient
- Post-processing workflows for sampled atomic configurations (as supported)

## Inputs & Outputs
- **Input formats**: Quantum ESPRESSO calculations and kubocalc inputs as described in the repository
- **Output data types**: Conductivity/transport tensors and derived coefficients; tabulated outputs

## Interfaces & Ecosystem
- Quantum ESPRESSO oriented

## Limitations & Known Constraints
- Requires converged sampling and careful numerical settings for Kubo-Greenwood calculations.

## Verification & Sources
**Primary sources**:
1. Source repository: https://github.com/janbbeck/kubocalc

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: PUBLIC
