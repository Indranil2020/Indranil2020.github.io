# Kondo

## Official Resources
- Homepage: Not applicable - refers to methods/solvers, not standalone software
- Documentation: Various implementations in other codes
- Source Repository: Multiple implementations exist
- License: Varies by implementation

## Overview
"Kondo" in the DMFT context refers to Kondo problem solvers or Kondo physics implementations rather than a specific standalone software package named "Kondo." The Kondo problem is a fundamental quantum impurity problem, and various codes implement Kondo problem solvers (e.g., NRG, exact diagonalization, QMC methods). The master list entry likely refers to Kondo problem solver implementations in research codes.

**Scientific domain**: Kondo physics, quantum impurity problems, DMFT  
**Target user community**: Researchers studying Kondo effect and heavy fermion systems

## Theoretical Methods
- Kondo problem (physics, not software)
- Quantum impurity models
- Kondo exchange interactions
- Heavy fermion physics
- Implemented via: NRG, ED, QMC, NCA methods

## Capabilities (CRITICAL)
**Note**: "Kondo" is a METHOD/PROBLEM, not specific software

- Kondo problem solvers exist in various codes:
  - NRGLjubljana (NRG for Kondo models)
  - Exact diagonalization codes
  - QMC impurity solvers
  - NCA/OCA solvers
  
- Kondo physics can be studied using multiple software packages

**Sources**: Master list notes: "METHOD - Kondo problem solvers are ubiquitous; specific code name likely refers to a research implementation"

## Inputs & Outputs
**Not applicable**: Depends on specific implementation in various codes

## Interfaces & Ecosystem
- **NRGLjubljana**: Dedicated Kondo problem solver
- **DMFT impurity solvers**: Can handle Kondo physics
- **Research implementations**: Various groups have Kondo solvers

## Limitations & Known Constraints
- **CRITICAL**: Not a specific standalone software
- "Kondo" refers to physics/methods, not a code
- Multiple implementations exist across different packages
- No single "Kondo software" to install
- Users should specify which Kondo solver implementation they need

## Verification & Sources
**Primary sources**:
1. Master list: "METHOD - Kondo problem solvers are ubiquitous"
2. Kondo problem is physics concept, not software
3. Multiple codes implement Kondo solvers

**Secondary sources**:
1. NRGLjubljana for NRG-based Kondo solutions
2. Various QMC codes handle Kondo models
3. Kondo physics literature

**Confidence**: UNCERTAIN - Master list explicitly marks as "METHOD"

**Verification status**: ⚠️ METHOD, NOT SOFTWARE
- Status: Physics method/problem, not standalone code
- Multiple implementations: CONFIRMED
- Specific "Kondo" software: DOES NOT EXIST
- Recommendation: Use NRGLjubljana or other impurity solvers for Kondo problems
