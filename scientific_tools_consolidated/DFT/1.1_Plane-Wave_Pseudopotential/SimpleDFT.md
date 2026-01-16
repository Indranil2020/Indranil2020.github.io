# SimpleDFT

## Official Resources
- Repository: Included in eminus repo / PyPI
- License: Apache License 2.0

## Overview
SimpleDFT is a minimalist implementation of plane-wave Density Functional Theory in Python. It serves as the pedagogical prototype for the larger `eminus` project. It is stripped down to the absolute essentials to demonstrate the logic of a DFT code in the fewest possible lines of readable Python.

**Scientific domain**: Education
**Target user community**: Absolute beginners to DFT implementation

## Key Features
- **Minimalism**: Focus on code brevity and clarity.
- **Basis**: Plane-waves.
- **Potentials**: GTH Pseudopotentials.
- **Potentials**: GTH Pseudopotentials.
- **Dependency**: Only NumPy.

## Theoretical Methods
- **Basis**: Plane-wave expansion.
- **Algorithm**: Direct minimization or simple SCF mixing.
- **Grid**: Uniform real-space grid for dual-space operations.

## Computational Cost
- **High**: Not optimized for speed; purely for demonstration.
- **Scaling**: Limited to tiny systems (1-10 atoms) due to $O(N^3)$ diagonalization in pure Python.

## Best Practices
- **Do not use for production**: This code is strictly for understanding the internal mechanics of a DFT code.
- **Compare**: Run alongside `eminus` or standard literature values to verify understanding.

## Usage
- primarily used as a learning resource to read alongside DFT textbooks.
- Can run simple atoms (He/H) and small molecules.

## Relation to eminus
- `eminus` is the production-ready(er) version with features and optimizations. `SimpleDFT` is the "skeleton" version.

## Verification
## Community and Support
- **Source**: Part of the eminus ecosystem.
- **Support**: via eminus GitHub issues.

## Verification & Sources
**Primary sources**:
1. eminus Repository
2. PyPI: https://pypi.org/project/SimpleDFT/

**Confidence**: VERIFIED - Exists as a teaching aid.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: Education
- Key Feature: 100% Python
