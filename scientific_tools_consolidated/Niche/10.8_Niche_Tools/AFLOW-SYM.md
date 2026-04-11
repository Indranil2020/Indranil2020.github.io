# AFLOW-SYM

## Official Resources
- Homepage: http://aflow.org/
- Documentation: http://aflow.org/aflow-sym/
- Source Repository: Part of AFLOW codebase
- License: GPL v3

## Overview
AFLOW-SYM is a robust symmetry analysis tool integrated into the AFLOW framework. It determines the symmetry of crystals (space group, Pearson symbol, Wyckoff positions) and performs symmetry-based operations like finding the primitive cell or standardization. It is designed to handle noisy experimental data and high-throughput calculations robustly.

**Scientific domain**: Crystallography, symmetry analysis  
**Target user community**: Materials scientists, crystallographers

## Capabilities (CRITICAL)
- **Symmetry Determination**: Space group, point group, crystal system.
- **Tolerance**: Adaptive tolerance handling for distorted structures.
- **Standardization**: Conversion to standard conventional/primitive cells.
- **Wyckoff Positions**: Identification of site symmetries.
- **Brillouin Zone**: Calculation of high-symmetry k-paths for band structures.

**Sources**: AFLOW website, Chem. Mater. 22, 585 (2010)

## Inputs & Outputs
- **Input formats**: POSCAR, AFLOW input
- **Output data types**: JSON, text report

## Interfaces & Ecosystem
- **AFLOW**: Part of the suite.
- **VASP**: Used for standardizing inputs.

## Workflow and Usage
1. `aflow --sym --poscar=POSCAR`
2. Returns symmetry report.

## Performance Characteristics
- High speed
- Robust against numerical noise

## Application Areas
- High-throughput database generation
- Structure classification
- Band structure path generation

## Community and Support
- Developed by Curtarolo Group (Duke)

## Verification & Sources
**Primary sources**:
1. Homepage: http://aflow.org/
2. Publication: R. H. Taylor et al., Comp. Mater. Sci. 151, 104 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL)
- Development: ACTIVE
- Applications: Symmetry analysis
