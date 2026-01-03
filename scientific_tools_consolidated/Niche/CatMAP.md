# CatMAP

## Official Resources
- Homepage: https://catmap.readthedocs.io/
- Documentation: https://catmap.readthedocs.io/
- Source Repository: https://github.com/SUNCAT-Center/CatMAP
- License: GPL v3

## Overview
CatMAP is a software package for thermodynamic and kinetic modeling of catalytic reactions. It allows users to create microkinetic models based on DFT-calculated energies. CatMAP automates the solution of the mean-field rate equations to predict turnover frequencies (TOF), coverages, and reaction rates as a function of temperature and pressure.

**Scientific domain**: Heterogeneous catalysis, microkinetic modeling  
**Target user community**: Catalysis researchers

## Capabilities (CRITICAL)
- **Microkinetics**: Solves steady-state rate equations.
- **Scaling Relations**: Uses linear scaling relations/BEP relations to estimate energies across materials.
- **Parser**: Reads energy inputs from text files (table format).
- **Analysis**: Rate control analysis, degree of rate control, coverage maps.
- **Phase Diagrams**: Surface coverage phase diagrams.

**Sources**: CatMAP documentation, Comp. Phys. Comm. 204, 206 (2016)

## Inputs & Outputs
- **Input formats**: Python setup script, energy table (txt)
- **Output data types**: Rates, coverages, plots

## Interfaces & Ecosystem
- **ASE**: Used for some thermodynamic utilities.
- **Matplotlib**: For plotting results.

## Workflow and Usage
1. Define reaction mechanism (elementary steps).
2. Provide formation energies of intermediates (from DFT).
3. `model = ReactionModel(setup_file='setup.mkm')`
4. `model.run()`
5. Analyze TOF vs T/P.

## Performance Characteristics
- Fast solution of algebraic/differential equations.
- Bottleneck is usually gathering DFT data.

## Application Areas
- Catalyst screening
- Understanding reaction mechanisms
- Volcano plots

## Community and Support
- Developed by SUNCAT (Stanford/SLAC)
- Active user base in catalysis

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/SUNCAT-Center/CatMAP
2. Publication: A. J. Medford et al., Comp. Phys. Comm. 204, 206 (2016)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Microkinetic modeling, catalysis
