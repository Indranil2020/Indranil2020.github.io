# MAISE (Module for Ab Initio Structure Evolution)

## Official Resources
- Homepage: https://github.com/tamercan/MAISE
- Documentation: https://github.com/tamercan/MAISE/wiki
- Source Repository: https://github.com/tamercan/MAISE
- License: Open-source

## Overview
MAISE is a package for evolutionary structure prediction that emphasizes the use of neural network potentials to accelerate the search. By training interatomic potentials on-the-fly during the evolutionary search, MAISE reduces the number of expensive ab-initio calculations required to find the global minimum.

**Scientific domain**: Evolutionary structure prediction, machine learning potentials, materials discovery  
**Target user community**: Computational materials scientists

## Theoretical Methods
- Evolutionary Algorithm
- Neural Network Potentials (NNP)
- On-the-fly machine learning
- Density Functional Theory (verification)
- Structure fingerprinting

## Capabilities (CRITICAL)
- Accelerated structure prediction using NNPs
- Automated training set generation and potential fitting
- Crystal structure evolution
- Interface with VASP for data generation and final verification
- Dimensionality support (clusters, crystals)

**Sources**: MAISE documentation, Phys. Rev. Lett. 110, 245501 (2013) (Reference to method)

## Inputs & Outputs
- **Input formats**: Configuration files, VASP inputs
- **Output data types**: Predicted structures, trained potentials

## Interfaces & Ecosystem
- **VASP**: Primary DFT engine
- **Neural Networks**: Internal implementation or interface

## Performance Characteristics
- Significantly reduces DFT calls compared to standard EA
- ML potential overhead is small compared to DFT

## Application Areas
- Complex unit cells
- Systems where DFT is too expensive for direct global search
- Metastable phase search

## Community and Support
- Open-source
- Developed by Tamer Can / Stony Brook University

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/tamercan/MAISE

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Repository: ACCESSIBLE
- Method: Evolutionary algorithm + Neural Networks
- Applications: Accelerated structure prediction, ML potentials
