# MACE (Multi-Atomic Cluster Expansion)

## Official Resources
- Homepage: https://mace-docs.readthedocs.io/
- Documentation: https://mace-docs.readthedocs.io/
- Source Repository: https://github.com/ACEsuit/mace
- License: MIT License

## Overview
MACE creates fast and accurate machine learning interatomic potentials using higher-order equivariant message passing. It combines the strengths of the Atomic Cluster Expansion (ACE) with message passing neural networks (MPNNs). MACE achieves state-of-the-art accuracy and is designed to be scalable for large simulations.

**Scientific domain**: Machine learning potentials, equivariant neural networks  
**Target user community**: MD users, ML researchers

## Capabilities (CRITICAL)
- **Higher Order Message Passing**: Captures many-body interactions efficiently.
- **Foundation Models**: Pre-trained "universal" potentials (MACE-MP-0) available for the periodic table.
- **Speed**: Optimized for GPU execution.
- **LAMMPS/ASE**: Interfaces for MD.

**Sources**: MACE GitHub, arXiv:2206.07697

## Inputs & Outputs
- **Input formats**: XYZ training data
- **Output data types**: PyTorch models

## Interfaces & Ecosystem
- **PyTorch**: Backend.
- **ASE**: Calculator.
- **LAMMPS**: Production MD.

## Workflow and Usage
1. Download pre-trained model or train your own.
2. `calc = MACECalculator(model_path='MACE.model', device='cuda')`
3. Run MD with ASE or LAMMPS.

## Performance Characteristics
- Very fast inference for an equivariant model.
- High accuracy across diverse chemical spaces.

## Application Areas
- General purpose MD
- Chemistry (reactions)
- Materials discovery

## Community and Support
- Developed by Kovacs, Batatia, et al. (Cambridge/EPFL)
- Rapidly growing community

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ACEsuit/mace
2. Publication: I. Batatia et al., NeurIPS 2022

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: SOTA ML potentials
