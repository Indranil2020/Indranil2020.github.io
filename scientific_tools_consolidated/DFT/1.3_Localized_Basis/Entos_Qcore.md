# Entos Qcore

## Official Resources
- Homepage: https://entos.ai/
- Documentation: https://docs.entos.ai/ (if available) / Academic papers
- Source Repository: Proprietary / Free for Academic
- License: Free Academic License / Commercial

## Overview
Entos Qcore is a modern quantum chemistry software package developed by Entos, Inc. It emphasizes the integration of physics-based methods with machine learning (ML), specifically featuring Molecular Orbital Based Machine Learning (MOB-ML). It is designed for efficiency and ease of use, providing standard DFT and wavefunction methods alongside innovative ML-accelerated approaches.

**Scientific domain**: Quantum chemistry, Machine Learning, DFT  
**Target user community**: Academic researchers, pharmaceutical industry, materials design

## Theoretical Methods
- Density Functional Theory (DFT)
- Hartree-Fock (HF)
- Molecular Orbital Based Machine Learning (MOB-ML)
- Semi-empirical methods (xTB interface)
- Geometry optimization
- Ab initio molecular dynamics
- Composite methods

## Capabilities (CRITICAL)
- Fast DFT calculations
- MOB-ML for high-accuracy energies at low cost
- Reaction path finding
- Conformational search
- Geometry optimization
- Vibrational frequency analysis
- Modern C++ architecture
- Python API (EntosQ)

## Key Strengths

### MOB-ML:
- Machine learning on molecular orbitals
- Transferable accuracy
- Correct physics scaling
- Reduced training data needs
- Accelerates high-level methods

### Modern Design:
- Efficient C++ core
- User-friendly input
- Structured output (JSON)
- Clean Python integration
- robust optimizers

### Performance:
- Optimized integral engines
- Fast SCF convergence
- Efficient threading
- Scalable algorithms

## Inputs & Outputs
- **Input formats**:
  - Qcore input format (structured)
  - XYZ / PDB coordinates
  - Python API calls
- **Output data types**:
  - JSON structured output
  - Energies, Gradients
  - Properties
  - Log files

## Interfaces & Ecosystem
- **Python**: EntosQ python library
- **OrbNet**: Integration with Entos AI tools
- **Visualization**: Standard tools via output formats

## Advanced Features

### Machine Learning Integration:
- Training MOB-ML models
- Using pre-trained models
- Active learning workflows
- Uncertainty estimation

### Geometry Optimization:
- Robust coordinate systems
- Transition state search
- Constrained optimization
- Reaction path following

## Performance Characteristics
- **Speed**: Competitive with modern C++ codes
- **Accuracy**: DFT and ML-corrected accuracies
- **System size**: Medium to large molecules (ML accelerated)
- **Memory**: Efficient handling
- **Parallelization**: SMP/OpenMP

## Computational Cost
- **DFT**: Standard scaling O(N^3-N^4)
- **MOB-ML**: Significantly faster than target method (e.g., CCSD(T))
- **Training**: One-time cost for ML models
- **Inference**: Very fast

## Limitations & Known Constraints
- **License**: Proprietary (Free Academic)
- **Source**: Not open source
- **Methods**: Focused on DFT/ML, less legacy method coverage
- **Documentation**: Access may depend on license

## Comparison with Other Codes
- **vs Gaussian/ORCA**: Entos emphasizes ML integration (MOB-ML)
- **vs Psi4**: Entos has proprietary ML features
- **vs TeraChem**: Both commercial, Entos focuses on ML/CPU, TeraChem on GPU
- **Unique strength**: MOB-ML for high-accuracy acceleration

## Application Areas

### Drug Discovery:
- **Binding Affinity**: High-throughput calculation of accurate protein-ligand binding energies
- **Conformational Analysis**: Rapid screening of low-energy conformers using ML potentials
- **pKa Prediction**: Fast and accurate free energy calculations
- **Virtual Screening**: Filtering large compound libraries with physics-based accuracy

### Catalysis & Reaction Engineering:
- **Transition States**: Efficient location of transition structures for reaction mechanisms
- **Barrier Heights**: Accurate activation energies using ML-corrected DFT
- **Catalyst Screening**: High-throughput evaluation of catalyst performance
- **Mechanism Exploration**: Interactive path searching on ML surfaces

## Best Practices
### Licensing & Setup:
- **Academic Access**: Register for the free academic license to access full features
- **Installation**: Use provided binary packages or containers for immediate deployment
- **Environment**: Set up creating isolated conda/python environments for EntosQ

### Machine Learning Workflows:
- **Training Data**: Ensure training sets cover the relevant chemical space of your problem
- **Validation**: Always validate ML predictions against full DFT for a subset of systems
- **Active Learning**: Use uncertainty queries to iteratively improve models

### Calculation Strategy:
- **Geometry Opt**: Use robust internal coordinates for flexible molecules
- **SCF Convergence**: Utilize appropriate DIIS or second-order convergence accelerators
- **Python API**: Leverage EntosQ for complex, multi-step workflows not possible in single inputs

## Community and Support
- Entos, Inc. support
- Academic user community
- Caltech/Miller group origins
- Publications and webinars

## Verification & Sources
**Primary sources**:
1. Homepage: https://entos.ai/
2. Miller group publications (Caltech) on MOB-ML
3. Manby et al. publications

**Confidence**: VERIFIED
- Status: Active commercial/academic software
- Technology: MOB-ML published
- Developer: Entos, Inc.
