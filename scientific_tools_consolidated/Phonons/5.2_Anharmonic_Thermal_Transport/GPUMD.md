# GPUMD

## Official Resources
- Homepage: https://gpumd.org/
- Source Repository: https://github.com/brucefan1983/GPUMD
- Documentation: https://gpumd.org/
- License: GPL-3.0

## Overview
GPUMD (Graphics Processing Units Molecular Dynamics) is a highly efficient molecular dynamics package fully implemented on GPUs. It features the neuroevolution potential (NEP) approach for machine learning potentials, enabling accurate and fast simulations of thermal transport properties.

**Scientific domain**: Molecular dynamics, thermal transport, machine learning potentials  
**Target user community**: Researchers studying thermal transport using MD with ML potentials

## Theoretical Methods
- Classical molecular dynamics
- Neuroevolution potentials (NEP)
- Green-Kubo thermal conductivity
- Homogeneous non-equilibrium MD (HNEMD)
- Spectral decomposition of thermal conductivity
- Heat current autocorrelation

## Capabilities (CRITICAL)
- GPU-accelerated MD simulations
- NEP machine learning potentials
- Thermal conductivity calculations
- Spectral thermal conductivity
- HNEMD method
- Green-Kubo method
- Phonon participation ratio
- Modal analysis

## Key Strengths

### GPU Acceleration:
- Fully GPU-native implementation
- Orders of magnitude speedup
- Large system sizes feasible
- Efficient memory usage

### NEP Machine Learning:
- Neuroevolution potential training
- Near-DFT accuracy
- Fast evaluation
- Active learning support

### Thermal Transport:
- Multiple methods (GK, HNEMD)
- Spectral decomposition
- Modal contributions
- Accurate predictions

## Inputs & Outputs
- **Input formats**:
  - xyz structure files
  - NEP potential files
  - run.in control file
  
- **Output data types**:
  - Thermal conductivity
  - Heat current
  - Spectral properties
  - Trajectory files
  - Thermodynamic properties

## Interfaces & Ecosystem
- **NEP training**: Built-in potential training
- **LAMMPS**: Some compatibility
- **ASE**: Python interface available
- **calorine**: Python package for NEP


## Advanced Features
- **NEP training**: Built-in neuroevolution potential development
- **HNEMD method**: Homogeneous non-equilibrium MD for thermal conductivity
- **Spectral decomposition**: Frequency-resolved thermal conductivity
- **Modal analysis**: Phonon mode contributions to transport
- **Multi-GPU support**: Scalable to large systems
- **Active learning**: Efficient training data selection

## Computational Cost
- NEP training: Hours to days (one-time cost)
- MD simulations: Very fast on GPU (millions of atoms)
- Thermal conductivity: Minutes to hours depending on convergence
- Overall: Orders of magnitude faster than DFT-MD

## Performance Characteristics
- **Speed**: Extremely fast on GPU
- **System size**: Millions of atoms feasible
- **Accuracy**: NEP provides DFT-level accuracy
- **Parallelization**: Multi-GPU support

## Limitations & Known Constraints
- Requires NVIDIA GPU
- NEP training requires expertise
- Classical MD limitations apply
- Learning curve for NEP development

## Application Areas
- Thermal conductivity of complex materials
- Disordered and amorphous systems
- Nanostructured materials
- High-temperature properties
- Materials with strong anharmonicity

## Comparison with Other Codes
- **vs LAMMPS**: GPUMD is fully GPU-native; LAMMPS has GPU packages but CPU-centric design
- **vs Phono3py/ShengBTE**: GPUMD uses MD-based methods; others use perturbation theory
- **vs DeepMD-kit**: Both support ML potentials; GPUMD has built-in NEP, DeepMD uses DP
- **Unique strength**: Integrated NEP training + GPU MD + thermal transport in one package

## Best Practices

### NEP Training:
- Use diverse training dataset
- Include thermal expansion data
- Validate phonon dispersions
- Test force/energy predictions

### MD Simulations:
- Equilibrate system thoroughly
- Use sufficient system size
- Run long enough for convergence
- Monitor temperature stability

### Thermal Conductivity:
- Use both GK and HNEMD methods
- Check size convergence
- Validate with experiments
- Analyze spectral contributions

## Community and Support
- Open-source GPL-3.0
- Very active development (Zheyong Fan)
- Comprehensive documentation
- Tutorial examples included
- Growing user community
- Regular updates and new features

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/brucefan1983/GPUMD
2. Z. Fan et al., Comput. Phys. Commun. 218, 10 (2017)
3. Z. Fan et al., Phys. Rev. B 104, 104309 (2021) - NEP

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, GPL-3.0)
- Documentation: Comprehensive
- Active development: Very active
- Academic citations: >500
