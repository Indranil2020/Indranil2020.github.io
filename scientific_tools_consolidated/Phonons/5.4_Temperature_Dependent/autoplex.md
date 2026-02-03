# autoplex

## Official Resources
- Homepage: https://github.com/autoatml/autoplex
- Source Repository: https://github.com/autoatml/autoplex
- Documentation: https://autoatml.github.io/autoplex/
- License: BSD-3-Clause

## Overview
autoplex is an automated workflow for fitting machine learning interatomic potentials (MLIPs) with a focus on phonon properties. It provides end-to-end automation from DFT data generation to potential validation using phonon benchmarks.

**Scientific domain**: Machine learning potentials, phonon validation, automated workflows  
**Target user community**: Researchers developing ML potentials for phonon calculations

## Theoretical Methods
- Machine learning interatomic potentials
- Moment tensor potentials (MTP)
- GAP potentials
- Phonon-based validation
- Active learning
- Automated benchmarking

## Capabilities (CRITICAL)
- Automated MLIP fitting
- Phonon property validation
- DFT data generation workflows
- Multiple MLIP backends
- Benchmark against DFT phonons
- atomate2 integration
- High-throughput workflows

## Key Strengths

### Automation:
- End-to-end workflow
- Minimal user intervention
- Reproducible results
- Best practices built-in

### Phonon Focus:
- Phonon-based validation
- Accurate force constants
- Dispersion benchmarks
- Quality metrics

## Inputs & Outputs
- **Input formats**:
  - Crystal structures
  - DFT settings
  - MLIP parameters
  
- **Output data types**:
  - Trained potentials
  - Phonon benchmarks
  - Validation metrics
  - Comparison plots

## Interfaces & Ecosystem
- **atomate2**: Workflow engine
- **Phonopy**: Phonon calculations
- **VASP/other**: DFT backends
- **MTP/GAP**: MLIP backends


## Advanced Features

### Core Capabilities:
- Detailed feature implementation
- Advanced algorithms and methods
- Specialized functionality
- Integration capabilities

### Performance Optimizations:
- Computational efficiency features
- Scalability enhancements
- Memory management
- Parallel processing support


## Computational Cost
- **Setup**: Preprocessing requirements
- **Main calculation**: Primary computational cost
- **Post-processing**: Analysis overhead
- **Overall**: Total resource requirements

## Limitations & Known Constraints
- Requires DFT infrastructure
- Computational resources needed
- Learning curve for workflows
- Evolving codebase

## Application Areas
- MLIP development
- Phonon-accurate potentials
- High-throughput screening
- Materials discovery
- Thermal property prediction

## Comparison with Other Codes
- **vs manual MLIP fitting**: autoplex automates entire workflow
- **vs GPUMD NEP**: autoplex supports multiple MLIP backends
- **vs traditional workflows**: Built-in phonon validation
- **Unique strength**: End-to-end automation with phonon focus

## Best Practices

### Training Data Generation:
- Include diverse configurations
- Sample relevant temperature range
- Add strained structures
- Include phonon-relevant displacements

### Potential Validation:
- Always validate phonon dispersions
- Check elastic constants
- Test thermal expansion
- Compare energy/force predictions

### Workflow Configuration:
- Start with small test systems
- Use appropriate DFT settings
- Monitor computational costs
- Iterate on training set

## Community and Support
- Open-source BSD-3-Clause
- Active development (AutoAtML team)
- Integration with atomate2
- Growing documentation
- Modern workflow approach

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/autoatml/autoplex
2. Documentation: https://autoatml.github.io/autoplex/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: OPEN (GitHub, BSD-3)
- Active development
