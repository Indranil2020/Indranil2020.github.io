# ACEMD

## Official Resources
- Homepage: https://www.acellera.com/acemd
- Documentation: https://software.acellera.com/acemd/
- License: Commercial (free for academics)

## Overview
ACEMD is a high-performance molecular dynamics engine specifically designed for GPU acceleration. Built on OpenMM, it provides an optimized and user-friendly interface for biomolecular simulations with exceptional speed on NVIDIA GPUs.

**Scientific domain**: Biomolecular simulations, drug discovery, GPU-accelerated MD  
**Target user community**: Pharmaceutical researchers, computational biologists

## Theoretical Methods
- Classical molecular dynamics
- Langevin dynamics
- Multiple force fields (AMBER, CHARMM)
- PME electrostatics
- Implicit solvent (GBSA)
- Replica exchange

## Capabilities (CRITICAL)
- Ultra-fast GPU MD simulations
- AMBER/CHARMM force field support
- Implicit and explicit solvent
- Replica exchange MD
- Metadynamics integration
- HTMD workflow integration

## Key Strengths

### GPU Performance:
- Optimized for NVIDIA GPUs
- Microsecond timescales routine
- Multi-GPU support
- Exceptional speed

### Integration:
- HTMD workflow
- PlayMolecule platform
- Automated setup

## Inputs & Outputs
- **Input formats**:
  - PDB structures
  - AMBER prmtop
  - CHARMM PSF
  
- **Output data types**:
  - XTC trajectories
  - DCD trajectories
  - Restart files

## Interfaces & Ecosystem
- **HTMD**: High-throughput MD
- **PlayMolecule**: Web platform
- **OpenMM**: Backend engine

## Advanced Features
- **GPU optimization**: CUDA-optimized kernels
- **Metadynamics**: Enhanced sampling
- **Replica exchange**: REMD support
- **Adaptive sampling**: With HTMD
- **Free energy**: Alchemical methods

## Performance Characteristics
- Among fastest GPU MD codes
- Optimized memory usage
- Excellent for long simulations
- Multi-GPU scaling

## Computational Cost
- GPU provides 100x+ speedup
- Microseconds per day achievable
- Efficient for large systems
- Overall: Industry-leading GPU performance

## Best Practices
- Use latest NVIDIA GPUs
- Enable mixed precision
- Use HTMD for workflows
- Validate force field choice

## Limitations & Known Constraints
- Commercial license
- NVIDIA GPU required
- Less flexible than OpenMM
- Biomolecular focus

## Application Areas
- Drug discovery
- Protein dynamics
- Membrane simulations
- Long timescale dynamics
- Adaptive sampling

## Community and Support
- Commercial support (Acellera)
- Documentation
- Tutorials
- Email support

## Verification & Sources
**Primary sources**:
1. Website: https://www.acellera.com/acemd
2. M. Harvey et al., J. Chem. Theory Comput. 5, 1632 (2009)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Commercial with academic license
- Well-documented
