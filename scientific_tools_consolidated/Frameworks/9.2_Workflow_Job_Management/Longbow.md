# Longbow

## Official Resources
- Source Repository: https://github.com/CCPBioSim/Longbow
- Documentation: https://longbow.readthedocs.io/
- PyPI: https://pypi.org/project/longbow/
- License: Open source (BSD-3)

## Overview
**Longbow** is a tool for automating simulations on remote HPC machines. It is designed to mimic the normal way an application is run locally but sends simulations to powerful machines, supporting LAMMPS, GROMACS, NAMD, AMBER, and other simulation packages.

**Scientific domain**: HPC job automation, remote simulation management  
**Target user community**: Researchers running MD/DFT simulations on remote HPC clusters

## Theoretical Methods
- Remote HPC job submission
- Local-like execution model
- File staging and retrieval
- Multi-scheduler support (Slurm, SGE, Torque)
- Job monitoring and recovery

## Capabilities (CRITICAL)
- Remote HPC job submission
- Local-like execution interface
- File staging (upload/download)
- Multiple scheduler support
- Job monitoring
- Error recovery and resubmission

**Sources**: GitHub repository, ReadTheDocs

## Key Strengths

### Remote Execution:
- Local-like interface
- Automatic file staging
- Multiple schedulers
- Job monitoring

### Multi-Code:
- LAMMPS, GROMACS, NAMD, AMBER
- Any command-line application
- Flexible configuration
- Plugin architecture

### Easy to Use:
- Simple configuration
- One-command execution
- Automatic recovery
- Logging and debugging

## Inputs & Outputs
- **Input formats**: Local simulation files, HPC configuration
- **Output data types**: Retrieved results, log files

## Interfaces & Ecosystem
- **Slurm**: Job scheduler
- **SGE/Torque**: Alternative schedulers
- **SSH**: Remote access
- **Python**: Core language

## Performance Characteristics
- **Speed**: Job management (fast)
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **Framework**: Negligible
- **Simulations**: Hours (on HPC)

## Limitations & Known Constraints
- **SSH required**: Need HPC access
- **Configuration**: Per-cluster setup
- **No provenance**: No data tracking
- **Not DFT-specific**: General simulation tool

## Comparison with Other Codes
- **vs AiiDA**: Longbow is simpler, AiiDA has provenance
- **vs FireWorks**: Longbow is remote-focused, FireWorks is database-backed
- **vs Parsl**: Longbow is local-like, Parsl is Python-native
- **Unique strength**: Local-like remote HPC execution with automatic file staging and multi-scheduler support

## Application Areas

### HPC Simulations:
- Remote MD simulations
- DFT calculations on clusters
- Batch job management
- Multi-node execution

### Multi-Code:
- LAMMPS on HPC
- GROMACS on HPC
- Any CLI application remotely

## Best Practices

### Setup:
- Configure HPC connection
- Set up scheduler template
- Test with simple job
- Use recovery mode

## Community and Support
- Open source (BSD-3)
- PyPI installable
- CCP BioSim maintained
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/CCPBioSim/Longbow

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Local-like remote HPC execution with automatic file staging and multi-scheduler support
