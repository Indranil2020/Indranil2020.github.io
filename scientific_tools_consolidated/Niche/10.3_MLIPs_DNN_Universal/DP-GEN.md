# DP-GEN

## Official Resources
- Source Repository: https://github.com/deepmodeling/dpgen
- Documentation: https://docs.deepmodeling.com/projects/dpgen/
- License: Open source (LGPL-3.0)

## Overview
**DP-GEN** (Deep Potential GENerator) is an active learning workflow for generating deep learning interatomic potentials. It automates the cycle of training, exploration, and labeling to build high-quality DeePMD potentials with minimal human intervention.

**Scientific domain**: Active learning workflow for DeePMD potential generation  
**Target user community**: Researchers building DeePMD potentials via active learning

## Theoretical Methods
- Active learning loop (train-explore-label)
- Model disagreement-based selection
- Multiple DeePMD model training
- DFT labeling of selected structures
- Convergence monitoring

## Capabilities (CRITICAL)
- Automated active learning
- Multiple DFT code support (VASP, QE, ABACUS, Gaussian)
- Model disagreement selection
- Convergence monitoring
- HPC job management

**Sources**: GitHub repository, Nat. Commun. 11, 5713 (2020)

## Key Strengths

### Automated:
- Train-explore-label loop
- Model disagreement selection
- Convergence monitoring
- Minimal human intervention

### Multi-Code:
- VASP, QE, ABACUS for DFT
- Gaussian, CP2K for molecules
- DeePMD-kit for training
- LAMMPS for exploration

### Production:
- Well-tested at scale
- HPC job management
- Slurm/PBS support
- Checkpoint/restart

## Inputs & Outputs
- **Input formats**: Initial structures, DFT parameters
- **Output data types**: Trained DeePMD models, training data

## Interfaces & Ecosystem
- **DeePMD-kit**: Training
- **LAMMPS**: Exploration
- **VASP/QE/ABACUS**: DFT labeling
- **Python**: Core

## Performance Characteristics
- **Speed**: Limited by DFT calls
- **Accuracy**: Near-DFT (converged)
- **System size**: Any
- **Automation**: Full

## Computational Cost
- **DFT calls**: Expensive (but minimized)
- **Training**: Hours per iteration
- **Total**: Days to weeks

## Limitations & Known Constraints
- **DFT cost**: Still need DFT calculations
- **DeePMD only**: Only for DeePMD models
- **Complex setup**: Multi-step workflow
- **LGPL license**: Copyleft

## Comparison with Other Codes
- **vs FLARE**: DP-GEN is batch, FLARE is on-the-fly
- **vs Metatrain**: DP-GEN is active learning, Metatrain is training
- **vs manual fitting**: DP-GEN is automated
- **Unique strength**: Automated active learning workflow for DeePMD with multi-DFT-code support

## Application Areas

### Potential Development:
- Automated MLIP construction
- High-throughput potential generation
- Multi-component systems
- Phase diagram exploration

### Production:
- Converged DeePMD potentials
- Published potential generation
- Multi-fidelity training

## Best Practices
- Start with small training set
- Use model disagreement for selection
- Monitor convergence metrics
- Use HPC for DFT calculations

## Community and Support
- Open source (LGPL-3.0)
- DeepModeling community
- Comprehensive documentation
- Published in Nature Communications

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/deepmodeling/dpgen

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Automated active learning workflow for DeePMD with multi-DFT-code support
