# GRACE

## Official Resources
- Source Repository: https://github.com/ICAMS/grace-tensorpotential
- License: Open source

## Overview
**GRACE** (Graph Atomic Cluster Expansion) is a foundation model for interatomic potentials combining ACE features with graph neural network architecture. The GRACE-2L model (15.3M parameters) achieves state-of-art accuracy on OMat24 and MPTraj datasets.

**Scientific domain**: Graph ACE foundation model for universal potentials  
**Target user community**: Researchers needing high-accuracy universal potential with ACE features

## Theoretical Methods
- Graph Atomic Cluster Expansion
- ACE polynomial basis + GNN
- Foundation model approach
- 15.3M parameters (2L model)
- Active learning support

## Capabilities (CRITICAL)
- Foundation model (89+ elements)
- ACE + GNN architecture
- 15.3M parameters (2L)
- OMat24/MPTraj training
- Active learning integration

**Sources**: GitHub repository, npj Comput. Mater. 11, 50 (2025)

## Key Strengths

### ACE + GNN:
- Systematic ACE completeness
- GNN flexibility
- Best of both worlds
- Physically motivated

### Foundation:
- 89+ elements
- Pretrained model
- Fine-tuning support
- Active learning

## Inputs & Outputs
- **Input formats**: Structures (ASE)
- **Output data types**: Energies, forces, stresses

## Interfaces & Ecosystem
- **ASE**: Calculator
- **LAMMPS**: MD engine
- **Python**: Core

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: State-of-art
- **System size**: 1-10000+ atoms

## Computational Cost
- **MD**: ~1000x faster than DFT
- **Fine-tuning**: Hours

## Limitations & Known Constraints
- **New project**: Still maturing
- **GPU required**: For production
- **ICAMS maintained**: Small team

## Comparison with Other Codes
- **vs MACE**: GRACE has ACE basis, MACE is pure equivariant
- **vs ACE1pack**: GRACE is GNN+ACE, ACE1pack is pure ACE
- **Unique strength**: Graph ACE foundation model combining systematic ACE completeness with GNN flexibility

## Application Areas

### Universal MD:
- High-accuracy MD
- Phase stability
- Mechanical properties
- Active learning workflows

## Best Practices
- Use GRACE-2L for best accuracy
- Fine-tune for specific chemistry
- Validate with DFT benchmarks

## Community and Support
- Open source
- ICAMS maintained
- Published in npj Computational Materials

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/ICAMS/grace-tensorpotential

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- Specialized strength: Graph ACE foundation model combining systematic ACE completeness with GNN flexibility
