# ALIGNN-FF

## Official Resources
- Source Repository: https://github.com/usnistgov/alignn
- Documentation: https://alignn.readthedocs.io/
- PyPI: https://pypi.org/project/alignn/
- License: Open source (MIT)

## Overview
**ALIGNN-FF** (Atomistic Line Graph Neural Network Force Field) is a GNN potential that uses line graphs to capture bond angles and dihedrals. It covers 5-118 elements and achieves competitive accuracy for both property prediction and MD simulations.

**Scientific domain**: Line graph GNN potential with bond angle/dihedral awareness  
**Target user community**: Researchers needing GNN potential with explicit angular information

## Theoretical Methods
- Atomistic Line Graph Neural Network
- Bond angle and dihedral features
- Line graph message passing
- 5-118 elements coverage
- Property prediction + force field

## Capabilities (CRITICAL)
- Line graph architecture
- 5-118 elements coverage
- Property prediction (formation energy, bandgap)
- Force field for MD
- JARVIS integration

**Sources**: GitHub repository, npj Comput. Mater. 7, 185 (2021)

## Key Strengths

### Line Graph:
- Explicit bond angles
- Dihedral information
- Better angular features
- Physically motivated

### Broad Coverage:
- 5-118 elements
- Property prediction
- Force field
- JARVIS dataset integration

### NIST:
- Government maintained
- Well documented
- JARVIS ecosystem
- Reproducible

## Inputs & Outputs
- **Input formats**: Structures (ASE/cif)
- **Output data types**: Properties, energies, forces

## Interfaces & Ecosystem
- **JARVIS**: Dataset and tools
- **ASE**: Calculator
- **PyTorch**: Backend
- **Python**: Core

## Performance Characteristics
- **Speed**: Fast (GPU)
- **Accuracy**: Competitive
- **System size**: 1-10000+ atoms
- **Automation**: Full

## Computational Cost
- **Training**: Hours on GPU
- **MD**: Fast

## Limitations & Known Constraints
- **Not universal pretrained**: Needs training
- **Line graph overhead**: Slightly slower
- **NIST focus**: JARVIS benchmarks

## Comparison with Other Codes
- **vs CGCNN**: ALIGNN has line graph, CGCNN does not
- **vs M3GNet**: ALIGNN is line graph, M3GNet is 3-body
- **Unique strength**: Line graph GNN with explicit bond angle/dihedral features covering 5-118 elements

## Application Areas

### Property Prediction:
- Formation energy
- Bandgap
- Elastic properties
- JARVIS benchmarks

### MD:
- Materials simulation
- Structure relaxation
- Energy evaluation

## Best Practices
- Use JARVIS training data
- Start with pretrained models
- Fine-tune for target systems

## Community and Support
- Open source (MIT)
- NIST maintained
- PyPI installable
- ReadTheDocs documentation

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/usnistgov/alignn

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Source code: ACCESSIBLE (GitHub)
- PyPI: AVAILABLE
- Specialized strength: Line graph GNN with explicit bond angle/dihedral features covering 5-118 elements
