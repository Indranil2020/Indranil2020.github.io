# AlphaCrystal / AlphaCrystal-II

## Overview
AlphaCrystal is a contact map-based deep learning algorithm for crystal structure prediction, inspired by AlphaFold's approach to protein structure prediction. AlphaCrystal-II extends this with distance matrix-based prediction.

## Theoretical Basis
- Contact map prediction
- Distance matrix prediction (AlphaCrystal-II)
- Deep residual neural networks
- Structure reconstruction from distances
- Composition-to-structure mapping

## Key Capabilities
- Crystal structure prediction from composition
- Contact/distance map prediction
- Deep learning approach
- No DFT during prediction
- Fast inference

**Sources**: ACS Omega (2023), arXiv:2404.04810

## Key Strengths

### Methodology:
- AlphaFold-inspired
- Contact/distance maps
- Deep learning

### Speed:
- Fast inference
- No DFT required
- Large-scale applicable

### Innovation:
- Novel approach
- Protein-inspired
- Distance-based

## Inputs & Outputs
- **Input formats**: Chemical composition
- **Output data types**: Predicted crystal structures, distance maps

## Interfaces & Ecosystem
- **Framework**: TensorFlow/PyTorch
- **Tools**: MLatticeABC, Cryspnet
- **Validation**: DFT codes

## Workflow and Usage
1. Input chemical composition
2. Predict distance/contact map
3. Reconstruct 3D structure
4. Refine structure
5. Validate with DFT

## Performance Characteristics
- Fast prediction
- GPU-accelerated
- Good for screening

## Computational Cost
- Prediction: fast
- No DFT required
- Validation: DFT

## Best Practices
- Use for initial screening
- Validate top predictions
- Consider multiple candidates
- Check structural validity

## Limitations & Known Constraints
- Training data dependent
- May miss novel structures
- Requires validation

## Application Areas
- High-throughput screening
- Initial structure guessing
- Materials discovery
- Composition-to-structure

## Comparison with Other Codes
- **vs DiffCSP**: Different approach (maps vs diffusion)
- **vs CSPML**: Different methodology
- **Unique strength**: AlphaFold-inspired, distance maps

## Community and Support
- Open-source (GitHub)
- USC development
- Published methodology

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/usccolumbia/AlphaCrystal
2. AlphaCrystal: ACS Omega (2023)
3. AlphaCrystal-II: arXiv:2404.04810

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE (GitHub)
- Source: OPEN
- Development: ACTIVE
- Applications: Crystal structure prediction
