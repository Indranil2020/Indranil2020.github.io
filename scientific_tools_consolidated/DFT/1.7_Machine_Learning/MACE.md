# MACE

## Official Resources
- **Homepage**: https://mace-docs.readthedocs.io/
- **Source Repository**: https://github.com/ACEsuit/mace
- **Developer**: Gabor Csanyi Group (University of Cambridge) / Ilyes Batatia
- **License**: MIT License

## Overview
**MACE** (Machine learning Atomic Cluster Expansion) represents a new generation of **Machine Learning Interatomic Potentials (MLIPs)**. It unifies the accuracy of body-ordered equivariant features (like NequIP) with the efficiency of message passing. While technically a potential, it is so closely tied to reproducing DFT quality at scale that it is essential for modern "ML-Enhanced DFT" workflows.

**Scientific domain**: Machine Learning Potentials, Molecular Dynamics, Materials Discovery.
**Target user community**: Researchers needing DFT-level accuracy for millions of atoms.

## Theoretical Methods
- **Approach**: Multi-body Atomic Cluster Expansion (MACE).
- **Architecture**: Higher-order equivariant message passing neural networks.
- **Training**: Trained on DFT forces and energies.

## Capabilities
- **Scaling**: Linear scaling with system size.
- **Accuracy**: State-of-the-art performance on benchmarks (e.g., Matbench Discovery).
- **Foundation Models**: "MACE-MP-0" provides a general-purpose potential for 89 elements.

## Key Strengths
- **Generalization**: Excellent extrapolation to unseen structures.
- **Speed**: Faster than pure equivariant GNNs (like NequIP) while maintaining high accuracy.
- **Ecosystem**: Part of the ACEsuite of tools.
## Performance Characteristics
- **Pareto Optimal**: Often cited as having the best trade-off between accuracy and computational cost on modern benchmarks (e.g., Matbench Discovery).
- **Inference**: Significantly faster than NequIP for similar accuracy levels due to efficient message passing architecture.

## Limitations & Known Constraints
- **Training**: Memory intensive during training, though less so than full rank-3 equivariant models.
- **Scope**: Like all potentials, it cannot natively predict electronic properties (band gaps, DOS) unless specifically trained on them as auxiliary targets (unlike DeepH).

## Best Practices
- **Foundation Models**: Try the pre-trained "MACE-MP-0" model first before training your own; it covers 89 elements and is a strong baseline.
- **Hardware**: GPU acceleration is highly recommended for both training and MD inference.

## Community and Support
- **Support**: Active development by the ACEsuit organization on GitHub.
## Comparison with Other Codes
- **vs DeepH**: MACE predicts energy/forces (for MD); DeepH predicts the electronic Hamiltonian (for band structure).
- **vs NequIP**: MACE is generally faster and scales to more neighbors thanks to the cluster expansion formalism.

## Verification & Sources
**Primary sources**:
1.  **Repository**: [MACE GitHub](https://github.com/ACEsuit/mace)
2.  **Literature**: Batatia, I., et al. "MACE: Higher order equivariant message passing neural networks..." *NeurIPS* (2022).

**Verification status**: ✅ VERIFIED
