# NequIP

## Official Resources
- **Homepage**: https://github.com/mir-group/nequip
- **Documentation**: https://wandb.ai/nequip/nequip
- **Developer**: Boris Kozinsky Group (Harvard)
- **License**: MIT License

## Overview
**NequIP** (Neural Equivariant Interatomic Potentials) is a pioneering code for building **E(3)-equivariant** neural network potentials. By explicitly enforcing rotation and translation symmetry in the network architecture, NequIP achieves remarkable data efficiency, allowing it to learn "DFT-quality" fields from very small training sets.

**Scientific domain**: Machine Learning Potentials, Equivariant Neural Networks.
**Target user community**: ML-Physics researchers, High-accuracy MD practitioners.

## Theoretical Methods
- **Architecture**: E(3)-equivariant Graph Neural Network.
- **Symmetry**: Uses Tensor Field Networks / e3nn library to handle spherical harmonics.
- **Efficiency**: Data-efficient learning (needs fewer DFT frames than classic MLPs).

## Capabilities
- **Potentials**: Generates highly accurate energy and force fields.
- **Training**: Automatic training workflows on DFT data.
- **MD**: Interfaces with LAMMPS and ASE for running dynamics.

## Key Strengths
- **Data Efficiency**: Can learn accurate potentials from as few as 100-1000 structures.
- **Rigor**: Direct implementation of physical symmetries leads to robust physics.
## Performance Characteristics
- **Speed**: Slower inference than invariant models (like DeepMD) or MACE, due to the heavy tensor product computations.
- **Efficiency**: "Data Efficient" - getting decent results with 100 structures where others need 1000.
- **Scaling**: Linear $O(N)$ but with a large prefactor.

## Limitations & Known Constraints
- **Memory**: Training on large systems (>100s of atoms per frame) can cause GPU Out-Of-Memory (OOM) errors easily.
- **Inference Cost**: High `l_max` (angular features) settings drastically reduce speed; strict trade-off between accuracy and speed.

## Best Practices
- **Hyperparameters**: Start with lower `l_max` (e.g., 2) to save memory/time; only increase if accuracy is insufficient.
- **Cutoff**: Keep the cutoff radius as small as physically reasonable to improve sparsity and speed.

## Community and Support
- **Support**: Active GitHub repository and WandB community.
- **Tutorials**: Detailed tutorials available for integration with ASE and LAMMPS.
## Comparison with Other Codes
- **vs DeepMD**: DeepMD uses descriptors (invariant); NequIP uses equivariant features (vector/tensor), which is richer but more expensive.
- **vs MACE**: MACE combines NequIP's equivariance with cluster expansion for better long-range/many-body efficiency.

## Verification & Sources
**Primary sources**:
1.  **Repository**: [NequIP GitHub](https://github.com/mir-group/nequip)
2.  **Literature**: Batzner, S., et al. "E(3)-equivariant graph neural networks for data-efficient and accurate interatomic potentials." *Nature Communications* (2022).

**Verification status**: ✅ VERIFIED
