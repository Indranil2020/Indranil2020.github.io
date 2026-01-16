# DeepH (Deep Hamiltonian)

## Official Resources
- **Homepage**: https://deeph-pack.readthedocs.io/
- **Source Repository**: https://github.com/mzjb/DeepH-pack
- **Developers**: H. Li Group (Tsinghua University)
- **License**: MIT License

## Overview
**DeepH** is a state-of-the-art framework for **Machine Learning Enhanced DFT**. It bypasses the computationally expensive self-consistent field (SCF) iterations of traditional DFT by learning a mapping from atomic structures to the DFT Hamiltonian. Utilizing E(3)-equivariant neural networks, it can predict accurate Hamiltonians for large-scale material systems, generalizing from small unit cells to huge supercells with $O(N)$ inference cost.

**Scientific domain**: Material Science, Machine Learning, Large-scale Electronic Structure.
**Target user community**: Researchers studying twisted bilayers, large supercells, and complex material interfaces.

## Theoretical Methods
- **Approach**: Supervised learning of the DFT Hamiltonian matrix.
- **Architecture**: E(3)-equivariant graph neural networks (e.g., DeepH-E3).
- **Basis**: Compatible with localized basis sets (PAO/LCAO).
- **Workflow**: Train on small system DFT data -> Predict on large system -> Diagonalize/Calculate properties.

## Capabilities
- **Prediction**: Direct prediction of Hamiltonian and Overlap matrices.
- **Scalability**: Capable of handling 10,000+ atom systems where conventional DFT is prohibitive.
- **Interfaces**: Works with OpenMX, ABACUS, FHI-aims, and SIESTA.
- **Applications**: Moiré superlattices, defects, amorphous systems.

## Key Strengths
- **Efficiency**: Reduces calculation time from hours/days (DFT) to minutes (Inference).
- **Generalizability**: A model trained on small perturbations can predict large, unrelaxed structures.
- **Accuracy**: Retains near-DFT methodology accuracy (meV scale) unlike simpler ML potentials which only give energy/forces.
## Performance Characteristics
- **Speed**: Inference is 100-1000x faster than corresponding DFT calculations.
- **Scaling**: Linear scaling $O(N)$ with system size for inference.

## Limitations & Known Constraints
- **Transferability**: While better than simple ML potentials, preventing "catastrophic forgetting" or ensuring transferability to completely unseen chemical environments remans a challenge.
- **Training Cost**: Generating the DFT training database is the bottleneck; the ML model can only be as good as the DFT reference (inherits DFT errors).
- **Complexity**: Predicting a Hamiltonian matrix is dimensionally more complex than predicting a scalar energy, requiring significant GPU memory for training.

## Best Practices
- **Dataset**: Ensure the training set covers the structural deformations expected in the production run.
- **Basis**: Use a consistent localized basis (e.g., OpenMX PAOs) for both training and inference; you cannot mix basis sets.

## Community and Support
- **Support**: GitHub Issues and documentation on ReadTheDocs.
## Comparison with Other Codes
- **vs Conventional DFT**: DeepH is not a replacement but an accelerator/extender. It needs conventional DFT for training data.
- **vs ML Potentials (NequIP/MACE)**: Potentials predict energy/forces for MD; DeepH predicts the *electronic Hamiltonian*, effectively giving access to band structures and wavefunctions of large systems.

## Verification & Sources
**Primary sources**:
1.  **Repository**: [DeepH-pack GitHub](https://github.com/mzjb/DeepH-pack)
2.  **Literature**: Li, H., et al. "Deep-learning density functional theory Hamiltonian for large-scale electronic-structure calculation." *Nature Computational Science* (2022).

**Verification status**: ✅ VERIFIED
