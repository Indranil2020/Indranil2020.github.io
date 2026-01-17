# DeePTB

## Official Resources
- Homepage: https://github.com/deepmodeling/DeePTB
- Source Repository: https://github.com/deepmodeling/DeePTB
- License: GNU Lesser General Public License v3.0

## Overview
DeePTB (Deep Learning for Tight-Binding) is a powerful software package that creates highly accurate tight-binding (TB) Hamiltonians using deep learning techniques. While traditional TB models rely on analytical Slater-Koster rules or manual parameterization, DeePTB trains deep neural networks to predict the Hamiltonian matrix elements directly from the atomic structure. This allows it to achieve *ab initio* accuracy (matching DFT bands and energies) while retaining the low computational cost of tight-binding, enabling simulation of large-scale systems with quantum accuracy.

**Scientific domain**: Machine Learning, Tight-Binding, Electronic Structure
**Target user community**: Materials scientists, Physicists studying large/disordered systems

## Theoretical Methods
- **Deep Neural Networks (DNN)**: Mapping structure to Hamiltonian.
- **Refined Slater-Koster (DeePTB-SK)**: Using NN to correct/modulate SK parameters.
- **Equivariant Neural Networks (DeePTB-E3)**: Directly predicting Hamiltonian blocks using E3-equivariant representations (no scaling laws assumed).
- **Spectral Properties**: Calculation of bands, DOS, and Berry curvature.

## Capabilities (CRITICAL)
- **Hamiltonian Prediction**: Fast generation of $H(k)$ for any configuration.
- **Band Structure**: Highly accurate reproduction of DFT bands.
- **Electronic Properties**: Total energy, forces, atomic charges.
- **Transferability**: Models trained on small supercells transfer to large systems/defects.
- **Efficiency**: Orders of magnitude faster than DFT.

## Key Strengths

### Ab Initio Accuracy:
- Unlike standard TB, DeePTB captures subtle environmental effects and hybridization shifts.
- Can describe complex bonding (halides, oxides) difficult for simple SK models.

### Symmetry Preserving:
- **DeePTB-E3** guarantees that the predicted Hamiltonian transforms correctly under rotation and translation, ensuring physical validity without data augmentation.

### Efficient Workflow:
- Part of the DeepModeling ecosystem (DeePMD), allowing seamless integration with MD workflows.

## Inputs & Outputs
- **Inputs**:
  - Atomic config (POSCAR/XYZ).
  - Training Data: DFT Hamiltonians/Eigenvalues (from VASP/ABACUS/OpenMX).
- **Outputs**:
  - Predicted TB Hamiltonian (sparse format).
  - Band structures.
  - Electronic Free Energy and Forces.

## Interfaces & Ecosystem
- **DeepModeling**: Compatible with DeePMD-kit workflows.
- **DFT Codes**: Interfaces to read Hamiltonians from VASP, ABACUS, and OpenMX.
- **Wannier90**: Can interface with Wannierized Hamiltonians.

## Advanced Features
- **Force Training**: Can train on forces to enable stable Molecular Dynamics.
- **Self-Consistency**: Can be coupled with a self-consistent density loop (though primarily used in non-SCF mode for bands).

## Performance Characteristics
- **Speed**: Hamiltonian generation is fast (inference time); Diagonalization cost dominates for large systems.
- **Scalability**: Linear scaling generation; Eigensolver depends on method.
- **Training**: Requires GPU for efficient training of E3 networks.

## Computational Cost
- **Inference**: Very Low (milliseconds for small cells).
- **Training**: Moderate (requires sufficient DFT data).

## Limitations & Known Constraints
- **Data Dependence**: The model is only as good as the DFT data it saw; extrapolation to unknown phases can be physical (due to E3) but inaccurate.
- **Hamiltonian Size**: Requires defining a basis set size (e.g., orbitals per atom) consistent with reference.

## Comparison with Other Codes
- **vs DFTB+**: DFTB+ uses a fixed physical model (2-center integrals); DeePTB uses a flexible Neural Network model (multi-center effects). DeePTB is more accurate but requires training.
- **vs Wannier90**: Wannier90 is a post-processing tool for one structure; DeePTB *learns* the function $H(R)$ to predict new structures.
- **vs SchNetPack**: SchNetPack predicts scalars (Energy) or vectors (Forces); DeePTB predicts Matrices (Hamiltonian).
- **Unique strength**: The State-of-the-Art implementation of Equivariant Neural Networks for electronic Hamiltonians.

## Application Areas
- **Disordered Systems**: Alloys, amorphous semiconductors.
- **Defects**: Electronic levels of vacancies/dopants in large supercells.
- **Twisted Bilayers**: Moiré physics requiring huge unit cells (thousands of atoms).
- **Finite Temperature**: Band structure renormalization due to thermal vibrations.

## Best Practices
- **Data Diversity**: Include distorted structures in training to learn bond-length dependence.
- **Check Symmetry**: Use DeePTB-E3 for crystals to ensure band degeneracies are respected.
- **Basis Consistency**: Ensure all DFT training data uses the exact same orbital projection/basis definition.

## Community and Support
- **GitHub**: Active issue tracking.
- **DeepModeling Community**: Large user base in China and global.

## Verification & Sources
**Primary sources**:
1. Repository: https://github.com/deepmodeling/DeePTB
2. Publication: "DeePTB: A deep learning package for tight-binding hamiltonians", *Phys. Rev. B* (etc.).

**Verification status**: ✅ VERIFIED
- Source code: OPEN (LGPLv3)
- Maturity: Research grade, rapidly adopting.
