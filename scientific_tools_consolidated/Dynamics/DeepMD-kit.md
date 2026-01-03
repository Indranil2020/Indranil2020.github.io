# DeepMD-kit

## Official Resources
- Homepage: http://www.deepmd.org/
- Documentation: https://docs.deepmd.org/
- Source Repository: https://github.com/deepmodeling/deepmd-kit
- License: GNU Lesser General Public License v3.0

## Overview
DeepMD-kit is a deep learning package for many-body potential energy representation and molecular dynamics. It allows users to train a deep neural network potential from ab initio data (DFT) and then use it to perform molecular dynamics simulations with ab initio accuracy but at a cost comparable to classical empirical potentials. It is a core component of the DeepModeling ecosystem.

**Scientific domain**: Machine learning potentials, deep learning, molecular dynamics  
**Target user community**: Computational chemists, materials scientists, ML-physics researchers

## Theoretical Methods
- Deep Potential Molecular Dynamics (DPMD)
- Deep Neural Networks (DNN) for PES
- Local atomic environment descriptors
- End-to-end symmetry preserving architecture
- Smoothness and continuity of PES
- Active learning (DP-GEN)

## Capabilities (CRITICAL)
- Training deep potentials from DFT data (VASP, QE, CP2K, etc.)
- Highly efficient MD interface with LAMMPS
- Accuracy comparable to DFT (typically < 1 meV/atom error)
- Linear scaling O(N) with system size
- GPU acceleration (CUDA/ROCm)
- Active learning workflow support
- Model compression for faster inference

**Sources**: DeepMD-kit documentation, Comp. Phys. Comm. 228, 178 (2018)

## Inputs & Outputs
- **Input formats**: Training data (coordinates, forces, energies, virials) in NumPy/HDF5 format
- **Output data types**: Trained model (.pb), Training logs, Validation metrics

## Interfaces & Ecosystem
- **LAMMPS**: Primary MD engine interface
- **i-PI**: Socket interface
- **ASE**: Python calculator interface
- **DP-GEN**: Active learning workflow manager
- **GROMACS/OpenMM**: Experimental/Third-party interfaces

## Workflow and Usage
1. Data generation: Run ab-initio calculations (VASP/QE)
2. Data prep: Convert to DeepMD format
3. Train: `dp train input.json`
4. Freeze: `dp freeze -o graph.pb`
5. Run MD: Use `pair_style deepmd` in LAMMPS

## Performance Characteristics
- Millions of atoms/day on GPUs
- Significant speedup over DFT (1000x-10000x)
- Slower than simple empirical potentials but much more accurate
- Optimized for NVIDIA GPUs

## Application Areas
- Water and ice phase diagrams
- High-entropy alloys
- Chemical reactions and catalysis
- Battery materials (electrolytes)
- Warm dense matter

## Community and Support
- Open-source (LGPL v3)
- Very active GitHub community
- Developer conferences (DeepModeling)
- Detailed documentation and tutorials

## Verification & Sources
**Primary sources**:
1. Homepage: http://www.deepmd.org/
2. GitHub: https://github.com/deepmodeling/deepmd-kit
3. Publication: Wang et al., Comp. Phys. Comm. 228, 178 (2018)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE (DeepModeling Community)
- Applications: Deep learning potentials, DPMD, LAMMPS interface, active learning
