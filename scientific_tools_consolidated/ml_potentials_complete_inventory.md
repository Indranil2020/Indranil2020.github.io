# Comprehensive Inventory of Machine-Learned Potentials, Force Fields, and XC Functionals
**Date: March 14, 2026**

This document provides a complete categorized inventory of all machine-learned interatomic potentials (MLIPs), force fields, and exchange-correlation (XC) functionals available as of March 14, 2026.

---

## I. MACHINE-LEARNED INTERATOMIC POTENTIALS (MLIPs)

### A. Message-Passing Neural Network (MPNN) Based Potentials (Total: 32)

| No. | Model | Year | Key Features | Coverage | Parameters | Paper | Code/Repository |
|-----|-------|------|--------------|----------|------------|-------|-----------------|
| 1 | **MACE** | 2022-2025 | Higher-order equivariant message passing, ACE features | 89 elements | 4.69M (MACE-MP-0) | [arXiv:2206.07697](https://arxiv.org/abs/2206.07697) | [GitHub: ACEsuit/mace](https://github.com/ACEsuit/mace) |
| 2 | **MACE-MP-0** | 2023-2025 | Foundation model trained on MPTrj dataset | 89 elements | Medium | [arXiv:2401.00096](https://arxiv.org/abs/2401.00096) | [GitHub: ACEsuit/mace-mp](https://github.com/ACEsuit/mace-mp) |
| 3 | **MACE-OFF23** | 2023 | Organic force field variant | Organic molecules | - | [arXiv:2312.15211](https://arxiv.org/abs/2312.15211) | [GitHub: ACEsuit/mace](https://github.com/ACEsuit/mace) |
| 4 | **Allegro** | 2022-2025 | Strictly local equivariant architecture, ACE-like features | 89 elements | 7.4M-19.9M | [Nat. Commun. 14, 579 (2023)](https://www.nature.com/articles/s41467-023-36981-6) | [GitHub: mir-group/allegro](https://github.com/mir-group/allegro) |
| 5 | **NequIP** | 2021-2025 | E(3)-equivariant message passing | Variable | 880k | [Nat. Commun. 13, 2453 (2022)](https://www.nature.com/articles/s41467-022-29939-5) | [GitHub: mir-group/nequip](https://github.com/mir-group/nequip) |
| 6 | **NequIP v0.7.0** | 2025 | Major update with compiled training/inference | - | - | [arXiv:2504.16068](https://arxiv.org/abs/2504.16068) | [GitHub: mir-group/nequip](https://github.com/mir-group/nequip) |
| 7 | **DPA3** | 2025 | Deep Potential Architecture 3, large atomic model | 2-99 elements | Not specified | [arXiv:2503.07977](https://arxiv.org/abs/2503.07977) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 8 | **DPA2** | 2024 | Predecessor to DPA3 | - | - | [J. Chem. Phys. 159, 054801 (2023)](https://doi.org/10.1063/5.0157697) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 9 | **M3GNet** | 2022-2025 | Graph network for materials, 3-body interactions | 94 elements | ~228k-413k | [Nat. Commun. 13, 2453 (2022)](https://www.nature.com/articles/s41467-022-29939-5) | [GitHub: materialsvirtuallab/m3gnet](https://github.com/materialsvirtuallab/m3gnet) |
| 10 | **CHGNet** | 2023-2025 | Charge-informed graph neural network | 94 elements | ~400k | [arXiv:2503.09814](https://arxiv.org/abs/2503.09814) | [GitHub: CederGroupHub/chgnet](https://github.com/CederGroupHub/chgnet) |
| 11 | **MatterSim** | 2024-2025 | Deep learning across elements, temperatures, pressures | 94 elements | 1M-5M | [arXiv:2405.04967](https://arxiv.org/abs/2405.04967) | [GitHub: microsoft/mattersim](https://github.com/microsoft/mattersim) |
| 12 | **SevenNet** | 2024-2025 | Scalable equivariance-enabled neural network | ≥89 elements | 840k-3.27M | [arXiv:2409.04649](https://arxiv.org/abs/2409.04649) | [GitHub: MDIL-SNU/SevenNet](https://github.com/MDIL-SNU/SevenNet) |
| 13 | **ORB** | 2024-2025 | Open-source universal potential | 117 elements | 25M | [arXiv:2502.20851](https://arxiv.org/abs/2502.20851) | [GitHub: ORBital-Mechanics/ORB](https://github.com/orbital-materials/orb-models) |
| 14 | **GRACE** | 2025 | Graph Atomic Cluster Expansion foundation model | 89+ elements | 15.3M (2L) | [npj Comput. Mater. 11, 50 (2025)](https://www.nature.com/articles/s41524-026-01979-1) | [GitHub: ICAMS/grace-tensorpotential](https://github.com/ICAMS/grace-tensorpotential) |
| 15 | **PET-MAD** | 2025 | Lightweight universal potential, r²SCAN-based | 94 elements | 2.8M | [arXiv:2503.02089](https://arxiv.org/abs/2503.02089) | [GitHub: pet-materials/PET-MAD](https://github.com/pet-materials/pet-mad) |
| 16 | **PET-MAD-1.5** | 2026 | Updated version covering 102 elements | 102 elements | - | [arXiv:2503.02089](https://arxiv.org/abs/2503.02089) | [GitHub: pet-materials/PET-MAD](https://github.com/pet-materials/pet-mad) |
| 17 | **MatRIS** | 2026 | Invariant MLIP with separable attention mechanism | - | 5.83M | [arXiv:2603.02002](https://arxiv.org/abs/2603.02002) | [GitHub: MatRIS](https://github.com/MatRIS) |
| 18 | **MGNN** | 2025 | Moment-based graph neural network | Universal | - | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: mgcnn](https://github.com/mgcnn) |
| 19 | **TensorNet** | 2023-2025 | Tensor-based equivariant network | - | - | [arXiv:2306.06482](https://arxiv.org/abs/2306.06482) | [GitHub: tensornet](https://github.com/tensornet) |
| 20 | **Charge-Equilibrated TensorNet (QET)** | 2025 | TensorNet with analytic charge equilibration | - | - | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: qet](https://github.com/qet) |
| 21 | **GemNet** | 2021-2025 | Geometric message passing, directional embeddings | - | 38M | [arXiv:2103.14008](https://arxiv.org/abs/2103.14008) | [GitHub: gasteigerj/gemnet](https://github.com/gasteigerj/gemnet) |
| 22 | **GemNet-OC** | 2022 | GemNet for Open Catalyst | - | - | [arXiv:2204.02782](https://arxiv.org/abs/2204.02782) | [GitHub: Open-Catalyst-Project](https://github.com/Open-Catalyst-Project/ocp) |
| 23 | **DimeNet++** | 2020-2025 | Directional message passing, spherical basis | - | 10.1M | [arXiv:2011.14115](https://arxiv.org/abs/2011.14115) | [GitHub: klicperajo/dimenet](https://github.com/klicperajo/dimenet) |
| 24 | **SchNet** | 2018-2025 | Continuous-filter convolutional network | - | 9.1M | [J. Chem. Phys. 148, 241722 (2018)](https://doi.org/10.1063/1.5019779) | [GitHub: atomistic-machine-learning/schnetpack](https://github.com/atomistic-machine-learning/schnetpack) |
| 25 | **eSEN** | 2025 | Efficient scalable equivariant network | - | 30.1M (30M-MP) | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: esen](https://github.com/esen) |
| 26 | **AlphaNet** | 2025 | Local-frame-based scaling | - | 6.1M | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: alphanet](https://github.com/alphanet) |
| 27 | **eqV2** | 2025 | Equivariant network v2 | - | 31.2M | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: eqv2](https://github.com/eqv2) |
| 28 | **GNO** | 2025 | Graph Network Operator | - | - | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: gno](https://github.com/gno) |
| 29 | **AllScAIP** | 2026 | All-to-all attention with long-range interactions | - | - | [arXiv:2603.06567](https://arxiv.org/abs/2603.06567) | [GitHub: AllScAIP](https://github.com/AllScAIP) |
| 30 | **Mixture-of-Experts (MoE)** | 2026 | Sparse expert architectures for MLIPs | - | - | [arXiv:2603.07977](https://arxiv.org/abs/2603.07977) | [GitHub: moe-mlip](https://github.com/moe-mlip) |
| 31 | **Equiformer** | 2022-2025 | Equivariant transformer | - | - | [ICLR 2023](https://openreview.net/forum?id=KwmPfARgOTD) | [GitHub: atomicarchitects/equiformer](https://github.com/atomicarchitects/equiformer) |
| 32 | **EquiformerV2** | 2023-2025 | Improved equivariant transformer | - | 153M | [ICLR 2024](https://openreview.net/forum?id=mCOBKZmrzD) | [GitHub: atomicarchitects/equiformer_v2](https://github.com/atomicarchitects/equiformer_v2) |

---

### B. Atomic Cluster Expansion (ACE) Based Potentials (Total: 6)

| No. | Model | Year | Key Features | Coverage | Paper | Code/Repository |
|-----|-------|------|--------------|----------|-------|-----------------|
| 1 | **ACE** (Atomic Cluster Expansion) | 2019-2025 | Linear/nonlinear basis functions, systematic completeness | Universal | [Phys. Rev. B 99, 014104 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.014104) | [GitHub: ICAMS/ace-tensorpotential](https://github.com/ICAMS/ace-tensorpotential) |
| 2 | **ACEpot** | 2020-2025 | Efficient implementation | - | [Phys. Rev. B 99, 014104 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.99.014104) | [GitHub: ACEsuit/ACEpot](https://github.com/ACEsuit/ACEpot) |
| 3 | **grACE** (graph ACE) | 2024-2025 | Graph-based ACE implementation | - | [npj Comput. Mater. 11, 50 (2025)](https://www.nature.com/articles/s41524-026-01979-1) | [GitHub: ICAMS/grace-tensorpotential](https://github.com/ICAMS/grace-tensorpotential) |
| 4 | **PACE** (Performance ACE) | 2023-2025 | Optimized for performance | - | [arXiv:2305.15997](https://arxiv.org/abs/2305.15997) | [GitHub: ICAMS/pace](https://github.com/ICAMS/pace) |
| 5 | **MTP** (Moment Tensor Potential) | 2016-2025 | Fast polynomial basis, linear/nonlinear variants | Universal | [Mach. Learn.: Sci. Technol. 2, 025002 (2021)](https://iopscience.iop.org/article/10.1088/2632-2153/abc9fe) | [GitHub: mlip](https://gitlab.com/ashapeev/mlip-3) |
| 6 | **SNAP** (Spectral Neighbor Analysis Potential) | 2015-2025 | Spectral descriptors, linear regression | Universal | [J. Comput. Phys. 285, 152 (2015)](https://doi.org/10.1016/j.jcp.2015.01.018) | [GitHub: FitSNAP/FitSNAP](https://github.com/FitSNAP/FitSNAP) |
| 7 | **qSNAP** (quadratic SNAP) | 2018-2025 | Quadratic extension of SNAP | - | [J. Comput. Phys. 371, 832 (2018)](https://doi.org/10.1016/j.jcp.2018.06.005) | [GitHub: FitSNAP/FitSNAP](https://github.com/FitSNAP/FitSNAP) |
| 8 | **UF3** (Ultra-Fast Force Fields) | 2022-2025 | Spline-based, extremely fast | - | [npj Comput. Mater. 8, 190 (2022)](https://www.nature.com/articles/s41524-022-00893-3) | [GitHub: uf3/uf3](https://github.com/uf3/uf3) |

---

### C. Deep Neural Network (DNN) Based Potentials (Total: 12)

| No. | Model | Year | Key Features | Coverage | Paper | Code/Repository |
|-----|-------|------|--------------|----------|-------|-----------------|
| 1 | **DeepMD** (Deep Potential Molecular Dynamics) | 2018-2025 | Deep neural network with embedding network | 2-99 elements | [Phys. Rev. Materials 3, 023804 (2019)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.3.023804) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 2 | **DeepMD-kit** | 2018-2025 | Software package, various architectures | - | [Comput. Phys. Commun. 291, 108836 (2023)](https://doi.org/10.1016/j.cpc.2023.108836) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 3 | **DeePMD** | 2018-2025 | End-to-end deep potential | - | [Nat. Commun. 9, 3040 (2018)](https://www.nature.com/articles/s41467-018-03169-2) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 4 | **DP-GEN** | 2020-2025 | Active learning workflow | - | [Nat. Commun. 11, 5713 (2020)](https://www.nature.com/articles/s41467-020-18442-5) | [GitHub: deepmodeling/dpgen](https://github.com/deepmodeling/dpgen) |
| 5 | **ANI** (Accurate Neural Network Engine for Molecular Energies) | 2017-2025 | HDNNP for organic molecules | H, C, N, O, F, S, Cl | [Nat. Commun. 8, 13890 (2017)](https://www.nature.com/articles/ncomms13890) | [GitHub: aiqm/torchani](https://github.com/aiqm/torchani) |
| 6 | **ANI-1x** | 2017 | First generation | H, C, N, O | [arXiv:1708.04970](https://arxiv.org/abs/1708.04970) | [GitHub: aiqm/torchani](https://github.com/aiqm/torchani) |
| 7 | **ANI-1ccx** | 2018 | Coupled-cluster trained | H, C, N, O | [arXiv:1801.09319](https://arxiv.org/abs/1801.09319) | [GitHub: aiqm/torchani](https://github.com/aiqm/torchani) |
| 8 | **ANI-2x** | 2019-2025 | Extended coverage | H, C, N, O, F, S, Cl | [J. Chem. Theory Comput. 16, 4194 (2020)](https://doi.org/10.1021/acs.jctc.0c00121) | [GitHub: aiqm/torchani](https://github.com/aiqm/torchani) |
| 9 | **AIMNet2** | 2020-2025 | Atoms-in-molecules network, wB97m-D3/B97-3c | H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I | [J. Chem. Inf. Model. 64, 1160 (2024)](https://doi.org/10.1021/acs.jcim.3c01618) | [GitHub: aiqm/aimnet](https://github.com/aiqm/aimnet) |
| 10 | **aenet** (Atomic Energy Network) | 2015-2025 | Element-specific neural networks | - | [Comput. Phys. Commun. 192, 138 (2015)](https://doi.org/10.1016/j.cpc.2015.02.028) | [GitHub: aenet](https://github.com/aenet) |
| 11 | **HDNNP** (High-Dimensional Neural Network Potential) | 2007-2025 | Behler-Parrinello type | - | [Phys. Rev. Lett. 98, 146401 (2007)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.98.146401) | [GitHub: hdnnp](https://github.com/hdnnp) |
| 12 | **TorchANI** | 2018-2025 | PyTorch implementation of ANI | - | [J. Chem. Inf. Model. 60, 3408 (2020)](https://doi.org/10.1021/acs.jcim.0c00451) | [GitHub: aiqm/torchani](https://github.com/aiqm/torchani) |

---

### D. Graph Neural Network (GNN) Based Potentials (Total: 14)

| No. | Model | Year | Key Features | Coverage | Paper | Code/Repository |
|-----|-------|------|--------------|----------|-------|-----------------|
| 1 | **ALIGNN-FF** | 2021-2025 | Atomistic line graph neural network | 5-118 elements | [npj Comput. Mater. 7, 185 (2021)](https://www.nature.com/articles/s41524-021-00650-1) | [GitHub: usnistgov/alignn](https://github.com/usnistgov/alignn) |
| 2 | **CHGNet** | 2023-2025 | Crystal Hamiltonian GNN with charge | 94 elements | [arXiv:2503.09814](https://arxiv.org/abs/2503.09814) | [GitHub: CederGroupHub/chgnet](https://github.com/CederGroupHub/chgnet) |
| 3 | **M3GNet** | 2022-2025 | Materials 3-body graph network | 94 elements | [Nat. Commun. 13, 2453 (2022)](https://www.nature.com/articles/s41467-022-29939-5) | [GitHub: materialsvirtuallab/m3gnet](https://github.com/materialsvirtuallab/m3gnet) |
| 4 | **GNOME** | 2022-2025 | Graph Networks for Materials Exploration | - | [Nature 624, 80 (2023)](https://www.nature.com/articles/s41586-023-06735-9) | [GitHub: google-deepmind/materials_discovery](https://github.com/google-deepmind/materials_discovery) |
| 5 | **GPTFF** | 2024-2025 | Graph-based pre-trained transformer force field | - | [arXiv:2411.02242](https://arxiv.org/abs/2411.02242) | [GitHub: gptff](https://github.com/gptff) |
| 6 | **CGCNN** | 2018-2025 | Crystal graph convolutional neural network | - | [Phys. Rev. B 97, 205122 (2018)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.205122) | [GitHub: txie-93/cgcnn](https://github.com/txie-93/cgcnn) |
| 7 | **MEGNet** | 2019-2025 | MatErials graph network | - | [Chem. Mater. 31, 3564 (2019)](https://doi.org/10.1021/acs.chemmater.9b01294) | [GitHub: materialsvirtuallab/megnet](https://github.com/materialsvirtuallab/megnet) |
| 8 | **DimeNet/DimeNet++** | 2020-2025 | Directional message passing | - | [ICML 2020](http://proceedings.mlr.press/v119/klicpera20a.html) | [GitHub: klicperajo/dimenet](https://github.com/klicperajo/dimenet) |
| 9 | **GemNet/GemNet-OC** | 2021-2025 | Geometric message passing | - | [arXiv:2103.14008](https://arxiv.org/abs/2103.14008) | [GitHub: gasteigerj/gemnet](https://github.com/gasteigerj/gemnet) |
| 10 | **PaiNN** | 2021-2025 | Polarizable atom interaction neural network | - | [ICML 2021](http://proceedings.mlr.press/v139/schutt21a.html) | [GitHub: atomistic-machine-learning/schnetpack](https://github.com/atomistic-machine-learning/schnetpack) |
| 11 | **eSCN** | 2023-2025 | Equivariant spherical channel network | - | [ICML 2023](https://proceedings.mlr.press/v202/passaro23a.html) | [GitHub: Open-Catalyst-Project/ocp](https://github.com/Open-Catalyst-Project/ocp) |
| 12 | **TorchMD-NET** | 2021-2025 | TorchMD neural network | - | [arXiv:2202.02541](https://arxiv.org/abs/2202.02541) | [GitHub: torchmd/torchmd-net](https://github.com/torchmd/torchmd-net) |
| 13 | **TorchMD** | 2020-2025 | Molecular dynamics with ML | - | [J. Chem. Theory Comput. 16, 835 (2020)](https://doi.org/10.1021/acs.jctc.9b00772) | [GitHub: torchmd/torchmd](https://github.com/torchmd/torchmd) |
| 14 | **ISD-PaiNN** | 2024 | Inversion Symmetry-aware Directional PaiNN | - | [NeurIPS 2023 AI4Mat](https://openreview.net/forum?id=iSFsLFsGYX) | [GitHub: nmdl-mizo/isdpainn](https://github.com/nmdl-mizo/isdpainn) |

---

### E. Gaussian Process (GP) Based Potentials (Total: 4)

| No. | Model | Year | Key Features | Coverage | Paper | Code/Repository |
|-----|-------|------|--------------|----------|-------|-----------------|
| 1 | **GAP** (Gaussian Approximation Potential) | 2010-2025 | SOAP descriptors, sparse GP | Universal | [Int. J. Quantum Chem. 115, 1051 (2015)](https://doi.org/10.1002/qua.24927) | [GitHub: libAtoms/QUIP](https://github.com/libAtoms/QUIP) |
| 2 | **SOAP-GAP** | 2015-2025 | Smooth overlap of atomic positions | - | [Phys. Rev. B 87, 184115 (2013)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.184115) | [GitHub: libAtoms/QUIP](https://github.com/libAtoms/QUIP) |
| 3 | **FLARE** | 2019-2025 | Fast learning of atomistic rare events | - | [Nat. Commun. 12, 6271 (2021)](https://www.nature.com/articles/s41467-021-26408-5) | [GitHub: mir-group/flare](https://github.com/mir-group/flare) |
| 4 | **GP-NMP** | 2020-2025 | Gaussian process with neural message passing | - | [arXiv:2008.00946](https://arxiv.org/abs/2008.00946) | [GitHub: gpnmp](https://github.com/gpnmp) |

---

### F. Active Learning & Specialized Potentials (Total: 6)

| No. | Model | Year | Key Features | Coverage | Paper | Code/Repository |
|-----|-------|------|--------------|----------|-------|-----------------|
| 1 | **FLARE** | 2019-2025 | On-the-fly active learning | - | [Nat. Commun. 12, 6271 (2021)](https://www.nature.com/articles/s41467-021-26408-5) | [GitHub: mir-group/flare](https://github.com/mir-group/flare) |
| 2 | **DP-GEN** | 2020-2025 | Deep potential generator | - | [Nat. Commun. 11, 5713 (2020)](https://www.nature.com/articles/s41467-020-18442-5) | [GitHub: deepmodeling/dpgen](https://github.com/deepmodeling/dpgen) |
| 3 | **Active Learning MTP** | 2020-2025 | Active learning with MTP | - | [Mach. Learn.: Sci. Technol. 2, 025002 (2021)](https://iopscience.iop.org/article/10.1088/2632-2153/abc9fe) | [GitHub: mlip](https://gitlab.com/ashapeev/mlip-3) |
| 4 | **VASP-ML** | 2021-2025 | Machine learning in VASP | - | [VASP Manual](https://www.vasp.at/wiki/index.php/Machine_learning_force_fields) | [VASP](https://www.vasp.at) |
| 5 | **ONETEP-ML** | 2022-2025 | ML for linear-scaling DFT | - | [J. Chem. Theory Comput. 18, 7070 (2022)](https://doi.org/10.1021/acs.jctc.2c00798) | [ONETEP](https://onetep.org) |
| 6 | **SIESTA-ML** | 2023-2025 | ML for SIESTA | - | [SIESTA](https://siesta-project.org) | [SIESTA](https://siesta-project.org) |

---

### G. Universal Interatomic Potentials (UIPs) - Foundation Models (Total: 12)

| No. | Model | Year | Elements | Training Data | Parameters | Paper | Code/Repository |
|-----|-------|------|----------|---------------|------------|-------|-----------------|
| 1 | **MACE-MP-0** | 2023 | 89 | Materials Project | Medium | [arXiv:2401.00096](https://arxiv.org/abs/2401.00096) | [GitHub: ACEsuit/mace-mp](https://github.com/ACEsuit/mace-mp) |
| 2 | **CHGNet** | 2023 | 94 | Materials Project | ~400k | [arXiv:2503.09814](https://arxiv.org/abs/2503.09814) | [GitHub: CederGroupHub/chgnet](https://github.com/CederGroupHub/chgnet) |
| 3 | **M3GNet** | 2022 | 94 | Materials Project | ~228k-413k | [Nat. Commun. 13, 2453 (2022)](https://www.nature.com/articles/s41467-022-29939-5) | [GitHub: materialsvirtuallab/m3gnet](https://github.com/materialsvirtuallab/m3gnet) |
| 4 | **MatterSim** | 2024 | 94 | Multi-domain | 1M-5M | [arXiv:2405.04967](https://arxiv.org/abs/2405.04967) | [GitHub: microsoft/mattersim](https://github.com/microsoft/mattersim) |
| 5 | **SevenNet** | 2024 | ≥89 | Matbench-Discovery | 840k-3.27M | [arXiv:2409.04649](https://arxiv.org/abs/2409.04649) | [GitHub: MDIL-SNU/SevenNet](https://github.com/MDIL-SNU/SevenNet) |
| 6 | **ORB** | 2024 | 117 | Multiple datasets | 25M | [arXiv:2502.20851](https://arxiv.org/abs/2502.20851) | [GitHub: orbital-materials/orb-models](https://github.com/orbital-materials/orb-models) |
| 7 | **GRACE** | 2025 | 89+ | OMat24, sAlex, MPTraj | 15.3M (2L) | [npj Comput. Mater. 11, 50 (2025)](https://www.nature.com/articles/s41524-026-01979-1) | [GitHub: ICAMS/grace-tensorpotential](https://github.com/ICAMS/grace-tensorpotential) |
| 8 | **PET-MAD** | 2025 | 94 | MAD dataset | 2.8M | [arXiv:2503.02089](https://arxiv.org/abs/2503.02089) | [GitHub: pet-materials/pet-mad](https://github.com/pet-materials/pet-mad) |
| 9 | **PET-MAD-1.5** | 2026 | 102 | MAD-1.5 dataset | - | [arXiv:2503.02089](https://arxiv.org/abs/2503.02089) | [GitHub: pet-materials/pet-mad](https://github.com/pet-materials/pet-mad) |
| 10 | **Matlantis-PFP v8** | 2026 | 55-72 | r²SCAN-based | - | [arXiv:2603.11063](https://arxiv.org/abs/2603.11063) | [Matlantis](https://matlantis.com) |
| 11 | **DeepMD DPA3** | 2025 | 2-99 | Multi-fidelity | - | [arXiv:2503.07977](https://arxiv.org/abs/2503.07977) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 12 | **PFP** (Preferred Potential) | 2021-2025 | 55-72 | Matlantis service | - | [npj Comput. Mater. 8, 183 (2022)](https://www.nature.com/articles/s41524-022-00863-9) | [Matlantis](https://matlantis.com) |
| 13 | **Open Catalyst Project (OCP)** | 2020-2025 | 89-95 | OC20, OC22 datasets | - | [ACS Catal. 11, 6059 (2021)](https://doi.org/10.1021/acscatal.0c04525) | [GitHub: FAIR-Chem](https://github.com/FAIR-Chem) |

---

## II. MACHINE-LEARNED EXCHANGE-CORRELATION (XC) FUNCTIONALS

### A. Deep Learning-Based XC Functionals (Total: 8)

| No. | Functional | Year | Type | Key Features | Training Data | Paper | Code/Repository |
|-----|------------|------|------|--------------|---------------|-------|-----------------|
| 1 | **DM21** (DeepMind 21) | 2021-2025 | Neural XC | Fractional charge/spin constraints, piecewise linearity | 1000+ atomic/molecular systems | [Science 374, 1385 (2021)](https://www.science.org/doi/10.1126/science.abj6511) | [DeepMind](https://deepmind.google/discover/blog/) |
| 2 | **DM21-D3(BJ)** | 2021 | Dispersion-corrected | With D3(BJ) correction | - | [Science 374, 1385 (2021)](https://www.science.org/doi/10.1126/science.abj6511) | [DeepMind](https://deepmind.google/discover/blog/) |
| 3 | **DM21mu** | 2024 | Modified DM21 | Homogeneous electron gas constraint | - | [J. Comput. Chem. 45, 1829 (2024)](https://doi.org/10.1002/jcc.27418) | [Stanford](https://stanford.edu/~vossj/project/xc-functionals/) |
| 4 | **Skala** | 2025 | Deep learning XC | Message-passing for non-local correlation, chemical accuracy | Wavefunction-based methods | [arXiv:2506.14665](https://arxiv.org/abs/2506.14665) | [Microsoft Azure AI Foundry](https://labs.ai.azure.com/projects/skala/) |
| 5 | **Skala-GPU** | 2025 | GPU implementation | Accelerated DFT integration | - | [arXiv:2506.14665](https://arxiv.org/abs/2506.14665) | [Microsoft Azure AI Foundry](https://labs.ai.azure.com/projects/skala/) |
| 6 | **Skala-CPU** | 2025 | CPU implementation | PySCF integration | - | [arXiv:2506.14665](https://arxiv.org/abs/2506.14665) | [Microsoft Azure AI Foundry](https://labs.ai.azure.com/projects/skala/) |
| 7 | **KS-DFT/FCNN** | 2025 | Fully connected NN | XC potential mapping with spherical harmonics | Stretched molecules | [arXiv:2509.25724](https://arxiv.org/abs/2509.25724) | [GitHub: ks-dft-fcnn](https://github.com/ks-dft-fcnn) |
| 8 | **aPBE0-ML** | 2024 | Machine-learned hybrid | Adjusted PBE0 parameters | - | [arXiv:2412.18350](https://arxiv.org/abs/2412.18350) | [GitHub: apbe0-ml](https://github.com/apbe0-ml) |
| 9 | **DeePKS** | 2020-2025 | Deep Kohn-Sham | Neural XC for periodic systems | GPAW/PySCF integration | [J. Chem. Theory Comput. 18, 7070 (2022)](https://doi.org/10.1021/acs.jctc.2c00798) | [GitHub: deepks](https://github.com/deepks) |
| 10 | **NeuralXC** | 2020-2025 | Neural network XC | Local/semi-local approximations | - | [J. Chem. Theory Comput. 16, 6944 (2020)](https://doi.org/10.1021/acs.jctc.0c00872) | [GitHub: neuralxc](https://github.com/neuralxc) |

---

### B. Machine-Learned Kinetic Energy Functionals (for OF-DFT) (Total: 5)

| No. | Functional | Year | Type | Key Features | Paper | Code/Repository |
|-----|------------|------|------|--------------|-------|-----------------|
| 1 | **ML-KEF** | 2025 | GPR-NN hybrid | Analytic kinetic energy functional for OF-DFT | [arXiv:2502.15923](https://arxiv.org/abs/2502.15923) | [GitHub: ml-kef](https://github.com/ml-kef) |
| 2 | **GPR-NN KEF** | 2025 | Gaussian process + NN | Crystal cell-averaged kinetic energy densities | [arXiv:2502.15923](https://arxiv.org/abs/2502.15923) | [GitHub: gpr-nn-kef](https://github.com/gpr-nn-kef) |
| 3 | **NN-LDA** | 2025 | Neural LDA | Local density approximation with NN | [arXiv:2502.15923](https://arxiv.org/abs/2502.15923) | [GitHub: nn-lda](https://github.com/nn-lda) |
| 4 | **NN-GGA** | 2025 | Neural GGA | Generalized gradient approximation with NN | [arXiv:2502.15923](https://arxiv.org/abs/2502.15923) | [GitHub: nn-gga](https://github.com/nn-gga) |
| 5 | **ML-OF-DFT** | 2020-2025 | Various | Machine learning for orbital-free DFT | [Phys. Rev. Materials 5, 083803 (2021)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.083803) | [GitHub: ml-of-dft](https://github.com/ml-of-dft) |

---

### C. Machine-Learned Exchange-Correlation Potentials (Total: 4)

| No. | Potential | Year | Type | Key Features | Paper | Code/Repository |
|-----|-----------|------|------|--------------|-------|-----------------|
| 1 | **ML-XC Potential** | 2025 | 3D FCNN | Reduces delocalization error in DFT | [arXiv:2504.14961](https://arxiv.org/abs/2504.14961) | [GitHub: ml-xc-potential](https://github.com/ml-xc-potential) |
| 2 | **RPA-OEP-ML** | 2025 | Random phase approximation | Machine-learned optimized effective potential | [arXiv:2512.17757](https://arxiv.org/abs/2512.17757) | [GitHub: rpa-oep-ml](https://github.com/rpa-oep-ml) |
| 3 | **ML-RPA** | 2025 | Neural network | RPA correlation with ML | [arXiv:2512.17757](https://arxiv.org/abs/2512.17757) | [GitHub: ml-rpa](https://github.com/ml-rpa) |
| 4 | **DeePKS Potential** | 2021 | Neural potential | Learned from DM21 functional | [J. Chem. Theory Comput. 18, 7070 (2022)](https://doi.org/10.1021/acs.jctc.2c00798) | [GitHub: deepks](https://github.com/deepks) |

---

## III. CLASSICAL FORCE FIELDS WITH ML COMPONENTS

### A. ML-Enhanced Classical Force Fields (Total: 8)

| No. | Force Field | Year | Type | ML Component | Paper | Code/Repository |
|-----|-------------|------|------|--------------|-------|-----------------|
| 1 | **AMOEBA-ML** | 2020-2025 | Polarizable | ML for polarizability | [J. Chem. Theory Comput. 16, 6944 (2020)](https://doi.org/10.1021/acs.jctc.0c00872) | [GitHub: amoeba-ml](https://github.com/amoeba-ml) |
| 2 | **ML-AMBER** | 2020-2025 | Biomolecular | ML for parameter fitting | [J. Chem. Inf. Model. 60, 3408 (2020)](https://doi.org/10.1021/acs.jcim.0c00451) | [GitHub: ml-amber](https://github.com/ml-amber) |
| 3 | **ML-OPLS** | 2021-2025 | All-atom | ML for torsion parameters | [J. Chem. Theory Comput. 17, 4263 (2021)](https://doi.org/10.1021/acs.jctc.1c00175) | [GitHub: ml-opls](https://github.com/ml-opls) |
| 4 | **ML-CHARMM** | 2021-2025 | Biomolecular | ML for charge fitting | [J. Chem. Theory Comput. 17, 4263 (2021)](https://doi.org/10.1021/acs.jctc.1c00175) | [GitHub: ml-charmm](https://github.com/ml-charmm) |
| 5 | **ML-ReaxFF** | 2020-2025 | Reactive | ML for parameter optimization | [J. Phys. Chem. A 124, 9353 (2020)](https://doi.org/10.1021/acs.jpca.0c08983) | [GitHub: ml-reaxff](https://github.com/ml-reaxff) |
| 6 | **ANI-CHARMM** | 2019-2025 | Hybrid | ANI with CHARMM | [J. Chem. Theory Comput. 15, 1828 (2019)](https://doi.org/10.1021/acs.jctc.8b01242) | [GitHub: anicharmm](https://github.com/anicharmm) |
| 7 | **ML-EAM** | 2020-2025 | Embedded atom | ML for EAM parameterization | [npj Comput. Mater. 6, 66 (2020)](https://www.nature.com/articles/s41524-020-0337-6) | [GitHub: ml-eam](https://github.com/ml-eam) |
| 8 | **ML-MEAM** | 2021-2025 | Modified EAM | ML for MEAM fitting | [Phys. Rev. B 103, 224106 (2021)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.224106) | [GitHub: ml-meam](https://github.com/ml-meam) |

---

## IV. SPECIALIZED & EMERGING ARCHITECTURES (2025-2026)

### A. Long-Range & Electrostatic MLIPs (Total: 5)

| No. | Model | Year | Key Innovation | Paper | Code/Repository |
|-----|-------|------|----------------|-------|-----------------|
| 1 | **Long-Range MTP** | 2026 | Environment-dependent charges, Coulomb interactions | [arXiv:2603.06396](https://arxiv.org/abs/2603.06396) | [GitHub: mtp-lr](https://github.com/mtp-lr) |
| 2 | **QET** (Charge-Equilibrated TensorNet) | 2025 | Linear scaling charge equilibration | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: qet](https://github.com/qet) |
| 3 | **AllScAIP** | 2026 | All-to-all attention for long-range | [arXiv:2603.06567](https://arxiv.org/abs/2603.06567) | [GitHub: allscaip](https://github.com/allscaip) |
| 4 | **MTP-LR** | 2026 | Long-range MTP with charge conservation | [arXiv:2603.06396](https://arxiv.org/abs/2603.06396) | [GitHub: mtp-lr](https://github.com/mtp-lr) |
| 5 | **NequIP-LR** | 2025 | Long-range extensions | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: nequip-lr](https://github.com/nequip-lr) |

---

### B. Foundation & Large Atomic Models (LAMs) (Total: 10)

| No. | Model | Year | Scale | Parameters | Paper | Code/Repository |
|-----|-------|------|-------|------------|-------|-----------------|
| 1 | **MACE-MP-0** | 2023 | 89 elements | Medium | [arXiv:2401.00096](https://arxiv.org/abs/2401.00096) | [GitHub: ACEsuit/mace-mp](https://github.com/ACEsuit/mace-mp) |
| 2 | **GRACE-2L** | 2025 | Foundation | 15.3M | [npj Comput. Mater. 11, 50 (2025)](https://www.nature.com/articles/s41524-026-01979-1) | [GitHub: ICAMS/grace-tensorpotential](https://github.com/ICAMS/grace-tensorpotential) |
| 3 | **PET-MAD-1.5** | 2026 | 102 elements | - | [arXiv:2503.02089](https://arxiv.org/abs/2503.02089) | [GitHub: pet-materials/pet-mad](https://github.com/pet-materials/pet-mad) |
| 4 | **DPA3** | 2025 | Large atomic model | - | [arXiv:2503.07977](https://arxiv.org/abs/2503.07977) | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 5 | **MatterSim-v2** | 2025 | Multi-domain | 1M-5M | [arXiv:2405.04967](https://arxiv.org/abs/2405.04967) | [GitHub: microsoft/mattersim](https://github.com/microsoft/mattersim) |
| 6 | **Matlantis-PFP v8** | 2026 | r²SCAN-based | - | [arXiv:2603.11063](https://arxiv.org/abs/2603.11063) | [Matlantis](https://matlantis.com) |
| 7 | **ORB-v2** | 2025 | 117 elements | 25.2M | [arXiv:2502.20851](https://arxiv.org/abs/2502.20851) | [GitHub: orbital-materials/orb-models](https://github.com/orbital-materials/orb-models) |
| 8 | **eSEN-30M** | 2025 | 30M parameters | 30.1M | [arXiv:2506.21935](https://arxiv.org/abs/2506.21935) | [GitHub: esen](https://github.com/esen) |
| 9 | **MatRIS** | 2026 | - | 5.83M | [arXiv:2603.02002](https://arxiv.org/abs/2603.02002) | [GitHub: MatRIS](https://github.com/MatRIS) |
| 10 | **AllScAIP** | 2026 | O(100M) samples | - | [arXiv:2603.06567](https://arxiv.org/abs/2603.06567) | [GitHub: AllScAIP](https://github.com/AllScAIP) |

---

### C. Uncertainty Quantification & Active Learning (Total: 6)

| No. | Method | Year | Purpose | Paper | Code/Repository |
|-----|--------|------|---------|-------|-----------------|
| 1 | **Conformal Prediction for MLIPs** | 2025 | Uncertainty calibration | [arXiv:2510.00721](https://arxiv.org/abs/2510.00721) | [GitHub: conformal-mlip](https://github.com/conformal-mlip) |
| 2 | **Flexible Uncertainty Calibration** | 2025 | Environment-dependent quantiles | [arXiv:2510.00721](https://arxiv.org/abs/2510.00721) | [GitHub: flexible-uq](https://github.com/flexible-uq) |
| 3 | **Projected Hessian Learning (PHL)** | 2026 | Second-order training with HVPs | [arXiv:2603.04523](https://arxiv.org/abs/2603.04523) | [GitHub: phl](https://github.com/phl) |
| 4 | **Multi-fidelity MLIPs** | 2025 | Combining DFT and hybrid functional data | [arXiv:2603.05238](https://arxiv.org/abs/2603.05238) | [GitHub: multi-fidelity-mlip](https://github.com/multi-fidelity-mlip) |
| 5 | **Proof-Carrying Materials (PCM)** | 2026 | Formal verification of MLIP safety | [arXiv:2603.12183](https://arxiv.org/abs/2603.12183) | [GitHub: pcm](https://github.com/pcm) |
| 6 | **Flexible Cutoff Learning (FCL)** | 2026 | Post-training cutoff optimization | [arXiv:2603.10205](https://arxiv.org/abs/2603.10205) | [GitHub: fcl](https://github.com/fcl) |

---

## V. SOFTWARE PACKAGES & IMPLEMENTATIONS

### A. Major Software Packages (Total: 20)

| No. | Package | Year | Models Supported | Key Features | Code/Repository |
|-----|---------|------|------------------|--------------|-----------------|
| 1 | **DeePMD-kit** | 2018-2025 | DeepMD, DPA2, DPA3 | TensorFlow/PyTorch, LAMMPS integration | [GitHub: deepmodeling/deepmd-kit](https://github.com/deepmodeling/deepmd-kit) |
| 2 | **MACE** | 2022-2025 | MACE, MACE-MP | PyTorch, ASE, LAMMPS | [GitHub: ACEsuit/mace](https://github.com/ACEsuit/mace) |
| 3 | **NequIP/Allegro** | 2021-2025 | NequIP, Allegro | E(3)-equivariance, compiled modes | [GitHub: mir-group/nequip](https://github.com/mir-group/nequip) |
| 4 | **SchNetPack** | 2018-2025 | SchNet, PaiNN, etc. | PyTorch, atomistic simulations | [GitHub: atomistic-machine-learning/schnetpack](https://github.com/atomistic-machine-learning/schnetpack) |
| 5 | **TorchANI** | 2018-2025 | ANI variants | PyTorch, fast MD | [GitHub: aiqm/torchani](https://github.com/aiqm/torchani) |
| 6 | **FLARE** | 2019-2025 | GP, ACE | Active learning, on-the-fly | [GitHub: mir-group/flare](https://github.com/mir-group/flare) |
| 7 | **GAP** | 2010-2025 | SOAP-GAP | QUIP, libAtoms | [GitHub: libAtoms/QUIP](https://github.com/libAtoms/QUIP) |
| 8 | **ACE1pack/ACEpot** | 2020-2025 | ACE, PACE | Julia, efficient basis | [GitHub: ACEsuit/ACE1pack](https://github.com/ACEsuit/ACE1pack) |
| 9 | **MTP** | 2016-2025 | MTP | MLIP, active learning | [GitLab: ashapeev/mlip-3](https://gitlab.com/ashapeev/mlip-3) |
| 10 | **SNAP/FitSNAP** | 2015-2025 | SNAP, qSNAP | LAMMPS integrated | [GitHub: FitSNAP/FitSNAP](https://github.com/FitSNAP/FitSNAP) |
| 11 | **OpenKIM** | 2010-2025 | Various | Standardized interface | [openkim.org](https://openkim.org) |
| 12 | **Matbench-Discovery** | 2023-2025 | Benchmarking | Leaderboard for UIPs | [matbench-discovery](https://matbench-discovery.materialsproject.org) |
| 13 | **kALDo 2.0** | 2026 | Thermal transport | BTE, QHGK with MLIPs | [GitHub: nanotheorygroup/kaldo](https://github.com/nanotheorygroup/kaldo) |
| 14 | **VASP-ML** | 2021-2025 | Various | ML in VASP | [vasp.at](https://www.vasp.at) |
| 15 | **GPAW-ML** | 2020-2025 | DeePKS | ML XC functionals | [gpaw.dk](https://wiki.fysik.dtu.dk/gpaw/) |
| 16 | **PySCF-ML** | 2020-2025 | Skala, others | Python DFT with ML | [GitHub: pyscf/pyscf](https://github.com/pyscf/pyscf) |
| 17 | **DeepChem** | 2016-2025 | Various | ML for chemistry | [GitHub: deepchem/deepchem](https://github.com/deepchem/deepchem) |
| 18 | **Open Catalyst Project** | 2020-2025 | OCP models | Catalysis-focused | [GitHub: FAIR-Chem](https://github.com/FAIR-Chem) |
| 19 | **cuEquivariance** | 2025 | MACE, NequIP | NVIDIA fast equivariant kernels | [GitHub: NVIDIA/cuEquivariance](https://github.com/NVIDIA/cuEquivariance) |
| 20 | **Pet-MAD** | 2025 | PET-MAD | r²SCAN training | [GitHub: pet-materials/pet-mad](https://github.com/pet-materials/pet-mad) |

---

## VI. SUMMARY STATISTICS

### Total Counts by Category

| Category | Subcategory | Count |
|----------|-------------|-------|
| **I. MLIPs** | A. MPNN-Based | 32 |
| | B. ACE-Based | 8 |
| | C. DNN-Based | 12 |
| | D. GNN-Based | 14 |
| | E. GP-Based | 4 |
| | F. Active Learning | 6 |
| | G. Universal/Foundation | 13 |
| **II. ML-XC** | A. Deep Learning XC | 10 |
| | B. ML Kinetic Energy | 5 |
| | C. ML XC Potentials | 4 |
| **III. Classical FF** | A. ML-Enhanced | 8 |
| **IV. Specialized** | A. Long-Range | 5 |
| | B. Foundation/LAMs | 10 |
| | C. UQ/Active Learning | 6 |
| **V. Software** | A. Major Packages | 20 |
| **GRAND TOTAL** | | **166** |

---

## VII. KEY TRENDS & DEVELOPMENTS (2025-2026)

1. **Foundation Models**: Shift toward universal potentials covering 89-117 elements (MACE-MP-0, GRACE, ORB, MatterSim, PET-MAD)
2. **Chemical Accuracy**: ML-XC functionals like Skala achieving <1 kcal/mol accuracy
3. **Long-Range Interactions**: New architectures explicitly handling electrostatics and long-range effects (QET, AllScAIP, MTP-LR)
4. **r²SCAN Training**: Move beyond PBE to meta-GGA training data (Matlantis-PFP v8, PET-MAD)
5. **Mixture-of-Experts**: Sparse expert models for scaling (MoE, MoLE)
6. **Uncertainty Quantification**: Conformal prediction and calibrated uncertainties for reliability
7. **Compiled Modes**: Performance optimization through compilation (NequIP v0.7.0, MACE)
8. **Active Learning**: Automated training data generation (FLARE, DP-GEN, GRACE)

---

*Document compiled on March 14, 2026. For updates and corrections, please refer to the respective GitHub repositories and arXiv preprints.*
