# OrbNet

## Official Resources
- Homepage: https://entos.ai/
- Documentation: Publications / Entos documentation
- Source Repository: Proprietary
- License: Commercial

## Overview
OrbNet is an AI-driven quantum chemistry method and software developed by Entos, Inc. It utilizes graph neural networks (Geometric Deep Learning) and domain-specific features (based on low-cost quantum calculations like GFN-xTB or semi-empirical methods) to predict high-level quantum chemical properties (like DFT or CCSD(T) energies) with high accuracy and vastly reduced computational cost.

**Scientific domain**: AI for Science, Machine Learning Potentials, Quantum Chemistry  
**Target user community**: Pharma, Materials Science, High-throughput screening

## Theoretical Methods
- Geometric Deep Learning (Graph Neural Networks)
- Symmetry-preserving representations
- Delta-learning (Correction to low-level method)
- Semi-empirical baselines (e.g., xTB)
- Density Functional Theory (Target usage)
- Coupled Cluster (Target usage)

## Capabilities (CRITICAL)
- Prediction of molecular energies and forces
- Geometry optimization using AI potentials
- 1000x speedup over DFT
- Accuracy comparable to DFT/CCSD(T) (within domain)
- Scalable to large molecules (proteins, supramolecular)
- GPU acceleration for inference

## Key Strengths

### Speed vs Accuracy:
- DFT accuracy at semi-empirical cost
- 3-4 orders of magnitude faster than DFT
- Enables dynamics on QC surfaces
- High-throughput compatible

### Physics-Aware ML:
- Uses quantum features (orbitals/density approx)
- Better generalization than pure geometry ML
- "Grey-box" approach combining physics with data
- Size extensive properties by design
- Rotational and translational invariance guaranteed
- Captures long-range electronic effects via attention mechanisms

### Scalability:
- Linear scaling inference
- Handles systems with thousands of atoms
- Protein-ligand binding capability
- Supramolecular systems

## Inputs & Outputs
- **Input**: Molecular Structure (XYZ/PDB)
- **Output**: Energy, Forces, Properties
- **Interface**: Python API, Entos platform

## Interfaces & Ecosystem
- **Entos Platform**: Integrated with Qcore
- **Python**: Pytorch-based backend likely
- **Cloud**: Often deployed as cloud service

## Advanced Features

### Graph Neural Networks:
- Node/Edge embeddings
- Message passing architecture
- Rotationally invariant
- Transferable features

### Training Data:
- Trained on large QC datasets
- Active learning support
- Domain adaptation

## Performance Characteristics
- **Speed**: Milliseconds to seconds per evaluation
- **Accuracy**: ~1-2 kcal/mol relative to high-level reference
- **System size**: Up to thousands of atoms
- **Hardware**: GPU accelerated inference

## Computational Cost
- **Inference**: Negligible compared to QC
- **Pre-calculation**: Cost of semi-empirical input (cheap)
- **Throughput**: Extremely high

## Limitations & Known Constraints
- **Domain applicability**: Valid within training production
- **Black/Grey box**: ML limitations on outliers
- **Proprietary**: Commercial access
- **Source**: Not open

## Comparison with Other Codes
- **vs ANI/SchNet**: OrbNet uses electronic features (more robust)
- **vs DFT**: OrbNet is approximation but 1000x faster
- **vs Force Fields**: OrbNet is reactive and electronic-aware
- **Unique strength**: Electronic-structure-aware Geometric Deep Learning

## Application Areas
### Pharmaceutical Research:
- **Lead Optimization**: rapid free energy perturbation (FEP) cycles
- **Docking Scoring**: Physics-based scoring function replacement
- **Library Enumeration**: Property prediction for millions of candidates
- **Toxicity Prediction**: Electronic descriptors for ADMET properties

### Material Design:
- **Polymer Properties**: Glass transition and mechanical property prediction
- **Screening**: High-throughput evaluation of organic electronics
- **Crystal Structure**: Ranking of polymorphs

### Biomolecular Simulation:
- **Reactions in Enzymes**: Modeling QM/MM regions with full QM speed
- **Protein Dynamics**: Long-timescale dynamics with electronic structure accuracy
- **Allosteric Effects**: Capturing subtle electronic changes in large systems

## Best Practices
### Model Validation:
- **Outlier Detection**: Flag structures with high uncertainty
- **Reference Checks**: Periodically verify against full DFT/CCSD(T)
- **Chemical Space**: Ensure target molecules fall within the training domain

### Hardware Optimization:
- **Batch Size**: Maximize GPU utilization with large batch sizes
- **Precision**: Use mixed precision (TF32/FP16) where supported for speed
- **Multi-GPU**: Distribute inference across available accelerators

### Workflow Integration:
- **Retraining**: Fine-tune models on project-specific data if available
- **Ensembling**: Use model ensembles to estimate prediction error
- **Hybrid Methods**: Use OrbNet for exploration, high-level QM for final confirmation

## Community and Support
- Entos, Inc.
- Commercial support

## Verification & Sources
**Primary sources**:
1. Entos.ai
2. "OrbNet: Deep Learning for Quantum Chemistry using Symmetry-Adapted Atomic-Orbital Features", J. Chem. Phys. (2020)
3. Anandkumar et al. publications

**Confidence**: VERIFIED
- Status: Active commercial technology
- Methodology: Published in high-impact journals
