# GNoME

## Overview
GNoME (Graph Networks for Materials Exploration) is Google DeepMind's deep learning tool for predicting the stability of inorganic crystal structures. It discovered 2.2 million new stable crystals, including 380,000 added to the Materials Project.

## Theoretical Basis
- Graph neural networks (GNN)
- Active learning
- Formation energy prediction
- Stability classification
- Convex hull analysis

## Key Capabilities
- Crystal stability prediction
- Large-scale materials discovery
- Formation energy prediction
- Active learning workflow
- Materials Project integration

**Sources**: Nature 2023, Google DeepMind

## Key Strengths

### Scale:
- 2.2 million new crystals
- 380,000 stable materials
- Periodic table coverage

### Methodology:
- Graph neural networks
- Active learning
- DFT validation

### Impact:
- Materials Project integration
- Public dataset
- Landmark discovery

## Inputs & Outputs
- **Input formats**: Crystal structures
- **Output data types**: Stability predictions, formation energies

## Interfaces & Ecosystem
- **Materials Project**: Data integration
- **Datasets**: Public GNoME database
- **Colab**: Example notebooks

## Workflow and Usage
1. Input crystal structure
2. GNN predicts formation energy
3. Assess stability vs convex hull
4. Identify stable candidates
5. Validate with DFT

## Performance Characteristics
- Fast inference
- High accuracy
- Large-scale applicable

## Computational Cost
- Inference: fast
- Training: significant (DeepMind scale)
- Validation: DFT

## Best Practices
- Use for stability screening
- Validate predictions with DFT
- Check against convex hull
- Consider synthesizability

## Limitations & Known Constraints
- Prediction only (not generative)
- Training data dependent
- May miss metastable phases

## Application Areas
- Materials stability screening
- High-throughput discovery
- Database expansion
- Materials exploration

## Comparison with Other Codes
- **vs MatterGen**: GNoME predictive, MatterGen generative
- **vs CDVAE**: Different purpose (prediction vs generation)
- **Unique strength**: Massive scale, DeepMind backing, Materials Project integration

## Community and Support
- Google DeepMind
- Public dataset
- Colab examples
- Nature publication

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/google-deepmind/materials_discovery
2. Publication: Nature (2023)
3. Blog: https://deepmind.google/discover/blog/millions-of-new-materials-discovered-with-deep-learning/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Source: OPEN (dataset)
- Development: DeepMind
- Applications: Materials stability prediction, discovery
