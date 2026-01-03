# m3gnet

## Official Resources
- Homepage: https://materialsvirtuallab.github.io/m3gnet/
- Documentation: https://materialsvirtuallab.github.io/m3gnet/
- Source Repository: https://github.com/materialsvirtuallab/m3gnet
- License: BSD 3-Clause License

## Overview
M3GNet is a graph neural network (GNN) potential trained on the massive Materials Project trajectory dataset (>187,000 materials). It serves as a "universal" interatomic potential for the periodic table, capable of relaxing structures and performing MD for diverse chemistries. It can also be used as a surrogate model for property prediction.

**Scientific domain**: Universal ML potentials, graph neural networks  
**Target user community**: Materials scientists needing quick relaxation/properties

## Capabilities (CRITICAL)
- **Universal Potential**: Covers 89 elements.
- **Relaxation**: Can relax arbitrary crystal structures (often replacing DFT relaxation).
- **MD**: Run molecular dynamics.
- **Surrogate**: Predicts formation energy, band gap, elastic moduli directly from structure.
- **Integration**: Works with Pymatgen and ASE.

**Sources**: M3GNet GitHub, Nat. Comput. Sci. 2, 718 (2022)

## Inputs & Outputs
- **Input formats**: Pymatgen Structure / ASE Atoms
- **Output data types**: Energy, forces, stress, properties

## Interfaces & Ecosystem
- **TensorFlow**: Backend.
- **Pymatgen**: Native object support.
- **ASE**: Calculator interface.

## Workflow and Usage
1. `model = M3GNet.load()`
2. `relaxer = Relaxer()`
3. `relax_results = relaxer.relax(structure)`

## Performance Characteristics
- Fast compared to DFT (seconds vs hours).
- Accuracy reasonably close to DFT PBE for stability prediction.

## Application Areas
- High-throughput screening (filtering unstable structures before DFT).
- Phonon calculations (approximate).
- MD of complex multi-component systems.

## Community and Support
- Developed by Ong Group (UCSD) / Materials Virtual Lab
- Active development (MatGL is the successor)

## Verification & Sources
**Primary sources**:
1. GitHub: https://github.com/materialsvirtuallab/m3gnet
2. Publication: C. Chen and S. P. Ong, Nat. Comput. Sci. 2, 718 (2022)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN (GitHub)
- Development: ACTIVE
- Applications: Universal potential, GNN
