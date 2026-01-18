# pytopomat

## Official Resources
- **Homepage**: https://github.com/ncfrey/pytopomat
- **Repository**: https://github.com/ncfrey/pytopomat
- **License**: MIT License

## Overview
**pytopomat** is a Python package designed for the high-throughput screening of topological materials within the **Materials Project** ecosystem. It provides a user-friendly interface to apply the theory of **Topological Quantum Chemistry (TQC)** to crystal structures, allowing users to automatically diagnose whether a material is topological or trivial based on its symmetry and valence electron count, often without expensive surface state calculations.

**Scientific domain**: Materials Informatics, High-Throughput Topology
**Target user community**: Data Scientists, Materials Project users

## Theoretical Methods
- **Topological Quantum Chemistry (TQC)**: Mapping rules between atomic limits and band representations.
- **Symmetry Indicators**: Detection of topological bands via symmetry eigenvalues.
- **Fu-Kane Parity**: $\mathbb{Z}_2$ classification for centrosymmetric crystals.
- **Spillage Criteria**: Quantifying the change in wavefunction character due to Spin-Orbit Coupling.

## Capabilities
- **Screening**:
  - `is_topological(structure)`: Boolean check for topological bands.
  - `get_symmetries()`: Extracts relevant space group info.
- **Database Integration**:
  - Fetches pre-computed topology data from the Materials Project API.
  - Integration with `atomate` workflows for VASP.
- **Symmetry Analysis**:
  - Verification of Wyckoff positions and stoichiometry against TQC rules.

## Key Strengths
- **Informatics Ready**: Built on top of `pymatgen`, making it the easiest tool for users already comfortable with the Materials Project Python stack.
- **Scale**: Optimized for screening thousands of compounds, rather than deep-diving into one.
- **Simplicity**: Abstracts the heavy group theory of TQC into simple function calls.

## Inputs & Outputs
- **Inputs**: `pymatgen.Structure` objects.
- **Outputs**: Boolean flags, symmetry indicator lists.

## Interfaces & Ecosystem
- **Dependencies**: `pymatgen`, `spglib`, `numpy`.
- **Integration**: Materials Project REST API.

## Performance Characteristics
- **Speed**: Instantaneous for lookups or symmetry checks.
- **Bottleneck**: Depends on backend DFT if running new calculations.

## Comparison with Other Codes
- **vs. TopMat**: TopMat handles Magnetic Space Groups (1651). pytopomat focuses on the standard Space Groups (230) and non-magnetic materials.
- **vs. IrRep**: IrRep is the low-level calculator of traces. pytopomat is the high-level screening logic.

## Application Areas
- **Mining Databases**: Finding new candidate Topological Insulators in ICSD.
- **Dataset Creation**: Generating labels for machine learning models of topology.

## Community and Support
- **Development**: Nathan C. Frey (UPenn / Berkeley).
- **Source**: GitHub.

## Verification & Sources
- **Repository**: [https://github.com/ncfrey/pytopomat](https://github.com/ncfrey/pytopomat)
- **Verification status**: ✅ VERIFIED
  - Valid informatics tool.
