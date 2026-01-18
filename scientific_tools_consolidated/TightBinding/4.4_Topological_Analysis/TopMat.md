# TopMat

## Official Resources
- **Homepage**: https://github.com/zjwang11/TopMat
- **Repository**: https://github.com/zjwang11/TopMat
- **License**: GPL-3.0

## Overview
**TopMat** is a specialized high-throughput coding framework designed for the search and classification of **Magnetic Topological Materials**. Based on the theory of **Magnetic Topological Quantum Chemistry (MTQC)**, it extends topological indicators to all 1651 Magnetic Space Groups (MSGs). It was the engine behind the construction of the Magnetic Topological Materials Database, enabling the identification of hundreds of new magnetic topological insulators and semimetals.

**Scientific domain**: Magnetic Topology, High-Throughput Screening
**Target user community**: Materials Scientists, Database Builders

## Theoretical Methods
- **Magnetic Topological Quantum Chemistry (MTQC)**: Rigorous classification of bands based on corepresentations of magnetic space groups.
- **Symmetry Indicators**: Calculation of indices (e.g., $\mathbb{Z}_2, \mathbb{Z}_4$) that diagnose topology from symmetry eigenvalues.
- **Compatibility Relations**: Checking connectivity of bands in magnetic BZ.

## Capabilities
- **Classification**:
  - Automatically identifies the Magnetic Space Group from a spin-polarized VASP structure.
  - Computes Symmetry Indicators (SI) for occupied bands.
  - Assigns topology labels: Magnetic Topological Insulator (MTI), Magnetic Weyl Semimetal (MWSM), etc.
- **Workflow**:
  - Manages VASP calculations (SCF $\to$ Band structure).
  - Interfaces with `irvsp` to get traces.
- **Database**:
  - Tools to generate and query large SQL/JSON databases of materials.

## Key Strengths
- **Magnetic Completeness**: Unlike earlier non-magnetic classifiers, TopMat handles the full complexity of Type III and Type IV magnetic space groups.
- **Automation**: Designed for "high-throughput" from the start—can process thousands of POSCARs with minimal manual intervention.
- **Rigorous**: Identifies "Obstructed Atomic Limits" (OAL) which are often overlooked but crucial for higher-order topology.

## Inputs & Outputs
- **Inputs**:
  - VASP `POSCAR` and `INCAR` (defining magnetic configuration).
  - VASP `WAVECAR` (for symmetry analysis).
- **Outputs**:
  - Topological diagnosis labels.
  - Symmetry indicator values.

## Interfaces & Ecosystem
- **Dependencies**: `irvsp` (included/required), Python, VASP.
- **Database**: The output feeds into the Billings/MTQC database.

## Performance Characteristics
- **Speed**: The classification step is instantaneous ($< 1$ sec). The bottleneck is the underlying DFT calculation.
- **Scalability**: trivially parallelizable over materials.

## Comparison with Other Codes
- **vs. pytopomat**: `pytopomat` is primarily for non-magnetic materials (Time-Reversal Symmetric). `TopMat` is the specialist for Magnetic groups.
- **vs. CheckTop**: `CheckTop` is a web-interface. `TopMat` is the backend code.

## Application Areas
- **Axion Insulators**: Searching for MnBi2Te4-like materials.
- **Weyl Magnets**: Identifying time-reversal breaking Weyl points in Co-based shandites.

## Community and Support
- **Development**: Z.J. Wang Group (IOP CAS).
- **Source**: GitHub.

## Verification & Sources
- **Primary Publication**: Gao et al., Phys. Rev. B 106, 035150 (2022).
- **Verification status**: ✅ VERIFIED
  - Valid research workflow code.
