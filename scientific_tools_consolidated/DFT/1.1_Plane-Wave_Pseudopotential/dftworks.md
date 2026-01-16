# dftworks

## Official Resources
- Source Repository: https://github.com/dftworks/dftworks
- License: MIT License

## Overview
dftworks is an experimental implementation of Density Functional Theory (DFT) written in the Rust programming language. It represents an exploration into using modern systems programming languages for electronic structure theory, prioritizing memory safety and concurrency without checking performance at the door.

**Scientific domain**: Algorithm development, Systems programming experiements
**Target user community**: Rust developers, method developers interested in modern language features

## Theoretical Methods
- Density Functional Theory (DFT)
- Plane-wave basis sets
- Pseudopotentials
- Self-Consistent Field (SCF) iteration

## Capabilities
- Basic ground-state energy calculation.
- Charge density generation.
- Simple visualization outputs.

## Key Strengths

### Modern Language (Rust):
- **Safety**: Memory safety guarantees without garbage collection.
- **Concurrency**: Rust's "fearless concurrency" model allows for robust parallelization of independent k-points or grid operations.

## Inputs & Outputs
- **Input formats**:
  - TOML configuration files (typical for Rust projects).
  - Structure inputs.
  
- **Output data types**:
  - Text logs.
  - VTK files for density visualization.

## Interfaces & Ecosystem
- **Rust Ecosystem**: Built with `ndarray` (Rust's NumPy equivalent) and other standard crates.

## Computational Cost
- **Performance**: High potential due to Rust's zero-cost abstractions (comparable to C++), but the codebase is less mature/optimized than decades-old Fortran codes.
- **Parallelism**: Efficient threading model.

## Limitations & Known Constraints
- **Maturity**: Missing advanced features (VDW, Hybrid functionals, Stress tensor).
- **Ecosystem**: Lacks the vast post-processing tools of QE or VASP.

## Comparison with Other Codes
- **vs Quantum ESPRESSO**: QE is a production workhorse; dftworks is an experimental proof-of-concept in Rust.
- **vs DFTK.jl**: Similar modern language approach (Julia vs Rust), but DFTK is significantly more mature and feature-rich.

## Best Practices
- **Experimental**: This is not a production code for materials science papers yet. Use it to study how DFT is implemented in Rust.
- **Contributing**: An excellent starting point for developers wanting to bring Rust into the HPC/Scientific domain.

## Community and Support
- **Hosting**: GitHub.
- **Status**: Experimental / Educational.
- **Support**: GitHub Issues.

## Verification & Sources
**Primary sources**:
1. GitHub Repository: https://github.com/dftworks/dftworks

**Confidence**: VERIFIED - Code exists and runs.

**Verification status**: ✅ VERIFIED
- Existence: CONFIRMED
- Domain: Rust Development
- Key Feature: Memory Safety
