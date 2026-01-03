# SpinW

## Official Resources
- **Homepage**: https://spinw.org/
- **Source Repository**: https://github.com/SpinW/spinw
- **Documentation**: https://spinw.org/spinwdoc/
- **License**: MIT License (v4.0+) / CC-BY-NC-SA (older versions)

## Overview
SpinW is a MATLAB library (with upcoming Python/C++ versions) for optimizing magnetic structures and calculating magnetic excitations (magnons) using Linear Spin Wave Theory (LSWT). It is widely used for fitting experimental inelastic neutron scattering (INS) data to spin Hamiltonian models.

**Scientific domain**: Magnetism, Neutron Scattering, Spin Dynamics
**Target user community**: Neutron scattering experimentalists, magnetism theorists

## Capabilities
- **Magnetic Structure**: Optimization of magnetic structures (simulated annealing).
- **Excitons**: Calculation of spin wave dispersion relations and spectral functions S(Q,w).
- **Fitting**: Fitting of model parameters (exchange J, anisotropy K, DMI) to experimental INS data.
- **Symmetry**: Analysis of magnetic symmetry groups.
- **Visualisation**: 3D plotting of magnetic structures and dispersion surfaces.

## Inputs & Outputs
- **Input**: Crystal structure (CIF), Hamiltonian parameters.
- **Output**: Dispersion plots, Neutron scattering cross-sections, formatted figures.

## Interfaces & Ecosystem
- **MATLAB**: Runs within the MATLAB environment.
- **Horace**: Interoperable with the Horace software for neutron data analysis.
- **Python**: `spinw` python interface (experimental).

## Verification & Sources
- **Primary Source**: [SpinW Homepage](https://spinw.org/)
- **Publication**: S. Toth and B. Lake, *J. Phys.: Condens. Matter* 27, 166002 (2015).
