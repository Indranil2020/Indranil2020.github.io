# ezSpectra

## 1. Ambiguity Resolution
There are two distinct software packages with the name **ezSpectra**:
1.  **ezSpectra Suite (Krylov Group)**: A comprehensive toolkit for modeling electronic spectroscopy (photoelectron, absorption).
2.  **ezSpectra (Mosey Group)**: A Python package for calculating vibrational spectra from MD trajectories.

## 2. ezSpectra Suite (Primary/Major)
- **Official Resources**:
  - **Homepage**: https://iopenshell.usc.edu/downloads/
  - **Reference**: *WIREs Comput. Mol. Sci.* 12, e1546 (2022).
  - **Developers**: Anna Krylov Group (USC).
- **Components**:
  - **ezFCF**: Franck-Condon factors for vibronic spectra.
  - **ezDyson**: Photoionization cross-sections using Dyson orbitals.
- **Capabilities**: Simulating photoelectron, photodetachment, and absorption spectra with vibronic structure.

## 3. ezSpectra (Mosey Group)
- **Official Resources**:
  - **Repository**: https://github.com/mosey-group/ezSpectra
  - **License**: MIT
- **Capabilities**:
  - Calculation of IR, Raman, and VCD spectra from dipole/polarizability trajectories obtained from molecular dynamics (e.g., CP2K).
  - Signal processing and plotting.

## 4. Recommendation
Users should identify which tool matches their domain (Electronic Structure/Quantum Chemistry vs. Molecular Dynamics/Vibrational Spectroscopy).
