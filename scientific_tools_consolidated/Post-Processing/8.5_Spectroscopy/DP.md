# DP (Dielectric Properties)

## 1. Description
**DP** is a code for calculating the linear response dielectric properties of periodic systems. It uses Time-Dependent Density Functional Theory (TDDFT) in the frequency domain with a plane-wave basis set. It is designed to compute the macroscopic dielectric function, Electron Energy Loss Spectra (EELS), and Inelastic X-ray Scattering (IXS) spectra.

## 2. Capability Verification
- **Functionality**:
  - Calculation of frequency-dependent dielectric function $\epsilon(\omega)$
  - Electron Energy Loss Spectroscopy (EELS)
  - Inelastic X-ray Scattering (IXS)
  - Local Field Effects (LFE) inclusion
  - Crystal local field effects
- **Key Features**:
  - Interface with ABINIT (reads WFK files)
  - Plane-wave basis set
  - Linear response TDDFT (RPA, ALDA, LRC kernels)

## 3. Authenticity & Usage
- **Official Website**: [http://dp-code.org/](http://dp-code.org/)
- **Source Code**: Available via website (GPL)
- **Developers**: The ETSF (European Theoretical Spectroscopy Facility) DP team (V. Olevano, L. Reining, G. Siny)

## 4. Technical Assessment
- **Status**: Active (Legacy/Stable)
- **Reliability**: High (Standard tool in the ETSF community)
- **Integration**: Strong coupling with ABINIT
