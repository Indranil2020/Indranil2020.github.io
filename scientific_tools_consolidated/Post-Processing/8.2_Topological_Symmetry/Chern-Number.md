# Chern-Number

## 1. Description
**Chern-Number** refers to a code for calculating the Chern number, a topological invariant, using the discretized Berry curvature on a grid. This specific entry corresponds to the tool developed by Stepan Tsirkin (University of Zurich), which is designed to post-process output from Wannier90 or tight-binding models.

## 2. Capability Verification
- **Functionality**:
  - Calculation of Berry curvature in the Brillouin zone
  - Integration of curvature to obtain the Chern number
  - Handling of metallic and insulating systems
  - Detection of topological phase transitions
- **Key Features**:
  - Grid-based integration method
  - Interface with Wannier90 output (`_hr.dat`)
  - Efficient Fortran/Python implementation

## 3. Authenticity & Usage
- **Source Code**: [https://github.com/stepan-tsirkin/chern-number](https://github.com/stepan-tsirkin/chern-number)
- **Developer**: Stepan Tsirkin
- **Related Tools**: `IrRep`, `BerryPI` (same developer ecosystem)

## 4. Technical Assessment
- **Status**: Research code
- **Reliability**: Verified for topological physics applications
- **Ease of Use**: Requires knowledge of topological band theory
