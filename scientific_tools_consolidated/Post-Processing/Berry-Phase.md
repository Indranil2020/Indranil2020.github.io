# Berry-Phase / berry

## 1. Ambiguity Resolution
"Berry-Phase" primarily refers to the **method** of calculating the geometric phase in electronic structure theory (Modern Theory of Polarization). However, there is also a specific software package named **berry** designed for this purpose.

## 2. Specific Software: **berry**
- **Official Resources**:
  - **Homepage**: https://ricardoribeiro-2020.github.io/berry/
  - **Source Repository**: https://github.com/ricardoribeiro-2020/berry
  - **Publication**: *Computer Physics Communications* 267, 108064 (2021).
- **Description**: A code for the differentiation of Bloch wavefunctions from DFT calculations. It calculates Berry connections, Berry curvature, and topological invariants (Chern numbers, Z2) by unwinding the phase of the Bloch states.
- **Capabilities**:
  - Second Harmonic Generation (SHG) conductivity
  - Anomalous Hall Conductivity
  - Topological invariants
  - Interface with Quantum ESPRESSO and Wannier90.

## 3. General Method Capabilities
Most major DFT codes implement Berry phase calculations natively for polarization:
- **VASP**: `LBERRY = .TRUE.` (Calculates polarization via Berry phase).
- **Quantum ESPRESSO**: `bp` calculation (Berry phase).
- **ABINIT**: Berry phase polarization.
- **Wannier90**: Calculates Berry curvature and anomalous Hall conductivity via Wannier interpolation.

## 4. Related Tools
- **BerryPI**: A Python wrapper for VASP to calculate Berry phases ([Link](https://github.com/PyBerry/BerryPI)).
- **PythTB**: Python Tight Binding (for model Hamiltonian Berry phase).
- **Z2Pack**: For topological invariant calculations.
