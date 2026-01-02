# Directory Structure Mapping to Master List

## Category 1: Ground-State DFT (80 codes)
**Directory**: `/scientific_tools_consolidated/DFT/`
**Files**: 80 markdown files

### 1.1 Plane-Wave / Pseudopotential Codes (23 tools)
- VASP.md
- Quantum-ESPRESSO.md
- ABINIT.md
- CASTEP.md
- CP2K.md *(also has TDDFT)*
- CPMD.md
- GPAW.md *(also has TDDFT)*
- JDFTx.md
- Qbox.md *(also has TDDFT)*
- PARSEC.md
- PARATEC.md
- SPARC.md
- RMGDFT.md
- ABACUS.md
- ATOMPAW.md *(dataset generator)*
- PROFESS.md
- MADNESS.md
- OpenAtom.md
- PWDFT.md
- PLATO.md
- NESSIE.md
- DFT-FE.md

### 1.2 All-Electron Codes (14 tools)
- WIEN2k.md
- Elk.md
- Fleur.md
- exciting.md *(also has TDDFT)*
- Questaal.md
- RSPt.md
- KKR.md
- JuKKR.md
- KKRnano.md
- KKRhost.md
- NRG-ETH.md *(DMFT impurity solver)*
- FPLO.md
- NRG-ETH-CSC.md *(DMFT impurity solver)*
- KKR-ASA.md

### 1.3 Localized Basis Sets (10 tools)
- FHI-aims.md
- SIESTA.md
- OpenMX.md
- CONQUEST.md
- ONETEP.md
- BigDFT.md
- CRYSTAL.md
- ADF.md
- DMol3.md *(commercial)*
- deMon2k.md

### 1.4 Quantum Chemistry Suites (26 tools)
- ORCA.md
- Gaussian.md
- PySCF.md
- PSI4.md
- Molpro.md
- NWChem.md *(also has TDDFT)*
- Turbomole.md
- Q-Chem.md
- GAMESS.md
- Dalton.md
- DIRAC.md
- CFOUR.md
- MRCC.md
- OpenMolcas.md
- BAGEL.md
- Columbus.md
- ACES.md
- ExaChem.md
- Quantum-Package.md
- CheMPS2.md
- SlowQuant.md
- BDF.md
- eT.md
- CC4S.md
- ACES-III.md

### 1.5 Tight-Binding DFT (5 tools)
- DFTB+.md
- xTB.md
- HOTBIT.md
- MOPAC.md
- AMS-DFTB.md

### 1.6 Specialized (2 tools)
- FlapwMBPT.md
- cmpy.md *(Python utilities library)*

---

## Category 2: TDDFT & Excited-State (19 unique codes)
**Directory**: `/scientific_tools_consolidated/TDDFT/`
**Files**: 19 markdown files

### TDDFT-Specific Codes (19 tools)
- Octopus.md
- SALMON.md
- Yambo.md
- turboTDDFT.md *(legacy QE module)*
- PyTDDFT.md *(research prototype)*
- TDAP.md *(uncertain)*
- BerkeleyGW.md
- WEST.md
- Spex.md
- SternheimerGW.md *(research code)*
- Fiesta.md
- molgw.md
- GreenX.md *(library)*
- SAX.md *(uncertain)*
- OCEAN.md
- NBSE.md *(uncertain)*
- DP-Code.md *(research)*
- DP-4.md *(research)*
- pyGWBSE.md *(Python educational)*

### Overlapping Codes (5 tools documented in DFT directory)
**Note**: These DFT codes also have TDDFT capabilities (master list #105-109)
- Qbox (in DFT/Qbox.md)
- NWChem (in DFT/NWChem.md)
- CP2K (in DFT/CP2K.md)
- GPAW (in DFT/GPAW.md)
- exciting (in DFT/exciting.md)

---

## Total Unique Codes: 99
- Category 1 unique: 80 codes
- Category 2 unique: 19 codes
- Overlapping (counted once): 5 codes in both categories

## Files Removed (Not Actual Codes)
From Category 1:
- Materials-Studio (GUI platform)
- Medea (GUI platform)
- FLAPW (method, not software)
- Molcas (superseded by OpenMolcas)
- DFT-F (typo/duplicate of DFT-FE)
- RMG (duplicate of RMGDFT)
- TURBOMOLE (case duplicate of Turbomole)
- PWPAW (not in master list)
- TBPW (not in master list)
- LMTO, LMTO-ASA (not in corrected master list)

## Next Category
**Category 3: DMFT & Many-Body (49 tools)**
- Directory to create: `/scientific_tools_consolidated/DMFT/`
- Subcategories:
  - DMFT Frameworks
  - Impurity Solvers
  - Many-Body Codes
  - Post-processing Tools
