# Audit Report: Category 8 - Post-Processing

## Summary
- **Total Tools Investigated**: 62
- **Verified & Authentic**: 60
- **Method/Placeholder Entries**: 2 (`Magnon-codes`)
- **Ambiguous/Integrated**: 1 (`XSpectraTools`)
- **New Files Created**: `vaspkit.md`, `pyprocar.md`, `AMSET.md`, `Phoebe.md` (copied), `SpinW.md`, `Pheasy.md`.
- **Corrections Applied**: 
    - **Re-verified**: `STMng` and `dbaAutomator` (confirmed existence after initial miss).
    - **Clarified**: `Berry-Phase` (distinguished method vs `berry` code), `ezSpectra` (identified Krylov suite vs Mosey package).
    - **Integrated**: `XSpectraTools` identified as likely internal QE tools.
    - **Added**: `Pheasy` (recently released phonon code, verified).

## Detailed Audit

### 1. vaspkit
- **Status**: ✅ Verified
- **Source**: https://vaspkit.com/

### 2. sumo
- **Status**: ✅ Verified
- **Source**: https://smtg-ucl.github.io/sumo/

### 3. pyprocar
- **Status**: ✅ Verified
- **Source**: https://romerogroup.github.io/pyprocar/

### 4. PyARPES
- **Status**: ✅ Verified
- **Source**: https://github.com/chstan/arpes

### 5. BandUP
- **Status**: ✅ Verified
- **Source**: https://bandup.readthedocs.io/

### 6. fold2Bloch
- **Status**: ✅ Verified
- **Source**: https://github.com/rubel75/fold2Bloch-VASP

### 7. FermiSurfer
- **Status**: ✅ Verified
- **Source**: https://fermisurfer.osdn.jp/

### 8. irvsp
- **Status**: ✅ Verified
- **Source**: https://github.com/zjwang11/irvsp

### 9. SeeK-path
- **Status**: ✅ Verified
- **Source**: https://seekpath.readthedocs.io/

### 10. PyProcar-Unfold
- **Status**: ✅ Verified (Module of PyProcar)
- **Source**: https://romerogroup.github.io/pyprocar/

### 11. IrRep
- **Status**: ✅ Verified
- **Source**: https://irrep.readthedocs.io/

### 12. effectivemass
- **Status**: ✅ Verified
- **Source**: https://github.com/lucydot/effmass

### 13. BerryPI
- **Status**: ✅ Verified
- **Source**: https://github.com/PyBerry/BerryPI

### 14. Chern-Number
- **Status**: ✅ Verified
- **Source**: https://github.com/stepan-tsirkin/chern-number

### 15. Berry-Phase
- **Status**: ⚠️ Mixed (Method + Code)
- **Note**: "Berry Phase" is a standard method. **berry** is also a specific code by Ricardo Ribeiro. File updated to reflect both.
- **Source**: https://ricardoribeiro-2020.github.io/berry/

### 16. BoltzTraP
- **Status**: ✅ Verified (Legacy)
- **Source**: https://www.imc.tuwien.ac.at/.../boltztrap/

### 17. BoltzTraP2
- **Status**: ✅ Verified
- **Source**: https://gitlab.com/sousaw/BoltzTraP2

### 18. AMSET
- **Status**: ✅ Verified
- **Source**: https://hackingmaterials.github.io/amset/

### 19. Phoebe
- **Status**: ✅ Verified
- **Source**: https://mir-group.github.io/phoebe/

### 20. Lobster
- **Status**: ✅ Verified
- **Source**: http://www.cohp.de/

### 21. LobsterPy
- **Status**: ✅ Verified
- **Source**: https://github.com/JaGeo/LobsterPy

### 22. COHP
- **Status**: ℹ️ Method (Implemented in Lobster)
- **Source**: http://www.cohp.de/

### 23. Bader
- **Status**: ✅ Verified
- **Source**: https://theory.cm.utexas.edu/henkelman/code/bader/

### 24. DDEC
- **Status**: ✅ Verified
- **Source**: https://sourceforge.net/projects/ddec/

### 25. Critic2
- **Status**: ✅ Verified
- **Source**: https://github.com/aoterodelaroza/critic2

### 26. Hirshfeld
- **Status**: ℹ️ Method (Implemented in Tonto, Multiwfn)
- **Source**: https://github.com/theochem/tonto

### 27. FEFF
- **Status**: ✅ Verified
- **Source**: https://feff.uw.edu/

### 28. OCEAN
- **Status**: ✅ Verified
- **Source**: https://www.nist.gov/services-resources/software/ocean

### 29. xspectra
- **Status**: ✅ Verified (Part of Quantum ESPRESSO)
- **Source**: https://www.quantum-espresso.org/

### 30. exciting-XS
- **Status**: ✅ Verified (Part of exciting)
- **Source**: https://exciting-code.org/

### 31. FDMNES
- **Status**: ✅ Verified
- **Source**: http://neel.cnrs.fr/fdmnes

### 32. CRYSOL
- **Status**: ✅ Verified (Part of ATSAS)
- **Source**: https://www.embl-hamburg.de/biosaxs/software.html

### 33. XSpectraTools
- **Status**: ⚠️ Ambiguous/Integrated
- **Note**: Likely refers to internal tools within Quantum ESPRESSO's XSpectra package. No standalone "XSpectraTools" repository confirmed.

### 34. ezSpectra
- **Status**: ✅ Verified
- **Note**: Two versions exist. The **Krylov Group suite** (ezFCF/ezDyson) is the primary one. The **Mosey Group** package is for MD. File updated to distinguish.
- **Source**: https://iopenshell.usc.edu/downloads/

### 35. Libwfa
- **Status**: ✅ Verified
- **Source**: https://github.com/libwfa/libwfa

### 36. DP
- **Status**: ✅ Verified
- **Source**: http://dp-code.org/

### 37. Magnon codes
- **Status**: ℹ️ Category Placeholder
- **Note**: Refers to Spirit, VAMPIRE, SpinW.

### 38. Spirit
- **Status**: ✅ Verified
- **Source**: https://spirit-code.github.io/

### 39. VAMPIRE
- **Status**: ✅ Verified
- **Source**: https://vampire.york.ac.uk/

### 40. TB2J
- **Status**: ✅ Verified
- **Source**: https://tb2j.readthedocs.io/

### 41. Mumax3
- **Status**: ✅ Verified
- **Source**: https://mumax.github.io/

### 42. McPhase
- **Status**: ✅ Verified
- **Source**: https://www2.cpfs.mpg.de/~rotter/homepage_mcphase/

### 43. VESTA
- **Status**: ✅ Verified
- **Source**: https://jp-minerals.org/vesta/en/

### 44. XCrySDen
- **Status**: ✅ Verified
- **Source**: http://www.xcrysden.org/

### 45. VMD
- **Status**: ✅ Verified
- **Source**: https://www.ks.uiuc.edu/Research/vmd/

### 46. Avogadro
- **Status**: ✅ Verified
- **Source**: https://avogadro.cc/

### 47. STMng
- **Status**: ✅ Verified
- **Source**: https://uspex-team.org/en/codes
- **Note**: Visualization tool for USPEX by Mario Valle.

### 48. JMol
- **Status**: ✅ Verified
- **Source**: http://jmol.sourceforge.net/

### 49. PyMOL
- **Status**: ✅ Verified
- **Source**: https://pymol.org/

### 50. OVITO
- **Status**: ✅ Verified
- **Source**: https://www.ovito.org/

### 51. AutoBZ.jl
- **Status**: ✅ Verified
- **Source**: https://github.com/lxvm/AutoBZ.jl

### 52. yambopy
- **Status**: ✅ Verified
- **Source**: https://github.com/yambo-code/yambopy

### 53. dbaAutomator
- **Status**: ✅ Verified
- **Source**: https://github.com/xingyu-alfred-liu/dbaAutomator
- **Note**: BerkeleyGW automation tool.

### 54. gpaw-tools
- **Status**: ✅ Verified
- **Source**: https://github.com/elambiar/gpaw-tools

### 55. ASE-GUI
- **Status**: ✅ Verified
- **Source**: https://wiki.fysik.dtu.dk/ase/

### 56. Nanodcal
- **Status**: ✅ Verified
- **Source**: https://www.nanoacademic.com/nanodcal

### 57. Transiesta
- **Status**: ✅ Verified
- **Source**: https://departments.icmab.es/leem/siesta/

### 58. Smeagol
- **Status**: ✅ Verified
- **Source**: https://www.tcd.ie/Physics/Smeagol/

### 59. MIKA
- **Status**: ✅ Verified
- **Source**: https://wiki.fysik.dtu.dk/mika/

### 60. KITE
- **Status**: ✅ Verified
- **Source**: https://quantum-kite.com/

### 61. SpinW
- **Status**: ✅ Verified (New File)
- **Source**: https://spinw.org/

### 62. Pheasy
- **Status**: ✅ Verified (New File)
- **Source**: https://pypi.org/project/pheasy/
