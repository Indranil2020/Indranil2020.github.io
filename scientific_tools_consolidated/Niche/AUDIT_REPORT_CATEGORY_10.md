# AUDIT REPORT: CATEGORY 10 (NICHE & MACHINE LEARNING)

## 1. Executive Summary
- **Category**: Niche & ML Tools
- **Total Files Audited**: 24 active files (plus moved/renamed items)
- **Status**: ✅ COMPLETED
- **Major Corrections**:
  - **Renamed**: `exactdiag` → `xdiag` (Confirmed active modern library).
  - **Clarified**: `Dual-fermions` updated to reflect it is a method with specific implementations (`opendf`, `DFermion`).
  - **Flagged**: `HubbardFermiMatsubara` marked as ⚠️ UNCERTAIN / RESEARCH CODE.
  - **Reorganized**: Moved `emmet`, `maggma`, `PyLada`, `MPWorks` to `Frameworks/` as they are core infrastructure.
  - **Moved**: `cmpy` to `DFT/1.6_Specialized`.

---

## 2. Detailed Audit & Corrections

### 2.1 Machine Learning Potentials & Tools

#### **MLIP** (Machine Learning Interatomic Potentials)
- **Status**: ✅ Exists
- **Path**: `Niche/MLIP.md`
- **Authenticity**: **VERIFIED** (GitLab: `shapeev/mlip-2`).
- **Accuracy**: Correctly identifies it as a Moment Tensor Potential (MTP) code.
- **Action**: Confirmed.

#### **n2p2** (Neural Network Potential Package)
- **Status**: ✅ Exists
- **Path**: `Niche/n2p2.md`
- **Authenticity**: **VERIFIED** (GitHub: `CompPhysVienna/n2p2`).
- **Accuracy**: Accurate description of Behler-Parrinello NN potentials and LAMMPS interface.
- **Action**: Confirmed.

#### **SIMPLE-NN**
- **Status**: ✅ Exists
- **Path**: `Niche/SIMPLE-NN.md`
- **Authenticity**: **VERIFIED** (GitHub: `MDIL-SNU/SIMPLE-NN`).
- **Accuracy**: Correctly identified as a user-friendly NN potential package.
- **Action**: Confirmed.

#### **AMP** (Atomistic Machine-learning Package)
- **Status**: ✅ Exists
- **Path**: `Niche/AMP.md`
- **Authenticity**: **VERIFIED** (Legacy/Stable).
- **Accuracy**: noted as legacy compared to newer equivariant models.
- **Action**: Confirmed.

#### **SchNetPack**
- **Status**: ✅ Exists
- **Path**: `Niche/SchNetPack.md`
- **Authenticity**: **VERIFIED** (GitHub: `atomistic-machine-learning/schnetpack`).
- **Accuracy**: Correctly identifies SchNet and PaiNN models.
- **Action**: Confirmed.

#### **MACE** (Multi-Atomic Cluster Expansion)
- **Status**: ✅ Exists
- **Path**: `Niche/MACE.md`
- **Authenticity**: **VERIFIED** (GitHub: `ACEsuit/mace`).
- **Accuracy**: Accurate description of higher-order equivariant message passing.
- **Action**: Confirmed.

#### **NequIP**
- **Status**: ✅ Exists
- **Path**: `Niche/NequIP.md`
- **Authenticity**: **VERIFIED** (GitHub: `mir-group/nequip`).
- **Accuracy**: Correctly identifies E(3)-equivariant nature and data efficiency.
- **Action**: Confirmed.

#### **Allegro**
- **Status**: ✅ Exists
- **Path**: `Niche/Allegro.md`
- **Authenticity**: **VERIFIED** (GitHub: `mir-group/allegro`).
- **Accuracy**: Correctly describes strict locality and scalability.
- **Action**: Confirmed.

#### **m3gnet**
- **Status**: ✅ Exists
- **Path**: `Niche/m3gnet.md`
- **Authenticity**: **VERIFIED** (GitHub: `materialsvirtuallab/m3gnet`).
- **Accuracy**: Accurate description of GNN capabilities and Materials Project integration.
- **Action**: Confirmed.

#### **Matbench**
- **Status**: ✅ Exists
- **Path**: `Niche/Matbench.md`
- **Authenticity**: **VERIFIED** (GitHub: `materialsproject/matbench`).
- **Accuracy**: Correctly identified as a benchmarking suite.
- **Action**: Confirmed.

#### **AFLOW-ML**
- **Status**: ✅ Exists
- **Path**: `Niche/AFLOW-ML.md`
- **Authenticity**: **VERIFIED** (AFLOW ecosystem).
- **Accuracy**: Describes API and property prediction capabilities.
- **Action**: Confirmed.

---

### 2.2 Specialized Physics Codes

#### **xdiag** (formerly exactdiag)
- **Status**: ✅ Exists (Renamed)
- **Path**: `Niche/xdiag.md`
- **Authenticity**: **VERIFIED** (GitHub: `awietek/xdiag`).
- **Correction**:
  - Original entry `exactdiag` pointed to a non-existent organization.
  - **FIXED**: Renamed to `xdiag` and updated links to the correct Alexander Wietek repository.
  - **Accuracy**: Validated features (Hubbard, t-J, Julia interface).

#### **EDRIXS**
- **Status**: ✅ Exists
- **Path**: `Niche/EDRIXS.md`
- **Authenticity**: **VERIFIED** (GitHub: `NSLS-II/edrixs`).
- **Accuracy**: Correctly describes RIXS/XAS simulation capabilities.
- **Action**: Confirmed.

#### **QMCPACK-addons** (Nexus)
- **Status**: ✅ Exists
- **Path**: `Niche/QMCPACK-addons.md`
- **Authenticity**: **VERIFIED** (Part of QMCPACK).
- **Accuracy**: Correctly describes workflow automation for QMC.
- **Action**: Confirmed.

---

### 2.3 Catalysis & Surface Science

#### **CatApp**
- **Status**: ✅ Exists
- **Path**: `Niche/CatApp.md`
- **Authenticity**: **VERIFIED** (SUNCAT web tool).
- **Action**: Confirmed.

#### **CatMAP**
- **Status**: ✅ Exists
- **Path**: `Niche/CatMAP.md`
- **Authenticity**: **VERIFIED** (GitHub: `SUNCAT-Center/CatMAP`).
- **Action**: Confirmed.

#### **GASpy**
- **Status**: ✅ Exists
- **Path**: `Niche/GASpy.md`
- **Authenticity**: **VERIFIED** (GitHub: `ulissigroup/GASpy`).
- **Action**: Confirmed.

---

### 2.4 Data & Repositories

#### **Zenodo**, **OSF**, **DataVerse**
- **Status**: ✅ Exist
- **Path**: `Niche/`
- **Authenticity**: **VERIFIED** (Global standard tools).
- **Action**: Confirmed.

---

### 2.5 Uncertain / Methodological Entries

#### **Dual-fermions**
- **Status**: ✅ Exists (Updated)
- **Path**: `Niche/Dual-fermions.md`
- **Issue**: "Dual-fermions" is a method, not a single code.
- **Correction**: Rewrote file to explain the method and point to specific implementations: **opendf** (CQMP) and **DFermion**.
- **Authenticity**: **VERIFIED** (as a method with valid implementations).

#### **HubbardFermiMatsubara**
- **Status**: ✅ Exists (Flagged)
- **Path**: `Niche/HubbardFermiMatsubara.md`
- **Issue**: No reliable public software package found.
- **Action**: Marked as ⚠️ **UNCERTAIN / RESEARCH CODE**. Advised users to look for alternatives.

---

### 2.6 Moved Files

- **emmet**, **maggma**, **PyLada**, **MPWorks**: Moved to `Frameworks/` to align with their function as infrastructure/frameworks.
- **cmpy**: Moved to `DFT/1.6_Specialized`.

---

## 3. Final Conclusion
The "Niche" category has been thoroughly audited. The Machine Learning section is particularly strong with up-to-date documentation on SOTA equivariant models (NequIP, MACE, Allegro). The ambiguity around generic names (`exactdiag`, `Dual-fermions`) has been resolved by identifying the specific modern libraries or clarifying the methodological nature of the entry.
