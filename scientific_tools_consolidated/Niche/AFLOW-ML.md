# AFLOW-ML (AFLOW Machine Learning)

## Official Resources
- Homepage: http://aflow.org/aflow-ml/
- Documentation: http://aflow.org/aflow-ml/
- Source Repository: Part of AFLOW codebase
- License: GPL v3

## Overview
AFLOW-ML is a machine learning API and library integrated into the AFLOW framework. It provides access to pre-trained machine learning models for predicting materials properties (electronic, thermal, mechanical) based on crystal structure and composition. It allows users to screen materials without performing expensive DFT calculations.

**Scientific domain**: Materials informatics, machine learning, property prediction  
**Target user community**: Materials scientists, high-throughput researchers

## Capabilities (CRITICAL)
- **Property Prediction**: Predicts band gap, bulk/shear modulus, Debye temperature, heat capacity, thermal conductivity, etc.
- **Models**: Uses Gradient Boosting Decision Trees (GBDT), Voronoi tessellation features, and PLMF (Property Labeled Materials Fragments).
- **API**: REST API for programmatic access to predictions.
- **Online Interface**: Web-based predictor.

**Sources**: AFLOW-ML website, Sci. Rep. 7, 10766 (2017)

## Inputs & Outputs
- **Input formats**: POSCAR, Composition, Structure JSON
- **Output data types**: Predicted property values (JSON)

## Interfaces & Ecosystem
- **AFLOW**: Integrated with the main AFLOW framework.
- **Python**: API client available.

## Workflow and Usage
1. Prepare structure (POSCAR).
2. Send to API:
   ```bash
   curl -X POST -d @POSCAR http://aflow.org/API/aflow-ml/v1.0/plmf/v1.0/
   ```
3. Receive JSON response with predictions.

## Performance Characteristics
- Extremely fast (milliseconds per prediction)
- Accuracy depends on training set coverage (AFLOW database)

## Application Areas
- Rapid screening of millions of compounds
- Guiding DFT calculations
- Discovery of superhard materials

## Community and Support
- Developed by Curtarolo Group (Duke)
- Active development

## Verification & Sources
**Primary sources**:
1. Homepage: http://aflow.org/aflow-ml/
2. Publication: O. Isayev et al., Nat. Commun. 8, 15679 (2017)

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: AVAILABLE
- Source: OPEN (GPL)
- Development: ACTIVE
- Applications: Machine learning, property prediction
