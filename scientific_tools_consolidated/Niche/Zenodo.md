# Zenodo

## Official Resources
- Homepage: https://zenodo.org/
- Documentation: https://help.zenodo.org/
- Source Repository: https://github.com/zenodo/zenodo
- License: Open Source (Code), Open Data

## Overview
Zenodo is a general-purpose open-access repository developed under the European OpenAIRE program and operated by CERN. It allows researchers to deposit data sets, research software, reports, and any other research related digital artifacts. For each submission, a persistent identifier (DOI) is minted, making the stored items easily citeable.

**Scientific domain**: General research repository, open science  
**Target user community**: All researchers

## Capabilities (CRITICAL)
- **DOI Minting**: Automatic DOI assignment for every upload.
- **GitHub Integration**: Automatically archives GitHub releases and assigns DOIs to code.
- **Versioning**: Supports versioning of datasets and software.
- **Communities**: Users can create "Communities" to curate collections (e.g., "Materials Science").
- **Long-term Retention**: Hosted at CERN Data Center.

**Sources**: Zenodo website

## Inputs & Outputs
- **Input formats**: Any file format (up to 50GB per record).
- **Output data types**: DOI, landing page.

## Interfaces & Ecosystem
- **GitHub**: Webhooks for auto-archiving.
- **REST API**: For automated uploading/downloading.
- **OAI-PMH**: Metadata harvesting.

## Workflow and Usage
1. Go to zenodo.org.
2. Upload files (data, code).
3. Fill metadata.
4. Publish -> Get DOI.
5. Cite DOI in paper.

## Performance Characteristics
- Reliable, long-term storage.
- Free.

## Application Areas
- Publishing supplementary data for papers.
- Archiving research software.
- Sharing training datasets for ML.

## Community and Support
- CERN / OpenAIRE.
- Global standard for general data sharing.

## Verification & Sources
**Primary sources**:
1. Homepage: https://zenodo.org/

**Confidence**: VERIFIED

**Verification status**: ✅ VERIFIED
- Website: ACTIVE
- Documentation: COMPREHENSIVE
- Source: OPEN
- Development: ACTIVE
- Applications: Data repository, DOI
