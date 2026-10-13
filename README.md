# Reduced-order CFD network model of filler-related visual loss

This repository contains the reproducibility package for a reduced-order computational fluid dynamics network model of filler-related visual loss.

## Contents

- `manuscript/`
  - Final manuscript DOCX.
- `supplement/`
  - Supplementary methods DOCX.
- `data/`
  - Baseline simulation outputs.
  - Scenario evaluations.
  - Stress-test outputs.
  - Choke-sensitivity outputs.
  - Combined supplementary data CSV.
  - Five-tab Excel supplementary data workbook.
- `code/`
  - Python reproducibility script used to generate the virtual cohort, solve the CFD network, create figures and export CSV outputs.
- `figures/`
  - Figure PNG files generated from the computational outputs.

## Reproducibility

The computational model uses a fixed random seed and generates 12,000 baseline intravascular exposure profiles. These profiles are re-evaluated across eight prespecified CFD scenarios, giving 96,000 scenario evaluations. Choke-resistance sensitivity is reported separately across upstream and downstream implementations.

The model is conditional on intravascular exposure and does not estimate clinical incidence.

## Recommended citation

Rahman E. A reduced-order CFD network model of filler-related visual loss. Repository release. DOI to be added after Zenodo archival.

## Data and code availability statement

Simulation outputs and the reproducibility script are provided in this repository. Before journal submission, create a GitHub release and archive it using Zenodo to obtain a permanent DOI.

