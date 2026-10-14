# Angiosome and Filler Related Blindness

Convergence-controlled reproducibility package for the reduced-order CFD study.
The model is conditional on intravascular exposure and is not calibrated to clinical incidence.

## Reproduce the numerical analysis

Use Python 3.12 and install the pinned dependencies in `code/requirements.txt`.
From this directory run:

```sh
python -m pip install -r code/requirements.txt
python code/cfd_reproducibility_script.py --output reproduced
```

The script generates the same 12,000 baseline profiles (seed 42), 96,000 S0-S7 evaluations,
120,000 fixed-index choke evaluations, 2,400,000 threshold profile-setting classifications,
numerical checks, summaries and Figures 1-6. It raises an error on nonconvergence or
inconsistent exported classifications. Results are written to `reproduced/data` and
`reproduced/figures`. Numerical comparisons should use a tolerance, rather than file hashes,
because CSV formatting and PDF metadata can differ between software environments.

## Reference outputs

S0 OA entry is 0.241500 and CRA exposure is 0.20533333333333334.
S0 retains sampled upstream choke modification with coefficient 0.35.
The fixed-index grid uses index 0.4461147459514779 and a different zero-constriction reference:
OA entry 0.2653333333333333 and CRA exposure 0.22783333333333333.
All grid changes use that reference. The original eight-update results are superseded;
`data/numerical_correction_record.json` documents the correction using identical inputs.

The five-tab Excel workbook contains all baseline inputs and outputs, all scenario records,
the eight scenario summaries, the ten choke summaries and all choke profile records.
The scenario and choke workbook tabs contain the essential numerical outputs; the compressed
CSV files retain every exported field. `profile_id` joins those records to the baseline inputs.
Scenario overrides are specified in the executable code and manuscript Table 2.
The full 25-combination threshold audit for all eight scenarios is supplied as CSV.
The combined CSV has 228,018 rows and a `record_type` column; these are records from distinct
populations and summaries, not 228,018 independent patients. Blank union-schema fields mean
that the field does not apply to that record type.

Large CSV files are gzip compressed without altering their contents. Decompress them using
a standard archive utility, or read them directly with `pandas.read_csv`.

## Verification and limits

Verification includes a 1,000-profile analytical constant-viscosity Poiseuille limit,
an independent Brent root solver in 200 baseline profiles, a 10^-8 versus 10^-10 tolerance
comparison, closure checks, volume conservation and classification reconciliation.
The input-parameter random-forest audit is held-out permutation screening, not Sobol analysis.
The separate threshold audit preserves scenario ordering across all 25 combinations.
No external physiological validation has been performed. The downstream implementation
changes a terminal-volume condition after the OA decision; it is not a coupled downstream
resistance model. Geometry, pressure coupling and several constitutive closures are
prescribed computational assumptions. No delayed venous route is solved.

## Repository contents

This repository contains code and numerical data only. Manuscript documents and figure images are not included. The code can generate figures locally when executed. Large outputs are gzip-compressed CSV files; the Excel workbook contains five tabs.

## Version provenance

The current branch contains the convergence-controlled calculation and supersedes the original fixed eight-cycle outputs. Previous GitHub releases and Zenodo records are archival versions and have not been replaced by this branch update. A new matching release and archive are needed before citing these updated data through a DOI.
