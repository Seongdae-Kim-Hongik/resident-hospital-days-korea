# Analysis code: deprivation, physical activity and resident inpatient days in Korean districts, 2011-2024

Code-only repository (version 1.0.0).

This repository contains the data-processing and analysis code for an ecological spatio-temporal analysis of sex-and-age-standardised inpatient days of residents in 229 districts of the Republic of Korea, 2011-2024. No data, results or manuscript files are included. All inputs are publicly released aggregate statistics.

## Software
Python 3.9 (pandas, numpy, geopandas, statsmodels, pymannkendall), R 4.x with INLA, sf and spdep.

## Environment variables

- `ST_REF_DATASET`: path to the reference SQLite dataset of district-year indicators (public statistics compiled by district and year)
- `ST_FBDZIP`: folder with resident registration population by district, sex and 5-year age group, one CSV per year
- `ST_NHIS_NATIONAL`: CSV of national National Health Insurance inpatient days by sex and 5-year age group (KOSIS TX_35001_A047, 2011-2024)
- `ST_SHP`: district boundary shapefile (229 districts)
- `ST_GRAPH`: INLA graph file of the 593-edge contiguity graph
- `ST_DATA_ROOT`: folder holding auxiliary public files referenced by the shared scripts
- `ST_THREADS`: INLA threads (1 used)
- `SK_RSLOTS`: parallel R jobs (2 used)

## Order of execution

| Step | Script | Purpose |
|---|---|---|
| 1 | `01_스크립트/99L_build_Estd_inpt.py` | Outcome (inpatient days of residents: in-district plus out-of-district admissions), structural-zero check, resident population by sex and 5-year age, indirectly standardised expected counts (national rates by sex and age, rescaled to observed annual totals) |
| 2 | `04_분석/INPT/gate/code/06_village_corrections.py` | Corrections of village sports facility counts with verified source-table row errors |
| 3 | `04_분석/INPT/build_pool_from_sqlite_inpt.py, 04_분석/INPT/build_pool_inpt_additions.py` | District-year covariate pool from the reference dataset |
| 4 | `04_분석/INPT/gate/code/10_inpt_cell_gate.py` | Cell-level checks (year jumps, boundary changes, unit breaks) |
| 5 | `04_분석/INPT/make_pool_sets_inpt.py` | Model settings sets_INPT_<version>.json and hypothesis table, hashed before fitting |
| 6 | `ST_VER=<version> python3 01_스크립트/80_build_disease_panel_v.py INPT` | Analytic panel for the main version (std3) and the eight sensitivity versions |
| 7 | `04_분석/INPT/make_nco_perio.py` | Negative-control outcome panel (periodontal care) used by 88g |
| 8 | `04_분석/INPT/run_inpt_main.sh, then 04_분석/INPT/run_inpt_rest.sh` | Main model, multiverse, sensitivity versions, period models and maps, emerging hot spots, sensitivity suite, posterior predictive checks, structural confounding |
| 8b | `04_분석/INPT/make_sets_posthoc_ma.py, then ST_VER=ma python3 01_스크립트/80_build_disease_panel_v.py INPT and ST_VER=ma ST_MV_N=2 ST_SENS_LITE=1 Rscript 01_스크립트/81_fit_inla_inpt.R INPT` | Post hoc sensitivity model adding Medical Aid type 1 patients per 100,000 as a structural covariate |
| 9 | `04_분석/INPT/90_inpt_summary_numbers.py` | Robustness verdicts, E-values, false discovery rates and the numbers file used for all reported values |

## Data sources
National Health Insurance district medical-use statistics (in-district and out-of-district inpatient days by district of residence), Korea Community Health Survey district estimates, resident registration population (Statistics Korea, KOSIS), national inpatient days by sex and age (KOSIS table TX_35001_A047) and other national district indicators.

## Korean identifiers
Folder names (`01_스크립트/` scripts, `02_패널/` panels created by the code, `04_분석/` analysis), some column and region names, variable labels and console messages are in Korean because the source tables and the shared pipeline are Korean; they are left unchanged so that the code runs on the original files. Comments are in English or were removed.

## Licence
MIT (see LICENSE).

## Citation
See CITATION.cff.
