# Code for: A transportability-aware triangulation of metabolic stroke evidence across the Shanghai Cooperation Organisation

Companion analysis code for the manuscript submitted to *European Journal of
Epidemiology* (Wang X, Li M, Zhang J, Zhang H, Ma Y).

The study triangulates three evidence streams for five metabolic risk factors of
stroke across the ten member states of the Shanghai Cooperation Organisation:

1. **Population stream** — GBD 2023 attributable fractions of stroke DALYs.
2. **Individual stream** — CHARLS, a Chinese national cohort (Fine–Gray
   competing-risk models).
3. **Genetic stream** — two-sample univariable and multivariable Mendelian
   randomisation in European, East Asian and South Asian ancestry.

The three streams are compared qualitatively under a prespecified rule set
(mechanisms M1–M4, with the final M2 rule); see `rules/`.

---

## Repository layout

```
code/
├── config.R                  # all paths are defined here
├── README.md
├── 01–12_*.R                 # GBD component
├── 20–34_*.R                 # CHARLS component
├── 40–46_*.R                 # Mendelian randomisation component
├── 50–53_*.R                 # cross-stream integration and rule application
├── 60–63_*.R                 # figures
└── rules/                    # M2 rule and mechanism archive
```

Scripts are numbered in the intended run order within each component. Scripts
in the CHARLS component share `21_charls_model_utils.R`; the integration
scripts consume outputs written by the earlier steps.

## Requirements

- R 4.5.0 (the version used for the manuscript).
- Packages: `tidyverse`, `data.table`, `ggplot2`, `ragg`, `patchwork`,
  `cowplot`, `ggrepel`, `ggnewscale`, `ggh4x`, `scales`, `reshape2`,
  `survival`, `cmprsk`, `rms`, `Hmisc`, `tableone`, `rstatix`, `gtsummary`,
  `pROC`, `xgboost`, `SHAPforxgboost`, `glmnet`, `caret`,
  `TwoSampleMR`, `MVMR`, `MendelianRandomization`,
  `cluster`, `fpc`, `mclust`, `dendextend`, `factoextra`, `ggdendro`,
  `compositions`, `vegan`, `spdep`, `geosphere`, `igraph`,
  `sf`, `rnaturalearth`, `ComplexHeatmap`, `circlize`, `mgcv`, `forestploter`.
  Several scripts install a missing package on first use; installing them
  up front is nevertheless recommended.

## Data

No data are distributed with this repository. The code expects the following
layout, which mirrors the authors' working directory:

```
<DATA_ROOT>/
├── rawdat/                  # GBD extracts, CHARLS files, intermediate results
└── GWAS/                    # GWAS summary statistics for the MR scripts
```

`<DATA_ROOT>` is set in `config.R` (default: `../data` relative to this
repository) or through the environment variable `GBDSTROKE_DATA_ROOT`.

Sources of the input data:

| Input | Source | Access |
|---|---|---|
| GBD 2023 estimates | GBD Results tool, Institute for Health Metrics and Evaluation | public, registration |
| CHARLS (Harmonized CHARLS, wave-1/3 blood subsample) | CHARLS data portal | application and approval |
| European exposure/outcome GWAS | OpenGWAS, FinnGen, UK Biobank–FinnGen meta-analyses | public (accession IDs in the scripts) |
| East Asian exposure/outcome GWAS | Taiwan Precision Medicine Initiative PheWeb, BioBank Japan PheWeb, China Kadoorie Biobank PheWeb | public |
| South Asian exposure/outcome GWAS | Pan-UK Biobank, MAGIC, GIGASTROKE | public |

The individual-level CHARLS data cannot be redistributed; the code documents
every file it reads so that the analysis can be repeated by researchers with
approved access. Accession identifiers for every GWAS dataset are given at the
top of the corresponding scripts.

## Configuration

All paths are resolved in `config.R`:

```r
CODE_DIR   # this repository
DATA_ROOT  # folder containing rawdat/ and GWAS/
RAWDAT     # DATA_ROOT/rawdat
GWAS_DIR   # DATA_ROOT/GWAS
```

Set the data root either by editing `config.R` or, preferably, in the
environment:

```bash
export GBDSTROKE_DATA_ROOT=/path/to/data
```

OpenGWAS requires a personal token. Export it before running the MR scripts;
never commit a token to the repository:

```bash
export OPENGWAS_JWT=<your token>
```

## Running the code

Run each script from this directory so that `config.R` is found:

```bash
cd code
Rscript 01_gbd_burden_trends.R
```

Recommended order:

| Step | Scripts |
|---|---|
| GBD | `01` → `02` → `03` → `04`; `05` → `06` → `07` → `08` → `09`; `10` → `11`; `12` |
| CHARLS | `20` → `21` (utility) → `22` → `23` → `24` → `25` → `26` → `27` → `28` → `29` → `30` → `31` → `32` → `33` → `34` |
| MR | `40`–`42` (univariable, by ancestry) → `43`–`45` (multivariable) → `46` (forest plot) |
| Integration | `50` or `51` (final version) → `52` → `53` |
| Figures | `60` → `61`; `62`; `63` |

`50_b5_integration.R` and `51_b5_integration_v9.R` are successive versions of
the same pipeline; **`51` is the final version** used for the manuscript. `50`
is retained for provenance.

Outputs are written next to the inputs under `<DATA_ROOT>` (for example
`rawdat/CHARLS/`, `rawdat/PAF/`, `B5_outputs/`), matching the paths used by the
scripts.

## Script index

| Repository file | Original analysis file | Purpose |
|---|---|---|
| `01_gbd_burden_trends.R` | `GBD-year.R` | Stroke burden trends, 2001–2023 |
| `02_gbd_trend_cluster_stability.R` | `GBD-year_stability.R` | Stability of burden-trajectory clusters |
| `03_gbd_trend_spatial_tests.R` | `GBD-year_空间检验.R` | Spatial tests of the burden clustering |
| `04_gbd_age_sex.R` | `GBD-age&sex.R` | Age- and sex-specific burden |
| `05_paf_attribution.R` | `PAF_v2.R` | Population attributable fractions and attribution clustering |
| `06_paf_cluster_stability.R` | `PAF_v2_stability.R` | Stability of attribution clusters |
| `07_paf_spatial_tests.R` | `PAF_v2_空间检验.R` | Spatial tests of attribution |
| `08_paf_subgroups.R` | `PAF_subgroup(age&sex).R` | PAF by age and sex |
| `09_paf_c2_sensitivity.R` | `PAF_v2_C2_sensitivity.R` | Sensitivity of the C1/C2 contrast |
| `10_subtype_composition.R` | `亚型占比.R` | Stroke-subtype composition |
| `11_subtype_cluster_stability.R` | `亚型占比_stability.R` | Stability of subtype clusters |
| `12_sdi_gam_frontier.R` | `前沿分析_GAM.R` | SDI–burden GAM and efficiency frontier |
| `20_charls_prepare_data.R` | `CHARLS整理改.R` | CHARLS cohort construction |
| `21_charls_model_utils.R` | `charls_model_utils.R` | Shared modelling utilities |
| `22_charls_baseline_table.R` | `CHARLS_汇总表.R` | Baseline characteristics (Table 2) |
| `23_charls_pairwise_smd.R` | `CHARLS_汇总pairwise.R` | Pairwise SMDs (Table S28) |
| `24_charls_fine_gray.R` | `Fine-Grey改.R` | Main Fine–Gray competing-risk models |
| `25_charls_rcs.R` | `RCS_FineGray.R` | Restricted cubic splines |
| `26_charls_subgroups.R` | `CHARLS_亚组_wholeSD.R` | Subgroup analyses |
| `27_charls_cause_specific_cox.R` | `CHARLS_敏感-cox.R` | Cause-specific Cox sensitivity analysis |
| `28_charls_dilution_prepare.R` | `CHARLS_稀释敏感整理.R` | Regression-dilution data preparation |
| `29_charls_dilution_fine_gray.R` | `Fine-Gre_稀释敏感.R` | Two-wave-averaged exposure analysis |
| `30_charls_cystatin_c.R` | `CHARLS_敏感_cys.R` | Creatinine–cystatin C eGFR sensitivity |
| `31_charls_shap.R` | `SHAP_改.R` | Gradient-boosted importance (XGBoost–SHAP) |
| `32_charls_cif.R` | `CHARLS_CIF.R` | Cumulative incidence functions (Fig. S18) |
| `33_charls_flowchart.R` | `CHARLS_flowchart.R` | Cohort flow chart (Fig. S17) |
| `34_charls_forest.R` | `CHARLS_forest.R` | SHR forest plot (Fig. 5) |
| `40_mr_eur.R` | `MR_EUR改.R` | Univariable MR, European ancestry |
| `41_mr_eas.R` | `MR_EAS改.R` | Univariable MR, East Asian ancestry |
| `42_mr_sas.R` | `MR_SAS改.R` | Univariable MR, South Asian ancestry |
| `43_mvmr_eur.R` | `MVMR_EUR改.R` | Multivariable MR, European ancestry |
| `44_mvmr_eas.R` | `MVMR_EAS改.R` | Multivariable MR, East Asian ancestry |
| `45_mvmr_sas.R` | `MVMR_SAS改.R` | Multivariable MR, South Asian ancestry |
| `46_mr_forest.R` | `MR_forest.R` | MR forest plot (Fig. 6) |
| `50_b5_integration.R` | `D5_full_pipeline.R` | Cross-stream integration (v8) |
| `51_b5_integration_v9.R` | `D5_full_pipeline_v9.R` | Cross-stream integration (final, v9) |
| `52_stability_check.R` | `stability_check.R` | Rule-application stability checks |
| `53_verify_bloc_adjacency.R` | `verify_bloc_adjacency.R` | Verification of bloc contiguity |
| `60_fig1a_bloc_selection.R` | `fig1a_bloc_selection.R` | Figure 1a |
| `61_fig1b_three_orderings.R` | `fig1b_three_orderings.R` | Figure 1b |
| `62_fig2b_map.R` | `Fig2B_map.R` | Figure 2b map |
| `63_fig_subtype_cluster_map.R` | `Subtype_ClusterMap.R` | Subtype-composition cluster map |

## Rule archive

`rules/` contains the prespecified decision rules used to classify each risk
factor. All four mechanisms (M1–M4) have equal standing; each is admitted only
where it predicts both the presence and the direction of a divergence and its
precondition holds.

- `rules/mechanisms.csv` — M1–M4 with the definition, the pattern each predicts
  and its precondition (Table S16). The M2 precondition lists its four
  conditions and states what it does not require, in the same form as the other
  mechanisms.
- `rules/classification.csv` — factor-by-factor application of the mechanisms,
  including the final class (C, SED, UD) (Table S17).
- `rules/README.md` — the three-stage decision procedure.

## Notes

- Figure scripts render TIFF at 600 dpi (300 dpi copies where provided) through
  `ragg`.
- The code contains no local paths, usernames or credentials. OpenGWAS tokens
  are read from the environment.
- Intermediate objects are written as CSV files so that each step can be
  inspected independently.
