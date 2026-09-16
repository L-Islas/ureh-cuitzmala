# Ecohydrological Response Units (UREH) in the Cuitzmala River Basin

Reproducible code, analysis-ready data, model specification, numerical outputs, spatial results, and manuscript figures supporting an ecohydrological regionalization of the Cuitzmala River Basin, Jalisco, Mexico.

The workflow develops **Ecohydrological Response Units (UREH)** by integrating local biophysical attributes, hydrological processes, and upstream water contribution within a Bayesian-network framework. The resulting regionalization differentiates hydrological function—**contribution, transfer, and reception**—and biophysical condition—**conserved and degraded**.

The repository also evaluates the sensitivity of the regionalization to controlled changes in the available evidence and examines its correspondence with an analytically independent regionalization based on multivariate similarity.

---

## Study system

The Cuitzmala River Basin is located on the southern coast of Jalisco, Mexico, extending from inland mountainous areas to the Pacific coast.

The analysis is conducted over **899 spatial units (microcatchments)** covering approximately **112,023 ha**.

The main analysis-ready spatial dataset is:

```text
data/input/ureh_units.gpkg
```

Coordinate reference system:

```text
EPSG:32613 — WGS 84 / UTM zone 13N
```

Each unit contains the spatial and hydrological attributes required by the analytical workflow, including land cover, soil, slope, precipitation, runoff, and upstream water contribution.

---

## Scope of reproducibility

This repository reproduces the analytical workflow **from the analysis-ready spatial units onward**.

It reproduces:

- preprocessing and discretization of analytical variables;
- construction of soft runoff evidence;
- parameterization of the Bayesian network;
- baseline inference (M0);
- posterior UREH probabilities and dominant classifications;
- posterior-support metrics;
- controlled sensitivity scenarios;
- Thornthwaite–Mather diagnostic analysis;
- multivariate clustering based on Gower distance and PAM;
- comparison of clustering configurations C1 and C2;
- descriptive comparison between multivariate clusters and UREH functional organization;
- final numerical tables;
- final spatial result layers;
- reproducibility regression checks.

The repository **does not regenerate the upstream geospatial preprocessing from the original raw remote-sensing and cartographic sources**.

In particular, it does not reproduce from raw data the original Landsat classification, DEM processing, drainage delineation, or the derivation of all continuous hydrological variables used to construct `ureh_units.gpkg`.

Those upstream procedures are documented in the associated manuscript. The repository begins from the consolidated analysis-ready spatial dataset used for the reported statistical and probabilistic analyses.

---

## Bayesian-network structure

The Bayesian network contains nine nodes:

```text
Soil
Land Cover
Slope
Precipitation
Upstream Water Contribution
Infiltration
Runoff
Hydrological Function
UREH
```

Its dependency structure is:

```text
Soil --------------------\
Land Cover ---------------\
Slope --------------------- > Infiltration ---> Runoff ----\
Precipitation ------------/        ^              ^          \
                            \       |              |           > Hydrological Function ---> UREH
                             \------|--------------/          /
Upstream Water Contribution -------------------------------/

Land Cover ------------------------------------------------> UREH
```

The operational model specification, state definitions, conditional probability tables, and provenance files are stored under:

```text
data/model/
```

---

## UREH classes

The final regionalization contains six categories resulting from the combination of hydrological function and biophysical condition:

| Hydrological function | Conserved | Degraded |
|---|---|---|
| Contribution | Contribution – Conserved | Contribution – Degraded |
| Transfer | Transfer – Conserved | Transfer – Degraded |
| Reception | Reception – Conserved | Reception – Degraded |

The dominant UREH is the maximum-a-posteriori (MAP) category for each spatial unit.

However, the complete posterior probability distribution is retained and reported. The dominant classification should therefore not be interpreted as a deterministic class assignment.

---

## Repository structure

```text
ureh-cuitzmala/
│
├── README.md
├── LICENSE
├── DATA_LICENSE.md
├── CITATION.cff
├── renv.lock
├── .Rprofile
├── .gitignore
├── ureh-cuitzmala.Rproj
├── run_all.R
│
├── data/
│   ├── README.md
│   │
│   ├── input/
│   │   ├── README.md
│   │   ├── ureh_units.gpkg
│   │   └── data_dictionary.csv
│   │
│   └── model/
│       ├── README.md
│       ├── model_states.csv
│       ├── cpt_*.csv
│       ├── source/
│       │   └── ureh_bayesian_network_model_specification.xlsx
│       └── provenance/
│
├── R/
│   ├── 00_config.R
│   ├── 01_prepare_data.R
│   ├── 02_soft_runoff_evidence.R
│   ├── 03_build_bayesian_network.R
│   ├── 04_infer_M0.R
│   ├── 05_summarize_M0.R
│   ├── 06_sensitivity_analysis.R
│   ├── 07_summarize_sensitivity.R
│   ├── 08_TM_diagnostic.R
│   ├── 09_prepare_clustering.R
│   ├── 10_gower_PAM.R
│   ├── 11_compare_C1_C2.R
│   ├── 12_compare_C2_UREH.R
│   └── 13_build_spatial_results.R
│
├── checks/
│   ├── expected_results.csv
│   ├── check_reproducibility.R
│   └── reproducibility_report.csv
│
├── results/
│   ├── tables/
│   ├── spatial/
│   │   └── ureh_results.gpkg
│   │
│   └── figures/
│       ├── README.md
│       ├── Fig01_study_area.png
│       ├── Fig02_bayesian_network.png
│       ├── Fig03_UREH_posterior_support.png
│       ├── Fig04_sensitivity_TM.png
│       └── Fig05_clustering_functional_composition.png
│
└── temp/
```

`temp/` contains intermediate computational products and is excluded from version control because its contents can be regenerated by the workflow.

---

## Analytical workflow

The analytical scripts are intended to be executed sequentially.

| Script | Purpose |
|---|---|
| `00_config.R` | Shared paths, parameters, state definitions, and workflow configuration |
| `01_prepare_data.R` | Validate input data, preprocess variables, and generate discrete analytical states |
| `02_soft_runoff_evidence.R` | Construct and regularize soft evidence for runoff |
| `03_build_bayesian_network.R` | Construct and parameterize the Bayesian network |
| `04_infer_M0.R` | Run baseline Bayesian inference and derive UREH posteriors |
| `05_summarize_M0.R` | Summarize UREH composition and posterior support |
| `06_sensitivity_analysis.R` | Run controlled evidence-removal sensitivity scenarios |
| `07_summarize_sensitivity.R` | Summarize MAP stability, Jensen–Shannon distance, and entropy changes |
| `08_TM_diagnostic.R` | Evaluate the effect of Thornthwaite–Mather runoff evidence |
| `09_prepare_clustering.R` | Prepare C1 and C2 multivariate clustering datasets |
| `10_gower_PAM.R` | Calculate Gower distances and evaluate PAM solutions |
| `11_compare_C1_C2.R` | Compare C1 and C2 clustering configurations |
| `12_compare_C2_UREH.R` | Describe correspondence between C2 clusters and UREH organization |
| `13_build_spatial_results.R` | Assemble the final spatial result layers |

The complete workflow is executed by:

```text
run_all.R
```

---

## Reproducing the analysis

Clone the repository:

```bash
git clone https://github.com/L-Islas/ureh-cuitzmala.git
cd ureh-cuitzmala
```

Open `ureh-cuitzmala.Rproj` in RStudio or start R from the repository root.

If the repository contains the finalized `renv.lock`, restore the package environment with:

```r
install.packages("renv")
renv::restore()
```

Then run the complete workflow:

```r
source("run_all.R")
```

`run_all.R` executes the analytical scripts sequentially and finishes by running the reproducibility checks.

The workflow is designed so that analytical steps communicate through explicitly generated files rather than depending on hidden objects remaining in the R session.

Each analytical script is executed in an isolated environment by the workflow controller.

---

## Reproducibility checks

The repository contains a frozen reference baseline:

```text
checks/expected_results.csv
```

and an automated regression-check script:

```text
checks/check_reproducibility.R
```

The checks compare regenerated outputs against **44 previously verified numerical and spatial controls**, covering:

- dataset dimensions and study-area extent;
- M0 UREH frequencies;
- posterior probability and uncertainty metrics;
- sensitivity results;
- Thornthwaite–Mather diagnostics;
- clustering results;
- spatial layer dimensions and coordinate reference system.

A successful complete run reports:

```text
Analytical scripts passed: 14 / 14
Reproducibility controls passed: 44 / 44
```

The 14 workflow steps consist of the 13 analytical scripts plus the final reproducibility-check script.

These checks assess **computational reproducibility and regression consistency** of the specified outputs.

They do **not** constitute independent or external validation of the Bayesian model or of the inferred UREH classes.

---

## Baseline model and sensitivity scenarios

### M0

M0 is the baseline inference using the complete evidence configuration defined for the study.

### M1

Removes local evidence for:

```text
Land Cover
Soil
Slope
```

while retaining precipitation, upstream contribution, and runoff evidence.

### M2

Removes:

```text
Precipitation
Runoff soft evidence
```

while retaining the remaining local structural evidence and upstream contribution.

### M3

Removes:

```text
Upstream Water Contribution
```

while retaining the remaining evidence.

Sensitivity is evaluated using complementary measures of posterior response, including:

- MAP-category stability;
- Jensen–Shannon distance between posterior distributions;
- change in normalized Shannon entropy.

The sensitivity analysis evaluates the response of the regionalization to controlled changes in evidence while keeping the model architecture, state definitions, conditional probability tables, priors, and inference structure fixed.

---

## Thornthwaite–Mather diagnostic

The Thornthwaite–Mather analysis is treated as a **diagnostic experiment**, not as an additional baseline model.

It compares network-derived runoff probabilities with runoff probabilities updated using the Thornthwaite–Mather evidence and evaluates how that update propagates toward hydrological function and UREH inference.

This diagnostic is reported separately from the principal M1–M3 sensitivity scenarios.

---

## Multivariate regionalization

An analytically independent regionalization based on multivariate similarity is used to examine the correspondence between structural similarity and inferred hydrological function.

Two configurations are considered:

### C1

Includes:

```text
Land Cover
Soil
Slope
Precipitation
Runoff
```

### C2

Includes the same variables plus:

```text
Upstream Water Contribution
```

Mixed-variable dissimilarity is calculated using **Gower distance**, followed by **Partitioning Around Medoids (PAM)**.

Candidate solutions from `k = 2` to `k = 20` are evaluated using silhouette statistics.

The selected solution contains:

```text
k = 12
```

for both C1 and C2.

C1 and C2 are highly similar, with only a small fraction of spatial units reassigned after upstream contribution is incorporated.

The subsequent comparison between C2 clusters and UREH is strictly descriptive.

It is intended to evaluate whether multivariate similarity and inferred hydrological function show spatial correspondence.

It should **not** be interpreted as classification accuracy, cluster purity, external validation, or evidence that one regionalization method is superior to the other.

---

## Main generated outputs

### Tables

Numerical results are written to:

```text
results/tables/
```

These include:

- M0 composition and posterior-support summaries;
- sensitivity summaries and unit-level diagnostics;
- Thornthwaite–Mather diagnostics;
- clustering evaluation and assignments;
- C1–C2 comparisons;
- C2–UREH functional and biophysical composition tables;
- spatial-layer manifest.

### Spatial results

The final generated spatial results are stored in:

```text
results/spatial/ureh_results.gpkg
```

with four layers:

```text
m0_ureh
sensitivity
tm_diagnostic
clustering_c1_c2
```

The GeoPackage is generated automatically by:

```text
R/13_build_spatial_results.R
```

and is therefore a reproducible output of the analytical workflow.

---

## Manuscript figures

Final manuscript figures are stored in:

```text
results/figures/
```

The five main figures are:

```text
Fig01_study_area.png
Fig02_bayesian_network.png
Fig03_UREH_posterior_support.png
Fig04_sensitivity_TM.png
Fig05_clustering_functional_composition.png
```

The numerical and spatial information underlying these figures is generated reproducibly by the analytical workflow.

However, the **final graphical composition is not reproduced pixel-for-pixel by `run_all.R`**.

Final figure preparation involved a combination of reproducible analytical outputs, GIS/cartographic processing, and graphical editing.

The repository therefore distinguishes between:

```text
reproducible analytical and spatial results
```

and:

```text
final publication-ready graphical composition
```

See:

```text
results/figures/README.md
```

for additional information.

---

## Data documentation

The analysis-ready input data are documented in:

```text
data/input/README.md
data/input/data_dictionary.csv
```

The Bayesian-network model files are documented in:

```text
data/model/README.md
```

Model provenance information is retained under:

```text
data/model/provenance/
```

These files document the operational model used by the reproducibility workflow and preserve traceability between the original model specification and the machine-readable conditional probability tables used in R.

---

## Software environment

The analysis is implemented in **R**.

Core packages used across the workflow include packages for:

```text
spatial data processing
Bayesian networks
data manipulation
clustering
statistical diagnostics
```

Exact package versions used for the archived reproducible environment are recorded in:

```text
renv.lock
```

The environment can be restored using:

```r
renv::restore()
```

---

## Citation

A permanent archival DOI will be added after release of the repository through Zenodo.

Until then, please cite the associated manuscript and this repository:

```text
L-Islas. Ecohydrological Response Units (UREH) in the Cuitzmala River Basin.
GitHub repository: https://github.com/L-Islas/ureh-cuitzmala
```

Machine-readable citation metadata are provided in:

```text
CITATION.cff
```

The citation information in this section and in `CITATION.cff` will be updated when the manuscript citation and Zenodo DOI are finalized.

---

## License

Source code in this repository is distributed under the terms specified in:

```text
LICENSE
```

Data, derived spatial products, and associated reuse conditions are described separately in:

```text
DATA_LICENSE.md
```

Data originating from external institutions or public cartographic sources remain subject to the terms and attribution requirements of their original providers.

---

## Repository status

The analytical workflow has been tested end-to-end from the analysis-ready input dataset through the final numerical and spatial products.

The verified workflow currently reproduces:

```text
14 / 14 workflow steps
44 / 44 reproducibility controls
```

The scientific analysis should be considered frozen for the archived manuscript version except where a documented error requires correction.

A tagged archival release (`v1.0.0`) and permanent Zenodo DOI will be created when the manuscript-associated repository version is formally frozen.