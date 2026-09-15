# Input data

This directory contains the analysis-ready spatial dataset used as the input for the UREH regionalization of the Cuitzmala River Basin.

## Main dataset

`ureh_units.gpkg`

- Layer: `ureh_units`
- Spatial units: 899 microcatchments
- Geometry: MULTIPOLYGON
- Coordinate reference system: WGS 84 / UTM zone 13N (EPSG:32613)
- Unique identifier: `id` (1–899)

The dataset was prepared from the spatial dataset used in the original analysis workflow. The original source identifier `ID_` was retained and renamed to `id`.

Only variables required to reproduce the Bayesian-network and multivariate analyses are retained.

## Variables

The GeoPackage contains:

- `id`: unique microcatchment identifier.
- `USyV`: land-use and land-cover category.
- `Suelos`: soil category.
- `Grad_Pend`: slope class.
- `P_sum`: precipitation variable.
- `Esc_sum`: local runoff variable.
- `Agua_D`: accumulated upstream water contribution.
- `geom`: microcatchment geometry.

Detailed variable descriptions are provided in `data_dictionary.csv`.

## Missing values

The analysis-ready input preserves missing values present in the source dataset:

- `P_sum`: 13 missing values.
- `Esc_sum`: 1 missing value.
- all other analytical variables: no missing values.

Missing values are not silently modified in this input file.

During data preparation, the single missing value in `Esc_sum` is explicitly converted to zero following the documented historical analysis workflow. Missing `P_sum` values are treated explicitly according to the requirements of each analysis.

## Derived variables

Variables such as discretized precipitation, runoff and upstream contribution classes, Bayesian-network states, `log1p(Agua_D)`, UREH classifications, posterior diagnostics, sensitivity metrics and clustering results are not stored in the input dataset.

They are generated reproducibly by the analysis scripts.

## Scope

Generation and validation of the land-cover product are outside the scope of this repository and will be documented separately.
