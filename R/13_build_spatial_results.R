# ============================================================
# 13. BUILD SPATIAL RESULTS
# ============================================================
#
# Purpose
# -------
# Assemble the reproducible analytical results with the
# authoritative microcatchment geometry.
#
# IMPORTANT
# ---------
# This script performs NO new statistical or probabilistic
# analysis.
#
# Geometry is always taken from:
#
#   data/input/ureh_units.gpkg
#
# Tabular analytical results are joined exclusively by `id`.
#
# Generated GeoPackage:
#
#   results/spatial/ureh_results.gpkg
#
# Layers:
#
#   m0_ureh
#   sensitivity
#   tm_diagnostic
#   clustering_c1_c2
#
# The source GeoPackage is NEVER modified.
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(sf)
library(dplyr)


results_spatial_dir <- here::here(
  "results",
  "spatial"
)

results_tables_dir <- here::here(
  "results",
  "tables"
)

dir.create(
  results_spatial_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

source_gpkg <- here::here(
  "data",
  "input",
  "ureh_units.gpkg"
)

output_gpkg <- here::here(
  "results",
  "spatial",
  "ureh_results.gpkg"
)


M0_file <- here::here(
  "temp",
  "M0_UREH_resultados_completos.csv"
)

sensitivity_file <- here::here(
  "results",
  "tables",
  "sensitivity_by_unit.csv"
)

TM_file <- here::here(
  "results",
  "tables",
  "tm_diagnostic_by_unit.csv"
)

clustering_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_assignments.csv"
)

C2_UREH_file <- here::here(
  "results",
  "tables",
  "clustering_C2_UREH_by_unit.csv"
)


required_files <- c(
  source_gpkg,
  M0_file,
  sensitivity_file,
  TM_file,
  clustering_file,
  C2_UREH_file
)


missing_files <- required_files[
  !file.exists(
    required_files
  )
]


if (
  length(
    missing_files
  ) > 0L
) {
  stop(
    "Required files not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 2. Validate source GeoPackage
# ------------------------------------------------------------

source_layers <- sf::st_layers(
  source_gpkg
)$name


if (!(
  "ureh_units" %in%
  source_layers
)) {
  stop(
    paste0(
      "Expected layer 'ureh_units' not found in:\n",
      source_gpkg,
      "\nAvailable layers: ",
      paste(
        source_layers,
        collapse = ", "
      )
    )
  )
}


units <- sf::st_read(
  source_gpkg,
  layer = "ureh_units",
  quiet = TRUE
)


# ------------------------------------------------------------
# 3. Validate authoritative geometry
# ------------------------------------------------------------

if (
  nrow(
    units
  ) !=
  n_expected_units
) {
  stop(
    "Authoritative geometry contains ",
    nrow(units),
    " units; expected ",
    n_expected_units,
    "."
  )
}


if (!(
  "id" %in%
  names(units)
)) {
  stop(
    "Authoritative geometry does not contain an id field."
  )
}


if (
  anyNA(
    units$id
  ) ||
  anyDuplicated(
    units$id
  )
) {
  stop(
    "Invalid ids in authoritative geometry."
  )
}


if (!identical(
  sort(
    as.integer(
      units$id
    )
  ),
  seq_len(
    n_expected_units
  )
)) {
  stop(
    "Authoritative geometry does not contain exactly ids 1:899."
  )
}


source_epsg <- sf::st_crs(
  units
)$epsg


if (
  is.na(
    source_epsg
  ) ||
  source_epsg !=
  expected_epsg
) {
  stop(
    "Unexpected CRS. Expected EPSG:",
    expected_epsg,
    "; found EPSG:",
    source_epsg,
    "."
  )
}


if (
  any(
    sf::st_is_empty(
      units
    )
  )
) {
  stop(
    "Authoritative geometry contains empty geometries."
  )
}


geometry_validity <- sf::st_is_valid(
  units
)


if (
  anyNA(
    geometry_validity
  ) ||
  any(
    !geometry_validity
  )
) {
  stop(
    paste0(
      "Authoritative geometry contains invalid features. ",
      "The script will not repair them silently."
    )
  )
}


# ------------------------------------------------------------
# 4. Calculate area from authoritative geometry
# ------------------------------------------------------------
#
# EPSG 32613 is projected in metres.
#
# ------------------------------------------------------------

units$area_ha <- as.numeric(
  sf::st_area(
    units
  )
) /
  10000


if (
  anyNA(
    units$area_ha
  ) ||
  any(
    !is.finite(
      units$area_ha
    )
  ) ||
  any(
    units$area_ha <= 0
  )
) {
  stop(
    "Invalid microcatchment areas calculated from geometry."
  )
}


# ------------------------------------------------------------
# 5. Helper: validate unit-level tabular data
# ------------------------------------------------------------

validate_unit_table <- function(
    x,
    object_name,
    require_exact_899 = TRUE
) {
  
  if (!(
    "id" %in%
    names(x)
  )) {
    stop(
      object_name,
      " does not contain an id field."
    )
  }
  
  
  if (
    anyNA(
      x$id
    )
  ) {
    stop(
      object_name,
      " contains missing ids."
    )
  }
  
  
  if (
    require_exact_899
  ) {
    
    if (
      nrow(x) !=
      n_expected_units
    ) {
      stop(
        object_name,
        " must contain ",
        n_expected_units,
        " rows."
      )
    }
    
    
    if (
      anyDuplicated(
        x$id
      )
    ) {
      stop(
        object_name,
        " contains duplicated ids."
      )
    }
    
    
    if (!identical(
      sort(
        as.integer(
          x$id
        )
      ),
      seq_len(
        n_expected_units
      )
    )) {
      stop(
        object_name,
        " does not contain exactly ids 1:899."
      )
    }
  }
  
  
  invisible(
    TRUE
  )
}


# ------------------------------------------------------------
# 6. Helper: join only new columns
# ------------------------------------------------------------
#
# Source geometry attributes remain authoritative.
#
# If a tabular output contains columns already present in the
# spatial source, those duplicated fields are not added again.
#
# ------------------------------------------------------------

join_new_columns <- function(
    spatial_data,
    tabular_data
) {
  
  columns_to_add <- setdiff(
    names(
      tabular_data
    ),
    names(
      spatial_data
    )
  )
  
  
  columns_to_add <- setdiff(
    columns_to_add,
    attr(
      spatial_data,
      "sf_column"
    )
  )
  
  
  tabular_subset <- tabular_data |>
    select(
      id,
      all_of(
        columns_to_add
      )
    )
  
  
  spatial_data |>
    left_join(
      tabular_subset,
      by = "id"
    )
}


# ------------------------------------------------------------
# 7. Read M0 results
# ------------------------------------------------------------

M0 <- read.csv(
  M0_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


validate_unit_table(
  M0,
  "M0 results"
)


if (!(
  "UREH_predicha" %in%
  names(M0)
)) {
  stop(
    "M0 results do not contain UREH_predicha."
  )
}


# ------------------------------------------------------------
# 8. Build M0 spatial layer
# ------------------------------------------------------------

m0_spatial <- join_new_columns(
  units,
  M0
)


m0_spatial <- m0_spatial |>
  mutate(
    
    funcion_UREH =
      case_when(
        
        grepl(
          "^aporte_",
          UREH_predicha
        ) ~ "aporte",
        
        grepl(
          "^transferencia_",
          UREH_predicha
        ) ~ "transferencia",
        
        grepl(
          "^receptora_",
          UREH_predicha
        ) ~ "recepcion",
        
        TRUE ~ NA_character_
      ),
    
    
    condicion_UREH =
      case_when(
        
        grepl(
          "_conservada$",
          UREH_predicha
        ) ~ "conservada",
        
        grepl(
          "_degradada$",
          UREH_predicha
        ) ~ "degradada",
        
        TRUE ~ NA_character_
      )
  )


if (
  anyNA(
    m0_spatial$funcion_UREH
  ) ||
  anyNA(
    m0_spatial$condicion_UREH
  )
) {
  stop(
    "Could not derive M0 function or condition from UREH labels."
  )
}


# ------------------------------------------------------------
# 9. Read sensitivity results
# ------------------------------------------------------------

sensitivity <- read.csv(
  sensitivity_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_sensitivity_fields <- c(
  "id",
  "escenario",
  "UREH_M0",
  "UREH_escenario",
  "MAP_estable",
  "distancia_JS",
  "delta_entropia"
)


missing_sensitivity_fields <- setdiff(
  required_sensitivity_fields,
  names(
    sensitivity
  )
)


if (
  length(
    missing_sensitivity_fields
  ) > 0L
) {
  stop(
    "Sensitivity table is missing fields: ",
    paste(
      missing_sensitivity_fields,
      collapse = ", "
    )
  )
}


expected_sensitivity_scenarios <- c(
  "M1",
  "M2",
  "M3",
  "Sin cobertura",
  "Sin suelo + pendiente"
)


if (!setequal(
  unique(
    sensitivity$escenario
  ),
  expected_sensitivity_scenarios
)) {
  stop(
    "Unexpected sensitivity scenarios."
  )
}


sensitivity_counts <- sensitivity |>
  count(
    escenario,
    name = "n"
  )


if (
  any(
    sensitivity_counts$n !=
    n_expected_units
  )
) {
  stop(
    "Each sensitivity scenario must contain ",
    n_expected_units,
    " units."
  )
}


if (
  anyDuplicated(
    sensitivity[
      ,
      c(
        "id",
        "escenario"
      )
    ]
  )
) {
  stop(
    "Duplicated id-scenario combinations in sensitivity results."
  )
}


# ------------------------------------------------------------
# 10. Build sensitivity spatial layer
# ------------------------------------------------------------
#
# Long format is intentional:
#
# one geometry is repeated for each sensitivity scenario.
#
# This makes scenario filtering and mapping straightforward.
#
# ------------------------------------------------------------

sensitivity_spatial <- units |>
  inner_join(
    sensitivity,
    by = "id"
  ) |>
  arrange(
    escenario,
    id
  )


expected_sensitivity_features <-
  n_expected_units *
  length(
    expected_sensitivity_scenarios
  )


if (
  nrow(
    sensitivity_spatial
  ) !=
  expected_sensitivity_features
) {
  stop(
    "Unexpected number of spatial sensitivity features."
  )
}


# ------------------------------------------------------------
# 11. Read TM diagnostic
# ------------------------------------------------------------

TM <- read.csv(
  TM_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


validate_unit_table(
  TM,
  "TM diagnostic"
)


required_TM_fields <- c(
  "id",
  "MAP_runoff_network",
  "MAP_runoff_TM",
  "JS_runoff",
  "delta_H_runoff",
  "MAP_UREH_network",
  "MAP_UREH_TM",
  "JS_UREH",
  "delta_H_UREH"
)


missing_TM_fields <- setdiff(
  required_TM_fields,
  names(
    TM
  )
)


if (
  length(
    missing_TM_fields
  ) > 0L
) {
  stop(
    "TM diagnostic is missing fields: ",
    paste(
      missing_TM_fields,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 12. Build TM spatial layer
# ------------------------------------------------------------

TM_spatial <- join_new_columns(
  units,
  TM
)


if (
  nrow(
    TM_spatial
  ) !=
  n_expected_units
) {
  stop(
    "Unexpected number of TM spatial features."
  )
}


# ------------------------------------------------------------
# 13. Read clustering assignments
# ------------------------------------------------------------

clustering <- read.csv(
  clustering_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


validate_unit_table(
  clustering,
  "C1/C2 clustering assignments"
)


required_clustering_fields <- c(
  "id",
  "cluster_C1",
  "cluster_C2_original",
  "cluster_C2",
  "changed_cluster",
  "silhouette_C1",
  "silhouette_C2"
)


missing_clustering_fields <- setdiff(
  required_clustering_fields,
  names(
    clustering
  )
)


if (
  length(
    missing_clustering_fields
  ) > 0L
) {
  stop(
    "Clustering assignments are missing fields: ",
    paste(
      missing_clustering_fields,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 14. Read C2 x UREH unit-level table
# ------------------------------------------------------------

C2_UREH <- read.csv(
  C2_UREH_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


validate_unit_table(
  C2_UREH,
  "C2-UREH unit table"
)


required_C2_UREH_fields <- c(
  "id",
  "cluster_C2",
  "UREH_predicha",
  "funcion",
  "condicion"
)


missing_C2_UREH_fields <- setdiff(
  required_C2_UREH_fields,
  names(
    C2_UREH
  )
)


if (
  length(
    missing_C2_UREH_fields
  ) > 0L
) {
  stop(
    "C2-UREH table is missing fields: ",
    paste(
      missing_C2_UREH_fields,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 15. Cross-check cluster_C2 between both outputs
# ------------------------------------------------------------

cluster_match <- clustering |>
  select(
    id,
    cluster_C2
  ) |>
  left_join(
    C2_UREH |>
      select(
        id,
        cluster_C2_check =
          cluster_C2
      ),
    by = "id"
  )


if (
  anyNA(
    cluster_match$cluster_C2_check
  ) ||
  any(
    cluster_match$cluster_C2 !=
    cluster_match$cluster_C2_check
  )
) {
  stop(
    "cluster_C2 differs between clustering and C2-UREH outputs."
  )
}


# ------------------------------------------------------------
# 16. Build clustering spatial layer
# ------------------------------------------------------------

clustering_spatial <- join_new_columns(
  units,
  clustering
)


C2_semantic_fields <- C2_UREH |>
  select(
    id,
    UREH_predicha,
    funcion,
    condicion
  )


clustering_spatial <- join_new_columns(
  clustering_spatial,
  C2_semantic_fields
)


if (
  nrow(
    clustering_spatial
  ) !=
  n_expected_units
) {
  stop(
    "Unexpected number of clustering spatial features."
  )
}


# ------------------------------------------------------------
# 17. Validate geometry consistency across layers
# ------------------------------------------------------------

validate_spatial_layer <- function(
    x,
    layer_name,
    expected_rows,
    expected_unique_ids = n_expected_units
) {
  
  if (
    nrow(x) !=
    expected_rows
  ) {
    stop(
      layer_name,
      ": unexpected feature count."
    )
  }
  
  
  if (
    length(
      unique(
        x$id
      )
    ) !=
    expected_unique_ids
  ) {
    stop(
      layer_name,
      ": unexpected number of unique ids."
    )
  }
  
  
  if (
    sf::st_crs(x)$epsg !=
    expected_epsg
  ) {
    stop(
      layer_name,
      ": unexpected CRS."
    )
  }
  
  
  if (
    any(
      sf::st_is_empty(
        x
      )
    )
  ) {
    stop(
      layer_name,
      ": empty geometries detected."
    )
  }
  
  
  if (
    any(
      !sf::st_is_valid(
        x
      )
    )
  ) {
    stop(
      layer_name,
      ": invalid geometries detected."
    )
  }
  
  
  invisible(
    TRUE
  )
}


validate_spatial_layer(
  m0_spatial,
  "m0_ureh",
  n_expected_units
)


validate_spatial_layer(
  sensitivity_spatial,
  "sensitivity",
  expected_sensitivity_features
)


validate_spatial_layer(
  TM_spatial,
  "tm_diagnostic",
  n_expected_units
)


validate_spatial_layer(
  clustering_spatial,
  "clustering_c1_c2",
  n_expected_units
)


# ------------------------------------------------------------
# 18. Write GeoPackage layers
# ------------------------------------------------------------
#
# delete_layer = TRUE replaces ONLY the named layer.
#
# The complete GeoPackage is never deleted with unlink().
#
# ------------------------------------------------------------

write_gpkg_layer <- function(
    x,
    layer_name
) {
  
  sf::st_write(
    obj =
      x,
    
    dsn =
      output_gpkg,
    
    layer =
      layer_name,
    
    driver =
      "GPKG",
    
    delete_layer =
      TRUE,
    
    quiet =
      TRUE
  )
}


message("")
message(
  "============================================"
)
message(
  "WRITING SPATIAL RESULTS"
)
message(
  "============================================"
)


message(
  "Writing layer: m0_ureh"
)

write_gpkg_layer(
  m0_spatial,
  "m0_ureh"
)


message(
  "Writing layer: sensitivity"
)

write_gpkg_layer(
  sensitivity_spatial,
  "sensitivity"
)


message(
  "Writing layer: tm_diagnostic"
)

write_gpkg_layer(
  TM_spatial,
  "tm_diagnostic"
)


message(
  "Writing layer: clustering_c1_c2"
)

write_gpkg_layer(
  clustering_spatial,
  "clustering_c1_c2"
)


# ------------------------------------------------------------
# 19. Validate GeoPackage layer inventory
# ------------------------------------------------------------

expected_output_layers <- c(
  "m0_ureh",
  "sensitivity",
  "tm_diagnostic",
  "clustering_c1_c2"
)


actual_output_layers <- sf::st_layers(
  output_gpkg
)$name


missing_output_layers <- setdiff(
  expected_output_layers,
  actual_output_layers
)


if (
  length(
    missing_output_layers
  ) > 0L
) {
  stop(
    "Generated GeoPackage is missing layers: ",
    paste(
      missing_output_layers,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 20. Read-back validation
# ------------------------------------------------------------

m0_check <- sf::st_read(
  output_gpkg,
  layer = "m0_ureh",
  quiet = TRUE
)

sensitivity_check <- sf::st_read(
  output_gpkg,
  layer = "sensitivity",
  quiet = TRUE
)

TM_check <- sf::st_read(
  output_gpkg,
  layer = "tm_diagnostic",
  quiet = TRUE
)

clustering_check <- sf::st_read(
  output_gpkg,
  layer = "clustering_c1_c2",
  quiet = TRUE
)


validate_spatial_layer(
  m0_check,
  "m0_ureh read-back",
  n_expected_units
)


validate_spatial_layer(
  sensitivity_check,
  "sensitivity read-back",
  expected_sensitivity_features
)


validate_spatial_layer(
  TM_check,
  "tm_diagnostic read-back",
  n_expected_units
)


validate_spatial_layer(
  clustering_check,
  "clustering_c1_c2 read-back",
  n_expected_units
)


# ------------------------------------------------------------
# 21. Strong M0 read-back control
# ------------------------------------------------------------

M0_original_order <- M0 |>
  arrange(
    id
  )


m0_check_order <- m0_check |>
  sf::st_drop_geometry() |>
  arrange(
    id
  )


if (!identical(
  as.character(
    m0_check_order$UREH_predicha
  ),
  as.character(
    M0_original_order$UREH_predicha
  )
)) {
  stop(
    "M0 UREH classification changed during spatial assembly."
  )
}


# ------------------------------------------------------------
# 22. Strong sensitivity read-back control
# ------------------------------------------------------------

sensitivity_original_control <- sensitivity |>
  arrange(
    escenario,
    id
  ) |>
  select(
    id,
    escenario,
    UREH_escenario,
    MAP_estable,
    distancia_JS,
    delta_entropia
  )


sensitivity_readback_control <-
  sensitivity_check |>
  sf::st_drop_geometry() |>
  arrange(
    escenario,
    id
  ) |>
  select(
    id,
    escenario,
    UREH_escenario,
    MAP_estable,
    distancia_JS,
    delta_entropia
  )


if (!identical(
  sensitivity_original_control$id,
  sensitivity_readback_control$id
)) {
  stop(
    "Sensitivity ids changed during spatial assembly."
  )
}


if (!identical(
  as.character(
    sensitivity_original_control$UREH_escenario
  ),
  as.character(
    sensitivity_readback_control$UREH_escenario
  )
)) {
  stop(
    "Sensitivity classifications changed during spatial assembly."
  )
}


numeric_sensitivity_difference <- max(
  abs(
    as.matrix(
      sensitivity_original_control[
        ,
        c(
          "distancia_JS",
          "delta_entropia"
        )
      ]
    ) -
      as.matrix(
        sensitivity_readback_control[
          ,
          c(
            "distancia_JS",
            "delta_entropia"
          )
        ]
      )
  )
)


if (
  numeric_sensitivity_difference >
  1e-12
) {
  stop(
    "Sensitivity metrics changed during GeoPackage write/read."
  )
}


# ------------------------------------------------------------
# 23. Spatial layer manifest
# ------------------------------------------------------------

spatial_layer_manifest <- data.frame(
  
  geopackage =
    rep(
      "results/spatial/ureh_results.gpkg",
      4
    ),
  
  layer = c(
    "m0_ureh",
    "sensitivity",
    "tm_diagnostic",
    "clustering_c1_c2"
  ),
  
  n_features = c(
    nrow(
      m0_check
    ),
    nrow(
      sensitivity_check
    ),
    nrow(
      TM_check
    ),
    nrow(
      clustering_check
    )
  ),
  
  n_unique_ids = c(
    length(
      unique(
        m0_check$id
      )
    ),
    length(
      unique(
        sensitivity_check$id
      )
    ),
    length(
      unique(
        TM_check$id
      )
    ),
    length(
      unique(
        clustering_check$id
      )
    )
  ),
  
  epsg = rep(
    expected_epsg,
    4
  ),
  
  purpose = c(
    "M0 UREH classification and posterior results",
    "Evidence-removal sensitivity scenarios in long format",
    "Network-only versus Thornthwaite-Mather diagnostic",
    "Final k=12 C1/C2 clustering assignments and C2-UREH attributes"
  ),
  
  stringsAsFactors = FALSE
)


manifest_file <- here::here(
  "results",
  "tables",
  "spatial_layer_manifest.csv"
)


write.csv(
  spatial_layer_manifest,
  manifest_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


manifest_check <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (
  nrow(
    manifest_check
  ) != 4L
) {
  stop(
    "Spatial layer manifest failed read-back validation."
  )
}


# ------------------------------------------------------------
# 24. Console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SPATIAL RESULTS COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "Source geometry:\n",
  "data/input/ureh_units.gpkg | layer: ureh_units\n\n",
  sep = ""
)


cat(
  "Generated GeoPackage:\n",
  "results/spatial/ureh_results.gpkg\n\n",
  sep = ""
)


print(
  spatial_layer_manifest,
  row.names = FALSE
)


cat(
  "\nM0 total area (ha): ",
  format(
    sum(
      m0_check$area_ha
    ),
    digits = 10
  ),
  "\n",
  sep = ""
)


cat(
  "Sensitivity features: ",
  nrow(
    sensitivity_check
  ),
  " = ",
  n_expected_units,
  " units x ",
  length(
    expected_sensitivity_scenarios
  ),
  " scenarios\n",
  sep = ""
)


cat(
  "Maximum sensitivity metric write/read difference: ",
  format(
    numeric_sensitivity_difference,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "CRS: EPSG:",
  expected_epsg,
  "\n"
)


cat(
  "Geometry source was not modified.\n"
)


cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 25. Completion
# ------------------------------------------------------------

message(
  "13_build_spatial_results.R completed successfully."
)

message(
  "Generated: results/spatial/ureh_results.gpkg"
)

message(
  "Layers: m0_ureh | sensitivity | tm_diagnostic | clustering_c1_c2"
)

message(
  "Generated: results/tables/spatial_layer_manifest.csv"
)

message(
  "Next: R/14_make_figures.R"
)