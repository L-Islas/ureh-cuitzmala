# ============================================================
# 05. SUMMARIZE M0
# ============================================================
#
# Purpose
# -------
# Summarize the final M0 UREH regionalization in terms of:
#
# 1. the six UREH classes;
# 2. the three hydrological functions;
# 3. the two relative biophysical conditions;
# 4. posterior-support diagnostics.
#
# Because microcatchments differ in area, both unit-based and
# area-based proportions are reported.
#
# Inputs
# ------
# data/input/ureh_units.gpkg
# results/tables/M0_UREH_clasificacion_dominante.csv
#
# Outputs
# -------
# results/tables/M0_UREH_composition.csv
# results/tables/M0_function_composition.csv
# results/tables/M0_condition_composition.csv
# results/tables/M0_posterior_support_summary.csv
# results/tables/M0_posterior_support_by_UREH.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(sf)
library(dplyr)

results_tables_dir <- here::here(
  "results",
  "tables"
)

dir.create(
  results_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Input paths
# ------------------------------------------------------------

classification_file <- here::here(
  "results",
  "tables",
  "M0_UREH_clasificacion_dominante.csv"
)

required_files <- c(
  file_ureh_units,
  classification_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required files not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ------------------------------------------------------------
# 2. Read spatial units
# ------------------------------------------------------------

ureh_sf <- st_read(
  file_ureh_units,
  layer = "ureh_units",
  quiet = TRUE
)

if (!inherits(
  ureh_sf,
  "sf"
)) {
  stop(
    "ureh_units.gpkg was not read as an sf object."
  )
}

if (nrow(ureh_sf) != n_expected_units) {
  stop(
    "Expected ",
    n_expected_units,
    " spatial units, found ",
    nrow(ureh_sf),
    "."
  )
}

if (
  is.na(st_crs(ureh_sf)$epsg) ||
  st_crs(ureh_sf)$epsg != expected_epsg
) {
  stop(
    "Unexpected CRS in ureH spatial input."
  )
}

if (anyNA(ureh_sf$id)) {
  stop(
    "Spatial input contains missing ids."
  )
}

if (anyDuplicated(ureh_sf$id)) {
  stop(
    "Spatial input contains duplicated ids."
  )
}


# ------------------------------------------------------------
# 3. Calculate microcatchment area
# ------------------------------------------------------------
#
# EPSG:32613 is projected in metres.
# sf::st_area() therefore returns square metres.
#
# 1 ha = 10,000 m2.
#
# ------------------------------------------------------------

ureh_sf$area_ha <-
  as.numeric(
    st_area(
      ureh_sf
    )
  ) /
  10000

if (
  anyNA(ureh_sf$area_ha) ||
  any(!is.finite(ureh_sf$area_ha)) ||
  any(ureh_sf$area_ha <= 0)
) {
  stop(
    "Invalid microcatchment areas detected."
  )
}

total_area_ha <- sum(
  ureh_sf$area_ha
)

if (
  !is.finite(total_area_ha) ||
  total_area_ha <= 0
) {
  stop(
    "Invalid total basin area."
  )
}


# ------------------------------------------------------------
# 4. Read M0 classification
# ------------------------------------------------------------

classification <- read.csv(
  classification_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_classification_fields <- c(
  "id",
  "UREH_predicha",
  "prob_UREH_predicha",
  "UREH_segunda",
  "prob_UREH_segunda",
  "margen_clasificacion",
  "entropia_UREH"
)

missing_fields <- setdiff(
  required_classification_fields,
  names(classification)
)

if (length(missing_fields) > 0L) {
  stop(
    "Missing fields in M0 classification: ",
    paste(
      missing_fields,
      collapse = ", "
    )
  )
}

if (
  nrow(classification) !=
  n_expected_units
) {
  stop(
    "M0 classification does not contain ",
    n_expected_units,
    " units."
  )
}

if (
  anyNA(classification$id) ||
  anyDuplicated(classification$id)
) {
  stop(
    "Invalid ids in M0 classification."
  )
}

if (!setequal(
  classification$id,
  ureh_sf$id
)) {
  stop(
    "Spatial input and M0 classification contain different ids."
  )
}


# ------------------------------------------------------------
# 5. Validate UREH states
# ------------------------------------------------------------

if (!all(
  classification$UREH_predicha %in%
  estados$ureh
)) {
  stop(
    "Unrecognized M0 UREH classes detected."
  )
}

if (!all(
  classification$UREH_segunda %in%
  estados$ureh
)) {
  stop(
    "Unrecognized second-ranked UREH classes detected."
  )
}


# ------------------------------------------------------------
# 6. Join area to M0 classification
# ------------------------------------------------------------

area_table <-
  st_drop_geometry(
    ureh_sf
  ) |>
  select(
    id,
    area_ha
  )

M0_summary_data <-
  classification |>
  left_join(
    area_table,
    by = "id"
  )

if (
  nrow(M0_summary_data) !=
  n_expected_units
) {
  stop(
    "Area join changed the number of M0 units."
  )
}

if (anyNA(
  M0_summary_data$area_ha
)) {
  stop(
    "Missing area values after joining M0 classification."
  )
}


# ------------------------------------------------------------
# 7. Derive hydrological function and condition
# ------------------------------------------------------------
#
# The six UREH states combine:
#
# function:
#   aporte / transferencia / recepcion
#
# condition:
#   conservada / degradada
#
# ------------------------------------------------------------

M0_summary_data <-
  M0_summary_data |>
  mutate(
    
    funcion = case_when(
      
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
    
    condicion = case_when(
      
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
  anyNA(M0_summary_data$funcion) ||
  anyNA(M0_summary_data$condicion)
) {
  stop(
    "Could not derive function or condition from one or more UREH states."
  )
}


# ------------------------------------------------------------
# 8. Generic composition summary
# ------------------------------------------------------------

summarize_composition <- function(
    data,
    grouping_variable
) {
  
  group_symbol <-
    rlang::sym(
      grouping_variable
    )
  
  result <-
    data |>
    group_by(
      !!group_symbol
    ) |>
    summarise(
      
      n_units = n(),
      
      area_ha =
        sum(
          area_ha
        ),
      
      .groups = "drop"
    ) |>
    mutate(
      
      unit_percentage =
        100 *
        n_units /
        sum(n_units),
      
      area_percentage =
        100 *
        area_ha /
        sum(area_ha)
    )
  
  result
}


# ------------------------------------------------------------
# 9. Six-class UREH composition
# ------------------------------------------------------------

M0_UREH_composition <-
  summarize_composition(
    M0_summary_data,
    "UREH_predicha"
  ) |>
  mutate(
    UREH_predicha = factor(
      UREH_predicha,
      levels =
        estados$ureh
    )
  ) |>
  arrange(
    UREH_predicha
  ) |>
  mutate(
    UREH_predicha =
      as.character(
        UREH_predicha
      )
  )

if (
  sum(
    M0_UREH_composition$n_units
  ) != n_expected_units
) {
  stop(
    "UREH composition does not sum to 899 units."
  )
}

if (
  abs(
    sum(
      M0_UREH_composition$area_percentage
    ) - 100
  ) >
  1e-10
) {
  stop(
    "UREH area percentages do not sum to 100."
  )
}


# ------------------------------------------------------------
# 10. Hydrological-function composition
# ------------------------------------------------------------

function_order <- c(
  "aporte",
  "transferencia",
  "recepcion"
)

M0_function_composition <-
  summarize_composition(
    M0_summary_data,
    "funcion"
  ) |>
  mutate(
    funcion = factor(
      funcion,
      levels =
        function_order
    )
  ) |>
  arrange(
    funcion
  ) |>
  mutate(
    funcion =
      as.character(
        funcion
      )
  )

if (
  sum(
    M0_function_composition$n_units
  ) != n_expected_units
) {
  stop(
    "Function composition does not sum to 899 units."
  )
}


# ------------------------------------------------------------
# 11. Relative biophysical-condition composition
# ------------------------------------------------------------

condition_order <- c(
  "conservada",
  "degradada"
)

M0_condition_composition <-
  summarize_composition(
    M0_summary_data,
    "condicion"
  ) |>
  mutate(
    condicion = factor(
      condicion,
      levels =
        condition_order
    )
  ) |>
  arrange(
    condicion
  ) |>
  mutate(
    condicion =
      as.character(
        condicion
      )
  )

if (
  sum(
    M0_condition_composition$n_units
  ) != n_expected_units
) {
  stop(
    "Condition composition does not sum to 899 units."
  )
}


# ------------------------------------------------------------
# 12. Posterior-support metrics
# ------------------------------------------------------------

support_metrics <- c(
  dominant_probability =
    "prob_UREH_predicha",
  
  classification_margin =
    "margen_clasificacion",
  
  normalized_entropy =
    "entropia_UREH"
)


# ------------------------------------------------------------
# 13. Overall posterior-support summary
# ------------------------------------------------------------

support_summary_list <- lapply(
  names(support_metrics),
  function(metric_name) {
    
    source_column <- support_metrics[[metric_name]]
    
    x <- M0_summary_data[[source_column]]
    
    if (
      anyNA(x) ||
      any(!is.finite(x))
    ) {
      stop(
        "Invalid values in posterior-support metric: ",
        metric_name,
        "."
      )
    }
    
    data.frame(
      
      metric =
        metric_name,
      
      n =
        length(x),
      
      mean =
        mean(x),
      
      sd =
        sd(x),
      
      min =
        min(x),
      
      q25 =
        unname(
          quantile(
            x,
            0.25
          )
        ),
      
      median =
        median(x),
      
      q75 =
        unname(
          quantile(
            x,
            0.75
          )
        ),
      
      max =
        max(x),
      
      stringsAsFactors = FALSE
    )
  }
)

M0_posterior_support_summary <-
  bind_rows(
    support_summary_list
  )


# ------------------------------------------------------------
# 14. Posterior support by MAP UREH
# ------------------------------------------------------------

M0_posterior_support_by_UREH <-
  M0_summary_data |>
  group_by(
    UREH_predicha
  ) |>
  summarise(
    
    n_units = n(),
    
    area_ha =
      sum(area_ha),
    
    dominant_probability_mean =
      mean(
        prob_UREH_predicha
      ),
    
    dominant_probability_median =
      median(
        prob_UREH_predicha
      ),
    
    classification_margin_mean =
      mean(
        margen_clasificacion
      ),
    
    classification_margin_median =
      median(
        margen_clasificacion
      ),
    
    normalized_entropy_mean =
      mean(
        entropia_UREH
      ),
    
    normalized_entropy_median =
      median(
        entropia_UREH
      ),
    
    .groups = "drop"
  ) |>
  mutate(
    
    area_percentage =
      100 *
      area_ha /
      total_area_ha,
    
    UREH_predicha = factor(
      UREH_predicha,
      levels =
        estados$ureh
    )
  ) |>
  arrange(
    UREH_predicha
  ) |>
  mutate(
    UREH_predicha =
      as.character(
        UREH_predicha
      )
  )


# ------------------------------------------------------------
# 15. Validate posterior-support ranges
# ------------------------------------------------------------

if (
  any(
    M0_summary_data$prob_UREH_predicha <
    -1e-10
  ) ||
  any(
    M0_summary_data$prob_UREH_predicha >
    1 + 1e-10
  )
) {
  stop(
    "Dominant posterior probabilities outside [0,1]."
  )
}

if (
  any(
    M0_summary_data$margen_clasificacion <
    -1e-10
  ) ||
  any(
    M0_summary_data$margen_clasificacion >
    1 + 1e-10
  )
) {
  stop(
    "Classification margins outside [0,1]."
  )
}

if (
  any(
    M0_summary_data$entropia_UREH <
    -1e-10
  ) ||
  any(
    M0_summary_data$entropia_UREH >
    1 + 1e-10
  )
) {
  stop(
    "Normalized entropy outside [0,1]."
  )
}


# ------------------------------------------------------------
# 16. Output paths
# ------------------------------------------------------------

ureh_composition_file <- here::here(
  "results",
  "tables",
  "M0_UREH_composition.csv"
)

function_composition_file <- here::here(
  "results",
  "tables",
  "M0_function_composition.csv"
)

condition_composition_file <- here::here(
  "results",
  "tables",
  "M0_condition_composition.csv"
)

support_summary_file <- here::here(
  "results",
  "tables",
  "M0_posterior_support_summary.csv"
)

support_by_ureh_file <- here::here(
  "results",
  "tables",
  "M0_posterior_support_by_UREH.csv"
)


# ------------------------------------------------------------
# 17. Write outputs
# ------------------------------------------------------------

write.csv(
  M0_UREH_composition,
  ureh_composition_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  M0_function_composition,
  function_composition_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  M0_condition_composition,
  condition_composition_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  M0_posterior_support_summary,
  support_summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  M0_posterior_support_by_UREH,
  support_by_ureh_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 18. Read-back validation
# ------------------------------------------------------------

ureh_composition_check <- read.csv(
  ureh_composition_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

function_composition_check <- read.csv(
  function_composition_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

condition_composition_check <- read.csv(
  condition_composition_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

support_summary_check <- read.csv(
  support_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

support_by_ureh_check <- read.csv(
  support_by_ureh_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(
    ureh_composition_check
  ) != 6L
) {
  stop(
    "Unexpected number of UREH classes after read-back."
  )
}

if (
  nrow(
    function_composition_check
  ) != 3L
) {
  stop(
    "Unexpected number of hydrological functions after read-back."
  )
}

if (
  nrow(
    condition_composition_check
  ) != 2L
) {
  stop(
    "Unexpected number of conditions after read-back."
  )
}

if (
  nrow(
    support_summary_check
  ) != 3L
) {
  stop(
    "Unexpected number of posterior-support metrics after read-back."
  )
}

if (
  nrow(
    support_by_ureh_check
  ) != 6L
) {
  stop(
    "Unexpected number of UREH support summaries after read-back."
  )
}


# ------------------------------------------------------------
# 19. Console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "M0 SUMMARY\n"
)

cat(
  "============================================\n"
)

cat(
  "Total basin area (ha):",
  total_area_ha,
  "\n\n"
)

cat(
  "Six UREH classes:\n"
)

print(
  M0_UREH_composition,
  row.names = FALSE
)

cat(
  "\nHydrological functions:\n"
)

print(
  M0_function_composition,
  row.names = FALSE
)

cat(
  "\nRelative biophysical condition:\n"
)

print(
  M0_condition_composition,
  row.names = FALSE
)

cat(
  "\nPosterior-support diagnostics:\n"
)

print(
  M0_posterior_support_summary,
  row.names = FALSE
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 20. Completion
# ------------------------------------------------------------

message(
  "05_summarize_M0.R completed successfully."
)

message(
  "Generated: results/tables/M0_UREH_composition.csv"
)

message(
  "Generated: results/tables/M0_function_composition.csv"
)

message(
  "Generated: results/tables/M0_condition_composition.csv"
)

message(
  "Generated: results/tables/M0_posterior_support_summary.csv"
)

message(
  "Generated: results/tables/M0_posterior_support_by_UREH.csv"
)