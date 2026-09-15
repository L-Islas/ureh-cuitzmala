# ============================================================
# 01. PREPARE ANALYSIS DATA
# ============================================================
#
# Purpose
# -------
# Read and validate the analysis-ready spatial dataset and generate
# the variables required by the Bayesian-network and multivariate
# analyses.
#
# Main operations
# ----------------
# 1. Read the authoritative microcatchment GeoPackage.
# 2. Verify identifiers, CRS and documented missing values.
# 3. Apply the documented Esc_sum NA -> 0 correction.
# 4. Generate log1p(Agua_D).
# 5. Discretize precipitation, runoff and upstream contribution.
# 6. Recode land cover, soil and slope to Bayesian-network states.
# 7. Apply the documented historical treatment of missing
#    precipitation categories.
# 8. Validate all Bayesian-network input states.
# 9. Write reproducible intermediate tables to temp/.
#
# Important
# ---------
# P_sum missing values are preserved in the continuous variable.
# Only the derived Bayesian-network precipitation state receives
# the documented historical assignment to "prec_baja".
#
# No spatial output is created here. Later mapping scripts should
# join analytical results to data/input/ureh_units.gpkg by `id`.
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(sf)
library(dplyr)
library(classInt)

dir.create(
  here::here("temp"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Read analysis-ready spatial input
# ------------------------------------------------------------

ureh <- st_read(
  file_ureh_units,
  layer = "ureh_units",
  quiet = TRUE
)


# ------------------------------------------------------------
# 2. Input integrity checks
# ------------------------------------------------------------

required_fields <- c(
  "id",
  "USyV",
  "Suelos",
  "Grad_Pend",
  "P_sum",
  "Esc_sum",
  "Agua_D"
)

missing_fields <- setdiff(
  required_fields,
  names(ureh)
)

if (length(missing_fields) > 0L) {
  stop(
    "Missing required input fields: ",
    paste(missing_fields, collapse = ", ")
  )
}

if (nrow(ureh) != n_expected_units) {
  stop(
    "Expected ",
    n_expected_units,
    " spatial units, but found ",
    nrow(ureh),
    "."
  )
}

if (anyNA(ureh$id)) {
  stop(
    "The input contains missing microcatchment identifiers."
  )
}

if (anyDuplicated(ureh$id)) {
  stop(
    "Microcatchment identifiers are not unique."
  )
}

if (!identical(
  sort(as.integer(ureh$id)),
  seq_len(n_expected_units)
)) {
  stop(
    "Microcatchment identifiers do not correspond to ",
    "the expected sequence 1:",
    n_expected_units,
    "."
  )
}

if (
  is.na(st_crs(ureh)$epsg) ||
  st_crs(ureh)$epsg != expected_epsg
) {
  stop(
    "Unexpected CRS. Expected EPSG:",
    expected_epsg,
    "."
  )
}


# ------------------------------------------------------------
# 3. Verify source categories
# ------------------------------------------------------------

expected_USyV <- c(
  "Bosque de Encino",
  "Bosque Tropical",
  "Sec Bosque Templado",
  "Sec Bosque Tropical",
  "Pastizal",
  "Agricultura",
  "Manglar",
  "Urbano"
)

expected_Suelos <- c(
  "Regosol",
  "Phaeozem",
  "Cambisol",
  "Luvisol",
  "Leptosol",
  "Arenosol",
  "Solonchank"
)

expected_Grad_Pend <- c(
  "6°-15°",
  "3°-6°",
  "15°-35°",
  "1°-3°",
  "<1°"
)

if (!setequal(
  unique(ureh$USyV),
  expected_USyV
)) {
  stop(
    "Unexpected categories detected in USyV."
  )
}

if (!setequal(
  unique(ureh$Suelos),
  expected_Suelos
)) {
  stop(
    "Unexpected categories detected in Suelos."
  )
}

if (!setequal(
  unique(ureh$Grad_Pend),
  expected_Grad_Pend
)) {
  stop(
    "Unexpected categories detected in Grad_Pend."
  )
}


# ------------------------------------------------------------
# 4. Verify documented missing-value pattern
# ------------------------------------------------------------

na_counts <- sapply(
  st_drop_geometry(ureh)[required_fields],
  function(x) sum(is.na(x))
)

cat("\nMissing values in source input:\n")
print(na_counts)

if (na_counts[["P_sum"]] != 13L) {
  stop(
    "Unexpected number of missing P_sum values. ",
    "Expected 13, found ",
    na_counts[["P_sum"]],
    "."
  )
}

if (na_counts[["Esc_sum"]] != 1L) {
  stop(
    "Unexpected number of missing Esc_sum values. ",
    "Expected 1, found ",
    na_counts[["Esc_sum"]],
    "."
  )
}

complete_fields <- c(
  "id",
  "USyV",
  "Suelos",
  "Grad_Pend",
  "Agua_D"
)

if (any(na_counts[complete_fields] != 0L)) {
  stop(
    "Unexpected missing values detected in fields ",
    "expected to be complete."
  )
}

message(
  "Input dataset verified: ",
  nrow(ureh),
  " microcatchments; EPSG:",
  st_crs(ureh)$epsg,
  "."
)

message(
  "Missing-value pattern matches the documented input dataset."
)


# ------------------------------------------------------------
# 5. Documented correction of Esc_sum
# ------------------------------------------------------------
#
# The historical analytical workflow documented that the single
# missing Esc_sum value corresponds to zero.
# This correction is applied explicitly before discretization.
# ------------------------------------------------------------

n_na_esc <- sum(is.na(ureh$Esc_sum))

ureh$Esc_sum[
  is.na(ureh$Esc_sum)
] <- 0

if (anyNA(ureh$Esc_sum)) {
  stop(
    "Esc_sum still contains missing values after correction."
  )
}

message(
  "Esc_sum: ",
  n_na_esc,
  " missing value replaced by 0 according to ",
  "the documented workflow."
)


# ------------------------------------------------------------
# 6. Upstream-contribution transformation
# ------------------------------------------------------------

ureh <- ureh |>
  mutate(
    Agua_D_log = log1p(Agua_D)
  )

if (
  anyNA(ureh$Agua_D_log) ||
  any(!is.finite(ureh$Agua_D_log))
) {
  stop(
    "Agua_D_log contains missing or non-finite values."
  )
}

message(
  "Agua_D_log generated successfully."
)


# ------------------------------------------------------------
# 7. Discretization thresholds
# ------------------------------------------------------------
#
# Precipitation:
#   quantile classification, 3 classes.
#
# Runoff:
#   Jenks natural breaks, 3 classes.
#
# Upstream contribution:
#   Jenks natural breaks on log1p(Agua_D), 4 classes.
#
# Missing P_sum values are explicitly excluded only when
# estimating the precipitation breaks. The original P_sum
# values remain missing in the analytical dataset.
# ------------------------------------------------------------

q_P <- classIntervals(
  ureh$P_sum[
    !is.na(ureh$P_sum)
  ],
  n = n_precipitacion,
  style = "quantile"
)

j_Esc <- classIntervals(
  ureh$Esc_sum,
  n = n_escorrentia,
  style = "jenks"
)

j_Agua <- classIntervals(
  ureh$Agua_D_log,
  n = n_aporte,
  style = "jenks"
)

cat("\nPrecipitation breaks:\n")
print(q_P$brks)

cat("\nRunoff breaks:\n")
print(j_Esc$brks)

cat("\nUpstream-contribution log breaks:\n")
print(j_Agua$brks)


# ------------------------------------------------------------
# 8. Discretized variables
# ------------------------------------------------------------

ureh <- ureh |>
  mutate(
    
    P_cat = cut(
      P_sum,
      breaks = q_P$brks,
      include.lowest = TRUE,
      labels = c(
        "Baja",
        "Media",
        "Alta"
      )
    ),
    
    Esc_cat = cut(
      Esc_sum,
      breaks = j_Esc$brks,
      include.lowest = TRUE,
      labels = c(
        "Baja",
        "Media",
        "Alta"
      )
    ),
    
    Agua_D_cat = cut(
      Agua_D_log,
      breaks = j_Agua$brks,
      include.lowest = TRUE,
      labels = c(
        "Muy Baja",
        "Baja",
        "Media",
        "Alta"
      )
    )
  )

discretization_na <- sapply(
  st_drop_geometry(ureh)[
    c(
      "P_cat",
      "Esc_cat",
      "Agua_D_cat"
    )
  ],
  function(x) sum(is.na(x))
)

cat("\nMissing values after discretization:\n")
print(discretization_na)

if (discretization_na[["P_cat"]] != 13L) {
  stop(
    "Expected 13 missing P_cat values after discretization, ",
    "but found ",
    discretization_na[["P_cat"]],
    "."
  )
}

if (discretization_na[["Esc_cat"]] != 0L) {
  stop(
    "Unexpected missing values in Esc_cat."
  )
}

if (discretization_na[["Agua_D_cat"]] != 0L) {
  stop(
    "Unexpected missing values in Agua_D_cat."
  )
}


# ------------------------------------------------------------
# 9. Recode variables to Bayesian-network states
# ------------------------------------------------------------

ureh <- ureh |>
  mutate(
    
    cobertura = recode(
      USyV,
      "Bosque de Encino" =
        "forestal_conservada",
      "Bosque Tropical" =
        "forestal_conservada",
      "Sec Bosque Templado" =
        "forestal_secundaria",
      "Sec Bosque Tropical" =
        "forestal_secundaria",
      "Pastizal" =
        "agropecuaria",
      "Agricultura" =
        "agropecuaria",
      "Manglar" =
        "hidrofila_humedal",
      "Urbano" =
        "antropica",
      .default = NA_character_
    ),
    
    suelo = recode(
      Suelos,
      "Regosol" =
        "suelo_infiltracion_baja",
      "Phaeozem" =
        "suelo_infiltracion_alta",
      "Cambisol" =
        "suelo_infiltracion_media",
      "Luvisol" =
        "suelo_infiltracion_media",
      "Leptosol" =
        "suelo_infiltracion_baja",
      "Arenosol" =
        "suelo_infiltracion_alta",
      "Solonchank" =
        "suelo_salino",
      .default = NA_character_
    ),
    
    pendiente = recode(
      Grad_Pend,
      "6°-15°" =
        "pend_media",
      "3°-6°" =
        "pend_media",
      "15°-35°" =
        "pend_alta",
      "1°-3°" =
        "pend_baja",
      "<1°" =
        "pend_baja",
      .default = NA_character_
    ),
    
    precipitacion = recode(
      as.character(P_cat),
      "Baja" =
        "prec_baja",
      "Media" =
        "prec_media",
      "Alta" =
        "prec_alta",
      .default = NA_character_
    ),
    
    escorrentia = recode(
      as.character(Esc_cat),
      "Baja" =
        "esc_baja",
      "Media" =
        "esc_media",
      "Alta" =
        "esc_alta",
      .default = NA_character_
    ),
    
    aporte = recode(
      as.character(Agua_D_cat),
      "Muy Baja" =
        "aporte_muy_bajo",
      "Baja" =
        "aporte_bajo",
      "Media" =
        "aporte_medio",
      "Alta" =
        "aporte_alto",
      .default = NA_character_
    )
  )


# ------------------------------------------------------------
# 10. Verify states before historical P_sum treatment
# ------------------------------------------------------------

network_variables <- c(
  "cobertura",
  "suelo",
  "pendiente",
  "precipitacion",
  "escorrentia",
  "aporte"
)

network_na_before <- sapply(
  st_drop_geometry(ureh)[network_variables],
  function(x) sum(is.na(x))
)

cat(
  "\nMissing Bayesian-network states before ",
  "precipitation reassignment:\n"
)

print(network_na_before)

expected_network_na_before <- c(
  cobertura = 0L,
  suelo = 0L,
  pendiente = 0L,
  precipitacion = 13L,
  escorrentia = 0L,
  aporte = 0L
)

if (!identical(
  as.integer(network_na_before),
  as.integer(expected_network_na_before)
)) {
  stop(
    "Unexpected missing-value pattern in ",
    "Bayesian-network states before precipitation reassignment."
  )
}


# ------------------------------------------------------------
# 11. Historical treatment of missing precipitation state
# ------------------------------------------------------------
#
# In the historical workflow, P_sum remained missing as a
# continuous variable, but the corresponding Bayesian-network
# precipitation state was assigned to "prec_baja".
# ------------------------------------------------------------

n_na_prec_state <- sum(
  is.na(ureh$precipitacion)
)

if (n_na_prec_state != 13L) {
  stop(
    "Expected 13 missing precipitation states before ",
    "historical reassignment, but found ",
    n_na_prec_state,
    "."
  )
}

ureh$precipitacion[
  is.na(ureh$precipitacion)
] <- "prec_baja"

message(
  n_na_prec_state,
  " missing precipitation states assigned to ",
  "'prec_baja' according to the documented ",
  "historical workflow."
)


# ------------------------------------------------------------
# 12. Validate Bayesian-network states
# ------------------------------------------------------------

if (!all(
  ureh$cobertura %in% estados$cobertura
)) {
  stop(
    "Invalid land-cover states detected."
  )
}

if (!all(
  ureh$suelo %in% estados$suelo
)) {
  stop(
    "Invalid soil states detected."
  )
}

if (!all(
  ureh$pendiente %in% estados$pendiente
)) {
  stop(
    "Invalid slope states detected."
  )
}

if (!all(
  ureh$precipitacion %in%
  estados$precipitacion
)) {
  stop(
    "Invalid precipitation states detected."
  )
}

if (!all(
  ureh$escorrentia %in%
  estados$escorrentia
)) {
  stop(
    "Invalid runoff states detected."
  )
}

if (!all(
  ureh$aporte %in% estados$aporte
)) {
  stop(
    "Invalid upstream-contribution states detected."
  )
}

network_na_after <- sapply(
  st_drop_geometry(ureh)[network_variables],
  function(x) sum(is.na(x))
)

cat(
  "\nMissing Bayesian-network states after preprocessing:\n"
)

print(network_na_after)

if (any(network_na_after != 0L)) {
  stop(
    "Missing values remain in Bayesian-network input states."
  )
}

message(
  "All Bayesian-network input states validated successfully."
)


# ------------------------------------------------------------
# 13. Verify continuous variables retained for clustering
# ------------------------------------------------------------

if (sum(is.na(ureh$P_sum)) != 13L) {
  stop(
    "P_sum missing values were unexpectedly modified."
  )
}

if (anyNA(ureh$Esc_sum)) {
  stop(
    "Esc_sum contains missing values."
  )
}

if (anyNA(ureh$Agua_D)) {
  stop(
    "Agua_D contains missing values."
  )
}

if (anyNA(ureh$Agua_D_log)) {
  stop(
    "Agua_D_log contains missing values."
  )
}

message(
  "Continuous variables for clustering preserved successfully."
)


# ------------------------------------------------------------
# 14. Create discretization-break table
# ------------------------------------------------------------

discretization_breaks <- bind_rows(
  
  data.frame(
    variable = "P_sum",
    method = "quantile",
    transformed = FALSE,
    break_order = seq_along(q_P$brks),
    break_value = as.numeric(q_P$brks)
  ),
  
  data.frame(
    variable = "Esc_sum",
    method = "jenks",
    transformed = FALSE,
    break_order = seq_along(j_Esc$brks),
    break_value = as.numeric(j_Esc$brks)
  ),
  
  data.frame(
    variable = "Agua_D_log",
    method = "jenks",
    transformed = TRUE,
    break_order = seq_along(j_Agua$brks),
    break_value = as.numeric(j_Agua$brks)
  )
)


# ------------------------------------------------------------
# 15. Write reproducible intermediate products
# ------------------------------------------------------------
#
# These files are generated products and belong in temp/.
# They should not be treated as primary input data.
# ------------------------------------------------------------

prepared_table <- st_drop_geometry(ureh)

write.csv(
  prepared_table,
  here::here(
    "temp",
    "ureh_prepared.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  discretization_breaks,
  here::here(
    "temp",
    "discretization_breaks.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 16. Final read-back checks
# ------------------------------------------------------------

prepared_check <- read.csv(
  here::here(
    "temp",
    "ureh_prepared.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

breaks_check <- read.csv(
  here::here(
    "temp",
    "discretization_breaks.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(prepared_check) != n_expected_units) {
  stop(
    "The written prepared dataset does not contain ",
    n_expected_units,
    " rows."
  )
}

if (anyDuplicated(prepared_check$id)) {
  stop(
    "Duplicate ids detected after writing the prepared dataset."
  )
}

if (nrow(breaks_check) !=
    length(q_P$brks) +
    length(j_Esc$brks) +
    length(j_Agua$brks)) {
  stop(
    "Unexpected number of discretization breaks after writing."
  )
}


# ------------------------------------------------------------
# 17. Completion message
# ------------------------------------------------------------

message(
  "01_prepare_data.R completed successfully."
)

message(
  "Prepared units: ",
  nrow(prepared_check)
)

message(
  "Generated: temp/ureh_prepared.csv"
)

message(
  "Generated: temp/discretization_breaks.csv"
)
