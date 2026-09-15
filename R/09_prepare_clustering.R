# ============================================================
# 09. PREPARE CLUSTERING DATA
# ============================================================
#
# Purpose
# -------
# Prepare two analytically independent multivariate-similarity
# configurations used as a reference for the Bayesian
# regionalization.
#
# C1 — local bio-physical similarity
#
#   USyV
#   Suelos
#   Grad_Pend
#   P_sum
#   Esc_sum
#
# C2 — C1 + upstream contribution
#
#   C1 variables
#   Agua_D_log = log1p(Agua_D)
#
# IMPORTANT
# ---------
# The clustering does NOT use:
#
# - UREH classifications;
# - hydrological-function posteriors;
# - any Bayesian posterior probabilities;
# - geographic coordinates.
#
# P_sum missing values are retained without imputation.
# Gower distance will use the variables available for each
# pairwise comparison.
#
# Inputs
# ------
# temp/ureh_prepared.csv
#
# Intermediate outputs
# --------------------
# temp/clustering_C1.rds
# temp/clustering_C2.rds
# temp/clustering_base.csv
#
# Final transparent output
# ------------------------
# results/tables/clustering_variable_specification.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(dplyr)

temp_dir <- here::here(
  "temp"
)

results_tables_dir <- here::here(
  "results",
  "tables"
)

dir.create(
  temp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Input
# ------------------------------------------------------------

prepared_file <- here::here(
  "temp",
  "ureh_prepared.csv"
)

if (!file.exists(
  prepared_file
)) {
  stop(
    "Prepared dataset not found:\n",
    prepared_file,
    "\nRun R/01_prepare_data.R first."
  )
}

base <- read.csv(
  prepared_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. Required source variables
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
  names(base)
)

if (
  length(
    missing_fields
  ) > 0L
) {
  stop(
    "Missing fields required for clustering: ",
    paste(
      missing_fields,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 3. General validation
# ------------------------------------------------------------

if (
  nrow(base) !=
  n_expected_units
) {
  stop(
    "Expected ",
    n_expected_units,
    " units, found ",
    nrow(base),
    "."
  )
}

if (
  anyNA(base$id) ||
  anyDuplicated(base$id)
) {
  stop(
    "Invalid ids in clustering input."
  )
}

if (!identical(
  sort(
    as.integer(
      base$id
    )
  ),
  seq_len(
    n_expected_units
  )
)) {
  stop(
    "Unexpected ids in clustering input."
  )
}


# ------------------------------------------------------------
# 4. Validate missing-value pattern
# ------------------------------------------------------------
#
# P_sum may contain NA.
#
# All other clustering variables must be complete.
#
# Esc_sum has already received the same analysis-ready treatment
# used in the main workflow: the single source NA was replaced
# by zero in R/01_prepare_data.R.
#
# ------------------------------------------------------------

fields_without_allowed_NA <- c(
  "id",
  "USyV",
  "Suelos",
  "Grad_Pend",
  "Esc_sum",
  "Agua_D"
)

NA_summary <- sapply(
  base[
    ,
    required_fields,
    drop = FALSE
  ],
  function(x) {
    sum(
      is.na(x)
    )
  }
)

if (
  any(
    NA_summary[
      fields_without_allowed_NA
    ] > 0
  )
) {
  stop(
    "Missing values found in variables where NA is not allowed: ",
    paste(
      names(
        NA_summary[
          fields_without_allowed_NA
        ]
      ),
      NA_summary[
        fields_without_allowed_NA
      ],
      sep = "=",
      collapse = "; "
    )
  )
}

n_missing_precipitation <-
  unname(
    NA_summary[
      "P_sum"
    ]
  )

message(
  "P_sum missing values retained: ",
  n_missing_precipitation,
  " units."
)


# ------------------------------------------------------------
# 5. Prepare nominal variables
# ------------------------------------------------------------

base$USyV <- factor(
  base$USyV
)

base$Suelos <- factor(
  base$Suelos
)

if (
  nlevels(
    base$USyV
  ) < 2L
) {
  stop(
    "USyV does not contain enough levels for clustering."
  )
}

if (
  nlevels(
    base$Suelos
  ) < 2L
) {
  stop(
    "Suelos does not contain enough levels for clustering."
  )
}


# ------------------------------------------------------------
# 6. Prepare ordinal slope
# ------------------------------------------------------------

slope_levels <- c(
  "<1°",
  "1°-3°",
  "3°-6°",
  "6°-15°",
  "15°-35°"
)

observed_slope_levels <- unique(
  as.character(
    base$Grad_Pend
  )
)

unrecognized_slope_levels <- setdiff(
  observed_slope_levels,
  slope_levels
)

if (
  length(
    unrecognized_slope_levels
  ) > 0L
) {
  stop(
    "Grad_Pend contains unrecognized classes: ",
    paste(
      unrecognized_slope_levels,
      collapse = ", "
    )
  )
}

base$Grad_Pend <- ordered(
  base$Grad_Pend,
  levels =
    slope_levels
)

if (
  anyNA(
    base$Grad_Pend
  )
) {
  stop(
    "Ordinal slope conversion generated missing values."
  )
}


# ------------------------------------------------------------
# 7. Prepare numerical variables
# ------------------------------------------------------------

base$P_sum <- suppressWarnings(
  as.numeric(
    base$P_sum
  )
)

base$Esc_sum <- suppressWarnings(
  as.numeric(
    base$Esc_sum
  )
)

base$Agua_D <- suppressWarnings(
  as.numeric(
    base$Agua_D
  )
)


# P_sum may contain the original 13 NA values.

if (
  sum(
    is.na(
      base$P_sum
    )
  ) !=
  n_missing_precipitation
) {
  stop(
    "Numeric conversion changed the missing-value pattern of P_sum."
  )
}


if (
  anyNA(
    base$Esc_sum
  ) ||
  any(
    !is.finite(
      base$Esc_sum
    )
  )
) {
  stop(
    "Invalid Esc_sum values after numeric conversion."
  )
}


if (
  anyNA(
    base$Agua_D
  ) ||
  any(
    !is.finite(
      base$Agua_D
    )
  )
) {
  stop(
    "Invalid Agua_D values after numeric conversion."
  )
}


if (
  any(
    base$Agua_D < 0
  )
) {
  stop(
    "Agua_D contains negative values; log1p cannot be used."
  )
}


# ------------------------------------------------------------
# 8. Upstream-contribution transformation
# ------------------------------------------------------------

base$Agua_D_log <- log1p(
  base$Agua_D
)

if (
  anyNA(
    base$Agua_D_log
  ) ||
  any(
    !is.finite(
      base$Agua_D_log
    )
  )
) {
  stop(
    "Invalid values generated by log1p(Agua_D)."
  )
}


# ------------------------------------------------------------
# 9. Construct C1
# ------------------------------------------------------------
#
# C1 represents local similarity only.
#
# ------------------------------------------------------------

clustering_C1 <- base |>
  select(
    id,
    USyV,
    Suelos,
    Grad_Pend,
    P_sum,
    Esc_sum
  )


# ------------------------------------------------------------
# 10. Construct C2
# ------------------------------------------------------------
#
# C2 differs from C1 only by the addition of upstream
# contribution.
#
# ------------------------------------------------------------

clustering_C2 <- base |>
  select(
    id,
    USyV,
    Suelos,
    Grad_Pend,
    P_sum,
    Esc_sum,
    Agua_D_log
  )


# ------------------------------------------------------------
# 11. Structural relationship between C1 and C2
# ------------------------------------------------------------

C1_variables <- setdiff(
  names(
    clustering_C1
  ),
  "id"
)

C2_variables <- setdiff(
  names(
    clustering_C2
  ),
  "id"
)

if (!setequal(
  C2_variables,
  c(
    C1_variables,
    "Agua_D_log"
  )
)) {
  stop(
    "C2 must correspond exactly to C1 + Agua_D_log."
  )
}


if (!identical(
  clustering_C1$id,
  clustering_C2$id
)) {
  stop(
    "C1 and C2 ids differ."
  )
}


# ------------------------------------------------------------
# 12. Validate variable classes
# ------------------------------------------------------------

expected_C1_classes <- c(
  USyV =
    "factor",
  
  Suelos =
    "factor",
  
  Grad_Pend =
    "ordered",
  
  P_sum =
    "numeric",
  
  Esc_sum =
    "numeric"
)


if (!is.factor(
  clustering_C1$USyV
)) {
  stop(
    "USyV must be a nominal factor."
  )
}

if (
  is.ordered(
    clustering_C1$USyV
  )
) {
  stop(
    "USyV must be nominal, not ordered."
  )
}


if (!is.factor(
  clustering_C1$Suelos
)) {
  stop(
    "Suelos must be a nominal factor."
  )
}

if (
  is.ordered(
    clustering_C1$Suelos
  )
) {
  stop(
    "Suelos must be nominal, not ordered."
  )
}


if (!is.ordered(
  clustering_C1$Grad_Pend
)) {
  stop(
    "Grad_Pend must be an ordered factor."
  )
}


if (!is.numeric(
  clustering_C1$P_sum
)) {
  stop(
    "P_sum must be numeric."
  )
}


if (!is.numeric(
  clustering_C1$Esc_sum
)) {
  stop(
    "Esc_sum must be numeric."
  )
}


if (!is.numeric(
  clustering_C2$Agua_D_log
)) {
  stop(
    "Agua_D_log must be numeric."
  )
}


# ------------------------------------------------------------
# 13. Confirm analytical independence of clustering inputs
# ------------------------------------------------------------
#
# C1 and C2 must contain only the variables explicitly defined
# for the multivariate-similarity analysis.
#
# No UREH classifications, Bayesian posterior probabilities,
# hydrological-function outputs, or spatial coordinates are
# allowed.
#
# ------------------------------------------------------------

expected_C1_fields <- c(
  "id",
  "USyV",
  "Suelos",
  "Grad_Pend",
  "P_sum",
  "Esc_sum"
)

expected_C2_fields <- c(
  "id",
  "USyV",
  "Suelos",
  "Grad_Pend",
  "P_sum",
  "Esc_sum",
  "Agua_D_log"
)


if (!identical(
  names(clustering_C1),
  expected_C1_fields
)) {
  stop(
    paste0(
      "C1 contains unexpected fields.\n",
      "Expected: ",
      paste(
        expected_C1_fields,
        collapse = ", "
      ),
      "\nFound: ",
      paste(
        names(clustering_C1),
        collapse = ", "
      )
    )
  )
}


if (!identical(
  names(clustering_C2),
  expected_C2_fields
)) {
  stop(
    paste0(
      "C2 contains unexpected fields.\n",
      "Expected: ",
      paste(
        expected_C2_fields,
        collapse = ", "
      ),
      "\nFound: ",
      paste(
        names(clustering_C2),
        collapse = ", "
      )
    )
  )
}


message(
  "Clustering inputs contain only the intended independent variables."
)

# ------------------------------------------------------------
# 14. Variable specification table
# ------------------------------------------------------------

clustering_variable_specification <- data.frame(
  
  variable = c(
    "USyV",
    "Suelos",
    "Grad_Pend",
    "P_sum",
    "Esc_sum",
    "Agua_D_log"
  ),
  
  type = c(
    "nominal",
    "nominal",
    "ordinal",
    "numeric",
    "numeric",
    "numeric"
  ),
  
  included_in = c(
    "C1 and C2",
    "C1 and C2",
    "C1 and C2",
    "C1 and C2",
    "C1 and C2",
    "C2 only"
  ),
  
  treatment = c(
    "Original categories; nominal factor",
    "Original categories; nominal factor",
    "Ordered factor; five original slope classes",
    "Original continuous scale; missing values retained",
    "Original continuous scale",
    "log1p(Agua_D)"
  ),
  
  gower_weight = rep(
    1,
    6
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 15. Missing-precipitation units
# ------------------------------------------------------------

missing_precipitation_units <- base |>
  filter(
    is.na(
      P_sum
    )
  ) |>
  select(
    id
  )


# ------------------------------------------------------------
# 16. Transparent clustering base
# ------------------------------------------------------------

clustering_base <- base |>
  select(
    id,
    USyV,
    Suelos,
    Grad_Pend,
    P_sum,
    Esc_sum,
    Agua_D,
    Agua_D_log
  )


# ------------------------------------------------------------
# 17. Output paths
# ------------------------------------------------------------

C1_file <- here::here(
  "temp",
  "clustering_C1.rds"
)

C2_file <- here::here(
  "temp",
  "clustering_C2.rds"
)

base_file <- here::here(
  "temp",
  "clustering_base.csv"
)

missing_precipitation_file <- here::here(
  "temp",
  "clustering_missing_P_sum.csv"
)

specification_file <- here::here(
  "results",
  "tables",
  "clustering_variable_specification.csv"
)


# ------------------------------------------------------------
# 18. Write outputs
# ------------------------------------------------------------

saveRDS(
  clustering_C1,
  C1_file
)

saveRDS(
  clustering_C2,
  C2_file
)

write.csv(
  clustering_base,
  base_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  missing_precipitation_units,
  missing_precipitation_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  clustering_variable_specification,
  specification_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 19. Read-back validation
# ------------------------------------------------------------

C1_check <- readRDS(
  C1_file
)

C2_check <- readRDS(
  C2_file
)

base_check <- read.csv(
  base_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

specification_check <- read.csv(
  specification_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (
  nrow(
    C1_check
  ) !=
  n_expected_units ||
  nrow(
    C2_check
  ) !=
  n_expected_units
) {
  stop(
    "Clustering RDS files failed read-back validation."
  )
}


if (!identical(
  C1_check$id,
  C2_check$id
)) {
  stop(
    "C1 and C2 ids differ after read-back."
  )
}


if (
  nrow(
    base_check
  ) !=
  n_expected_units
) {
  stop(
    "Clustering base failed read-back validation."
  )
}


if (
  nrow(
    specification_check
  ) != 6L
) {
  stop(
    "Variable specification failed read-back validation."
  )
}


if (!is.factor(
  C1_check$USyV
)) {
  stop(
    "USyV factor class was not preserved in RDS."
  )
}


if (!is.factor(
  C1_check$Suelos
)) {
  stop(
    "Suelos factor class was not preserved in RDS."
  )
}


if (!is.ordered(
  C1_check$Grad_Pend
)) {
  stop(
    "Grad_Pend ordered-factor class was not preserved in RDS."
  )
}


# ------------------------------------------------------------
# 20. Console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "CLUSTERING DATA PREPARATION\n"
)

cat(
  "============================================\n"
)

cat(
  "Units:",
  n_expected_units,
  "\n"
)

cat(
  "C1 variables:",
  paste(
    C1_variables,
    collapse = " | "
  ),
  "\n"
)

cat(
  "C2 variables:",
  paste(
    C2_variables,
    collapse = " | "
  ),
  "\n"
)

cat(
  "P_sum missing values retained:",
  n_missing_precipitation,
  "\n\n"
)


cat(
  "Variable classes in C1:\n"
)

print(
  sapply(
    clustering_C1[
      ,
      -1,
      drop = FALSE
    ],
    function(x) {
      paste(
        class(x),
        collapse = "/"
      )
    }
  )
)


cat(
  "\nVariable classes in C2:\n"
)

print(
  sapply(
    clustering_C2[
      ,
      -1,
      drop = FALSE
    ],
    function(x) {
      paste(
        class(x),
        collapse = "/"
      )
    }
  )
)


cat(
  "\nMissing values in C1:\n"
)

print(
  sapply(
    clustering_C1[
      ,
      -1,
      drop = FALSE
    ],
    function(x) {
      sum(
        is.na(x)
      )
    }
  )
)


cat(
  "\nMissing values in C2:\n"
)

print(
  sapply(
    clustering_C2[
      ,
      -1,
      drop = FALSE
    ],
    function(x) {
      sum(
        is.na(x)
      )
    }
  )
)


cat(
  "\nNumerical-variable summary:\n"
)

print(
  summary(
    clustering_C2[
      ,
      c(
        "P_sum",
        "Esc_sum",
        "Agua_D_log"
      )
    ]
  )
)


cat(
  "\n============================================\n\n"
)


# ------------------------------------------------------------
# 21. Completion
# ------------------------------------------------------------

message(
  "09_prepare_clustering.R completed successfully."
)

message(
  "Generated: temp/clustering_C1.rds"
)

message(
  "Generated: temp/clustering_C2.rds"
)

message(
  "Generated: temp/clustering_base.csv"
)

message(
  "Generated: temp/clustering_missing_P_sum.csv"
)

message(
  "Generated: results/tables/clustering_variable_specification.csv"
)

message(
  "Next: R/10_gower_PAM.R"
)