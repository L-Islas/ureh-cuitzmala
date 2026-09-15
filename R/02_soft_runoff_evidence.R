# ============================================================
# 02. SOFT RUNOFF EVIDENCE
# ============================================================
#
# Purpose
# -------
# Generate the soft evidence for the runoff node used in the
# Bayesian-network inference.
#
# The evidence is derived from continuous Esc_sum using the
# piecewise membership functions of the historical workflow.
#
# The raw membership vector contains three components:
#
#   ev_esc_baja
#   ev_esc_media
#   ev_esc_alta
#
# and sums to 1 for every microcatchment.
#
# A regularization parameter epsilon is subsequently applied:
#
#   p_adj = epsilon + (1 - K * epsilon) * p
#
# where K = 3 runoff states.
#
# This prevents exact zero weights while preserving a total
# evidence mass equal to 1.
#
# Inputs
# ------
# temp/ureh_prepared.csv
#
# Parameters are read from R/00_config.R:
#
# b1_esc
# b2_esc
# m1_esc
# m2_esc
# epsilon_esc
#
# Outputs
# -------
# temp/soft_runoff_evidence.csv
# temp/soft_runoff_parameters.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------

source("R/00_config.R")

temp_dir <- here::here("temp")

dir.create(
  temp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Validate required configuration parameters
# ------------------------------------------------------------

required_parameters <- c(
  "b1_esc",
  "b2_esc",
  "m1_esc",
  "m2_esc",
  "epsilon_esc"
)

missing_parameters <- required_parameters[
  !vapply(
    required_parameters,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_parameters) > 0L) {
  stop(
    "Missing runoff-evidence parameters in R/00_config.R: ",
    paste(
      missing_parameters,
      collapse = ", "
    )
  )
}

parameter_values <- c(
  b1_esc = b1_esc,
  b2_esc = b2_esc,
  m1_esc = m1_esc,
  m2_esc = m2_esc,
  epsilon_esc = epsilon_esc
)

if (
  anyNA(parameter_values) ||
  any(!is.finite(parameter_values))
) {
  stop(
    "Runoff-evidence parameters must be finite numeric values."
  )
}

if (!(
  b1_esc < b2_esc &&
  b2_esc < m1_esc &&
  m1_esc < m2_esc
)) {
  stop(
    "Runoff thresholds must satisfy: ",
    "b1_esc < b2_esc < m1_esc < m2_esc."
  )
}

if (
  epsilon_esc <= 0 ||
  epsilon_esc >= 1 / 3
) {
  stop(
    "epsilon_esc must satisfy 0 < epsilon_esc < 1/3."
  )
}

if (!identical(
  estados$escorrentia,
  c(
    "esc_baja",
    "esc_media",
    "esc_alta"
  )
)) {
  stop(
    "Unexpected runoff-state definition in R/00_config.R."
  )
}

n_runoff_states <- length(
  estados$escorrentia
)

message(
  "Runoff-evidence parameters validated."
)


# ------------------------------------------------------------
# 2. Read prepared analytical dataset
# ------------------------------------------------------------

prepared_file <- here::here(
  "temp",
  "ureh_prepared.csv"
)

if (!file.exists(prepared_file)) {
  stop(
    "Prepared dataset not found: ",
    prepared_file,
    "\nRun R/01_prepare_data.R first."
  )
}

ureh <- read.csv(
  prepared_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Validate prepared dataset
# ------------------------------------------------------------

required_fields <- c(
  "id",
  "Esc_sum",
  "Esc_cat",
  "escorrentia"
)

missing_fields <- setdiff(
  required_fields,
  names(ureh)
)

if (length(missing_fields) > 0L) {
  stop(
    "Missing required fields in prepared dataset: ",
    paste(
      missing_fields,
      collapse = ", "
    )
  )
}

if (nrow(ureh) != n_expected_units) {
  stop(
    "Expected ",
    n_expected_units,
    " microcatchments, but found ",
    nrow(ureh),
    "."
  )
}

if (anyNA(ureh$id)) {
  stop(
    "Missing microcatchment identifiers detected."
  )
}

if (anyDuplicated(ureh$id)) {
  stop(
    "Duplicate microcatchment identifiers detected."
  )
}

if (!identical(
  sort(as.integer(ureh$id)),
  seq_len(n_expected_units)
)) {
  stop(
    "Microcatchment ids do not correspond to the expected ",
    "sequence 1:",
    n_expected_units,
    "."
  )
}

if (!is.numeric(ureh$Esc_sum)) {
  stop(
    "Esc_sum must be numeric."
  )
}

if (
  anyNA(ureh$Esc_sum) ||
  any(!is.finite(ureh$Esc_sum))
) {
  stop(
    "Esc_sum contains missing or non-finite values."
  )
}

if (any(ureh$Esc_sum < 0)) {
  stop(
    "Esc_sum contains negative values."
  )
}

if (anyNA(ureh$Esc_cat)) {
  stop(
    "Esc_cat contains missing values."
  )
}

if (anyNA(ureh$escorrentia)) {
  stop(
    "The runoff-state variable contains missing values."
  )
}

if (!all(
  ureh$escorrentia %in%
  estados$escorrentia
)) {
  stop(
    "Invalid runoff states detected in the prepared dataset."
  )
}


# ------------------------------------------------------------
# 4. Verify Esc_cat -> runoff-state correspondence
# ------------------------------------------------------------

expected_runoff_state <- c(
  "Baja" = "esc_baja",
  "Media" = "esc_media",
  "Alta" = "esc_alta"
)

invalid_esc_cat <- setdiff(
  unique(ureh$Esc_cat),
  names(expected_runoff_state)
)

if (length(invalid_esc_cat) > 0L) {
  stop(
    "Unexpected Esc_cat values: ",
    paste(
      invalid_esc_cat,
      collapse = ", "
    )
  )
}

derived_runoff_state <- unname(
  expected_runoff_state[
    ureh$Esc_cat
  ]
)

if (!identical(
  derived_runoff_state,
  ureh$escorrentia
)) {
  stop(
    "Esc_cat and escorrentia are not internally consistent."
  )
}

message(
  "Prepared runoff variables validated."
)


# ------------------------------------------------------------
# 5. Historical runoff membership functions
# ------------------------------------------------------------
#
# Low runoff:
#
#   1                       x <= b1
#   (b2 - x)/(b2 - b1)      b1 < x < b2
#   0                       x >= b2
#
# ------------------------------------------------------------

mu_esc_baja <- function(x) {
  
  if (!is.finite(x)) {
    stop(
      "mu_esc_baja received a non-finite value."
    )
  }
  
  if (x <= b1_esc) {
    
    return(1)
    
  } else if (x < b2_esc) {
    
    return(
      (b2_esc - x) /
        (b2_esc - b1_esc)
    )
    
  } else {
    
    return(0)
  }
}


# ------------------------------------------------------------
# Medium runoff
# ------------------------------------------------------------
#
#   0                       x <= b1
#   (x - b1)/(b2 - b1)      b1 < x < b2
#   1                       b2 <= x <= m1
#   (m2 - x)/(m2 - m1)      m1 < x < m2
#   0                       x >= m2
#
# ------------------------------------------------------------

mu_esc_media <- function(x) {
  
  if (!is.finite(x)) {
    stop(
      "mu_esc_media received a non-finite value."
    )
  }
  
  if (x <= b1_esc) {
    
    return(0)
    
  } else if (x < b2_esc) {
    
    return(
      (x - b1_esc) /
        (b2_esc - b1_esc)
    )
    
  } else if (x <= m1_esc) {
    
    return(1)
    
  } else if (x < m2_esc) {
    
    return(
      (m2_esc - x) /
        (m2_esc - m1_esc)
    )
    
  } else {
    
    return(0)
  }
}


# ------------------------------------------------------------
# High runoff
# ------------------------------------------------------------
#
#   0                       x <= m1
#   (x - m1)/(m2 - m1)      m1 < x < m2
#   1                       x >= m2
#
# ------------------------------------------------------------

mu_esc_alta <- function(x) {
  
  if (!is.finite(x)) {
    stop(
      "mu_esc_alta received a non-finite value."
    )
  }
  
  if (x <= m1_esc) {
    
    return(0)
    
  } else if (x < m2_esc) {
    
    return(
      (x - m1_esc) /
        (m2_esc - m1_esc)
    )
    
  } else {
    
    return(1)
  }
}


# ------------------------------------------------------------
# 6. Validate membership functions at threshold boundaries
# ------------------------------------------------------------

boundary_values <- c(
  b1_esc,
  b2_esc,
  m1_esc,
  m2_esc
)

boundary_membership <- data.frame(
  x = boundary_values,
  
  esc_baja = vapply(
    boundary_values,
    mu_esc_baja,
    numeric(1)
  ),
  
  esc_media = vapply(
    boundary_values,
    mu_esc_media,
    numeric(1)
  ),
  
  esc_alta = vapply(
    boundary_values,
    mu_esc_alta,
    numeric(1)
  )
)

expected_boundary_membership <- rbind(
  c(1, 0, 0),
  c(0, 1, 0),
  c(0, 1, 0),
  c(0, 0, 1)
)

observed_boundary_membership <- as.matrix(
  boundary_membership[
    ,
    c(
      "esc_baja",
      "esc_media",
      "esc_alta"
    )
  ]
)

if (!isTRUE(
  all.equal(
    observed_boundary_membership,
    expected_boundary_membership,
    tolerance = 1e-12,
    check.attributes = FALSE
  )
)) {
  stop(
    "Runoff membership functions failed the threshold-boundary check."
  )
}

message(
  "Runoff membership functions validated at all thresholds."
)


# ------------------------------------------------------------
# 7. Calculate raw soft evidence
# ------------------------------------------------------------

ev_esc_baja <- vapply(
  ureh$Esc_sum,
  mu_esc_baja,
  numeric(1)
)

ev_esc_media <- vapply(
  ureh$Esc_sum,
  mu_esc_media,
  numeric(1)
)

ev_esc_alta <- vapply(
  ureh$Esc_sum,
  mu_esc_alta,
  numeric(1)
)

raw_evidence <- cbind(
  ev_esc_baja,
  ev_esc_media,
  ev_esc_alta
)

raw_evidence_sum <- rowSums(
  raw_evidence
)


# ------------------------------------------------------------
# 8. Validate raw evidence
# ------------------------------------------------------------

probability_tolerance_local <- 1e-10

if (
  anyNA(raw_evidence) ||
  any(!is.finite(raw_evidence))
) {
  stop(
    "Raw runoff evidence contains missing or non-finite values."
  )
}

if (
  any(
    raw_evidence <
    -probability_tolerance_local
  ) ||
  any(
    raw_evidence >
    1 + probability_tolerance_local
  )
) {
  stop(
    "Raw runoff evidence contains values outside [0, 1]."
  )
}

max_raw_sum_deviation <- max(
  abs(
    raw_evidence_sum - 1
  )
)

if (
  max_raw_sum_deviation >
  probability_tolerance_local
) {
  stop(
    "Raw runoff evidence does not sum to 1. ",
    "Maximum deviation = ",
    max_raw_sum_deviation,
    "."
  )
}

message(
  "Raw runoff evidence validated."
)


# ------------------------------------------------------------
# 9. Regularize soft evidence
# ------------------------------------------------------------
#
# Historical regularization:
#
# p_adj = epsilon + (1 - K * epsilon) * p
#
# K = 3 runoff states.
#
# With epsilon = 0.05:
#
# raw 0 -> adjusted 0.05
# raw 1 -> adjusted 0.90
#
# ------------------------------------------------------------

regularization_factor <-
  1 -
  n_runoff_states *
  epsilon_esc

ev_esc_baja_aj <-
  epsilon_esc +
  regularization_factor *
  ev_esc_baja

ev_esc_media_aj <-
  epsilon_esc +
  regularization_factor *
  ev_esc_media

ev_esc_alta_aj <-
  epsilon_esc +
  regularization_factor *
  ev_esc_alta

adjusted_evidence <- cbind(
  ev_esc_baja_aj,
  ev_esc_media_aj,
  ev_esc_alta_aj
)

adjusted_evidence_sum <- rowSums(
  adjusted_evidence
)


# ------------------------------------------------------------
# 10. Validate regularized evidence
# ------------------------------------------------------------

if (
  anyNA(adjusted_evidence) ||
  any(!is.finite(adjusted_evidence))
) {
  stop(
    "Regularized runoff evidence contains missing or ",
    "non-finite values."
  )
}

if (
  any(
    adjusted_evidence <
    epsilon_esc -
    probability_tolerance_local
  )
) {
  stop(
    "Regularized runoff evidence contains values below epsilon."
  )
}

maximum_adjusted_weight <-
  1 -
  (n_runoff_states - 1) *
  epsilon_esc

if (
  any(
    adjusted_evidence >
    maximum_adjusted_weight +
    probability_tolerance_local
  )
) {
  stop(
    "Regularized runoff evidence exceeds its theoretical maximum."
  )
}

max_adjusted_sum_deviation <- max(
  abs(
    adjusted_evidence_sum - 1
  )
)

if (
  max_adjusted_sum_deviation >
  probability_tolerance_local
) {
  stop(
    "Regularized runoff evidence does not sum to 1. ",
    "Maximum deviation = ",
    max_adjusted_sum_deviation,
    "."
  )
}

message(
  "Regularized runoff evidence validated."
)


# ------------------------------------------------------------
# 11. Assemble output table
# ------------------------------------------------------------

soft_runoff_evidence <- data.frame(
  
  id = ureh$id,
  
  Esc_sum = ureh$Esc_sum,
  
  Esc_cat = ureh$Esc_cat,
  
  escorrentia = ureh$escorrentia,
  
  ev_esc_baja = ev_esc_baja,
  
  ev_esc_media = ev_esc_media,
  
  ev_esc_alta = ev_esc_alta,
  
  suma_evidencia = raw_evidence_sum,
  
  ev_esc_baja_aj = ev_esc_baja_aj,
  
  ev_esc_media_aj = ev_esc_media_aj,
  
  ev_esc_alta_aj = ev_esc_alta_aj,
  
  suma_evidencia_aj = adjusted_evidence_sum,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Parameter provenance table
# ------------------------------------------------------------

soft_runoff_parameters <- data.frame(
  
  parameter = c(
    "b1_esc",
    "b2_esc",
    "m1_esc",
    "m2_esc",
    "epsilon_esc"
  ),
  
  value = c(
    b1_esc,
    b2_esc,
    m1_esc,
    m2_esc,
    epsilon_esc
  ),
  
  role = c(
    "Upper limit of full low-runoff membership",
    "Start of full medium-runoff membership",
    "End of full medium-runoff membership",
    "Start of full high-runoff membership",
    "Regularization applied to the three runoff evidence weights"
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 13. Final object-level validation
# ------------------------------------------------------------

if (
  nrow(soft_runoff_evidence) !=
  n_expected_units
) {
  stop(
    "Unexpected number of rows in the runoff-evidence table."
  )
}

if (
  anyDuplicated(
    soft_runoff_evidence$id
  )
) {
  stop(
    "Duplicate ids detected in the runoff-evidence table."
  )
}

if (
  anyNA(
    soft_runoff_evidence
  )
) {
  stop(
    "The final runoff-evidence table contains missing values."
  )
}


# ------------------------------------------------------------
# 14. Write outputs
# ------------------------------------------------------------

evidence_file <- here::here(
  "temp",
  "soft_runoff_evidence.csv"
)

parameters_file <- here::here(
  "temp",
  "soft_runoff_parameters.csv"
)

write.csv(
  soft_runoff_evidence,
  evidence_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  soft_runoff_parameters,
  parameters_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 15. Read-back verification
# ------------------------------------------------------------

evidence_check <- read.csv(
  evidence_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

parameters_check <- read.csv(
  parameters_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(evidence_check) !=
  n_expected_units
) {
  stop(
    "Read-back validation failed for soft_runoff_evidence.csv."
  )
}

if (
  anyDuplicated(
    evidence_check$id
  )
) {
  stop(
    "Duplicate ids detected after reading the evidence file."
  )
}

if (
  anyNA(evidence_check)
) {
  stop(
    "Missing values detected after reading the evidence file."
  )
}

readback_raw_sum <- rowSums(
  evidence_check[
    ,
    c(
      "ev_esc_baja",
      "ev_esc_media",
      "ev_esc_alta"
    )
  ]
)

readback_adjusted_sum <- rowSums(
  evidence_check[
    ,
    c(
      "ev_esc_baja_aj",
      "ev_esc_media_aj",
      "ev_esc_alta_aj"
    )
  ]
)

if (
  max(
    abs(
      readback_raw_sum - 1
    )
  ) >
  probability_tolerance_local
) {
  stop(
    "Raw evidence failed validation after CSV read-back."
  )
}

if (
  max(
    abs(
      readback_adjusted_sum - 1
    )
  ) >
  probability_tolerance_local
) {
  stop(
    "Regularized evidence failed validation after CSV read-back."
  )
}

if (
  nrow(parameters_check) !=
  length(required_parameters)
) {
  stop(
    "Unexpected parameter table after CSV read-back."
  )
}


# ------------------------------------------------------------
# 16. Diagnostic summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SOFT RUNOFF EVIDENCE SUMMARY\n"
)

cat(
  "============================================\n"
)

cat(
  "Microcatchments:",
  nrow(soft_runoff_evidence),
  "\n"
)

cat(
  "Esc_sum range:",
  paste(
    range(ureh$Esc_sum),
    collapse = " - "
  ),
  "\n"
)

cat(
  "Raw evidence range:",
  paste(
    range(raw_evidence),
    collapse = " - "
  ),
  "\n"
)

cat(
  "Regularized evidence range:",
  paste(
    range(adjusted_evidence),
    collapse = " - "
  ),
  "\n"
)

cat(
  "Maximum raw sum deviation:",
  format(
    max_raw_sum_deviation,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Maximum regularized sum deviation:",
  format(
    max_adjusted_sum_deviation,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "epsilon:",
  epsilon_esc,
  "\n"
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 17. Completion
# ------------------------------------------------------------

message(
  "02_soft_runoff_evidence.R completed successfully."
)

message(
  "Generated: temp/soft_runoff_evidence.csv"
)

message(
  "Generated: temp/soft_runoff_parameters.csv"
)