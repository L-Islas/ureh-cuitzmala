# ============================================================
# 04. INFER M0
# ============================================================
#
# Purpose
# -------
# Run the complete M0 inference for all microcatchments.
#
# M0 incorporates:
#
# - observed soil state;
# - observed land-cover state;
# - observed slope state;
# - observed precipitation state;
# - observed upstream-contribution state;
# - regularized soft evidence for runoff.
#
# Inference sequence
# ------------------
# 1. Infiltration conditional distribution.
# 2. Runoff distribution implied by the Bayesian network.
# 3. Runoff update using regularized soft evidence.
# 4. Hydrological-function distribution.
# 5. UREH posterior distribution.
# 6. MAP classification and posterior-support diagnostics.
#
# Inputs
# ------
# temp/ureh_prepared.csv
# temp/soft_runoff_evidence.csv
# temp/bayesian_network_model.rds
#
# Outputs
# -------
# temp/M0_UREH_resultados_completos.csv
#
# results/tables/M0_UREH_clasificacion_dominante.csv
# results/tables/M0_UREH_probabilidades.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(dplyr)
library(tidyr)

temp_dir <- here::here("temp")

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
# 1. Input paths
# ------------------------------------------------------------

prepared_file <- here::here(
  "temp",
  "ureh_prepared.csv"
)

evidence_file <- here::here(
  "temp",
  "soft_runoff_evidence.csv"
)

model_file <- here::here(
  "temp",
  "bayesian_network_model.rds"
)

required_files <- c(
  prepared_file,
  evidence_file,
  model_file
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
    ),
    "\nRun scripts 01-03 first."
  )
}


# ------------------------------------------------------------
# 2. Read inputs
# ------------------------------------------------------------

ureh <- read.csv(
  prepared_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

soft_runoff_evidence <- read.csv(
  evidence_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

model <- readRDS(
  model_file
)


# ------------------------------------------------------------
# 3. Validate model object
# ------------------------------------------------------------

required_model_objects <- c(
  "dag",
  "fit_ureh",
  "estados",
  "cpt_infiltracion",
  "cpt_escorrentia",
  "cpt_funcion_hidrologica",
  "cpt_ureh",
  "priors"
)

missing_model_objects <- setdiff(
  required_model_objects,
  names(model)
)

if (length(missing_model_objects) > 0L) {
  stop(
    "The Bayesian-network model object is incomplete. Missing: ",
    paste(
      missing_model_objects,
      collapse = ", "
    )
  )
}

estados_modelo <- model$estados

if (!identical(
  estados_modelo,
  estados
)) {
  stop(
    "Model states do not match R/00_config.R."
  )
}

cpt_infiltracion <-
  model$cpt_infiltracion

cpt_escorrentia <-
  model$cpt_escorrentia

cpt_funcion_hidrologica <-
  model$cpt_funcion_hidrologica

cpt_ureh <-
  model$cpt_ureh

cols_ureh <- estados$ureh

if (length(cols_ureh) != 6L) {
  stop(
    "M0 expects exactly six UREH states."
  )
}

message(
  "Bayesian-network model loaded and validated."
)


# ------------------------------------------------------------
# 4. Validate prepared dataset
# ------------------------------------------------------------

required_prepared_fields <- c(
  "id",
  "suelo",
  "cobertura",
  "pendiente",
  "precipitacion",
  "aporte"
)

missing_prepared_fields <- setdiff(
  required_prepared_fields,
  names(ureh)
)

if (length(missing_prepared_fields) > 0L) {
  stop(
    "Missing fields in prepared dataset: ",
    paste(
      missing_prepared_fields,
      collapse = ", "
    )
  )
}

if (nrow(ureh) != n_expected_units) {
  stop(
    "Expected ",
    n_expected_units,
    " prepared microcatchments, found ",
    nrow(ureh),
    "."
  )
}

if (anyNA(ureh$id)) {
  stop(
    "Prepared dataset contains missing ids."
  )
}

if (anyDuplicated(ureh$id)) {
  stop(
    "Prepared dataset contains duplicated ids."
  )
}

if (!identical(
  sort(as.integer(ureh$id)),
  seq_len(n_expected_units)
)) {
  stop(
    "Prepared-data ids do not correspond to 1:",
    n_expected_units,
    "."
  )
}


# ------------------------------------------------------------
# 5. Validate observed states
# ------------------------------------------------------------

validate_observed_states <- function(
    x,
    allowed_states,
    variable_name
) {
  
  x <- trimws(
    as.character(x)
  )
  
  if (anyNA(x)) {
    stop(
      variable_name,
      " contains missing states."
    )
  }
  
  invalid_states <- setdiff(
    unique(x),
    allowed_states
  )
  
  if (length(invalid_states) > 0L) {
    stop(
      variable_name,
      " contains invalid states: ",
      paste(
        invalid_states,
        collapse = ", "
      )
    )
  }
  
  invisible(TRUE)
}


validate_observed_states(
  ureh$suelo,
  estados$suelo,
  "suelo"
)

validate_observed_states(
  ureh$cobertura,
  estados$cobertura,
  "cobertura"
)

validate_observed_states(
  ureh$pendiente,
  estados$pendiente,
  "pendiente"
)

validate_observed_states(
  ureh$precipitacion,
  estados$precipitacion,
  "precipitacion"
)

validate_observed_states(
  ureh$aporte,
  estados$aporte,
  "aporte"
)

message(
  "Observed M0 states validated."
)


# ------------------------------------------------------------
# 6. Validate soft-runoff evidence
# ------------------------------------------------------------

required_evidence_fields <- c(
  "id",
  "ev_esc_baja_aj",
  "ev_esc_media_aj",
  "ev_esc_alta_aj"
)

missing_evidence_fields <- setdiff(
  required_evidence_fields,
  names(soft_runoff_evidence)
)

if (length(missing_evidence_fields) > 0L) {
  stop(
    "Missing fields in soft-runoff evidence: ",
    paste(
      missing_evidence_fields,
      collapse = ", "
    )
  )
}

if (
  nrow(soft_runoff_evidence) !=
  n_expected_units
) {
  stop(
    "Soft-runoff evidence does not contain ",
    n_expected_units,
    " microcatchments."
  )
}

if (
  anyNA(
    soft_runoff_evidence$id
  )
) {
  stop(
    "Soft-runoff evidence contains missing ids."
  )
}

if (
  anyDuplicated(
    soft_runoff_evidence$id
  )
) {
  stop(
    "Soft-runoff evidence contains duplicated ids."
  )
}

evidence_matrix <- as.matrix(
  soft_runoff_evidence[
    ,
    c(
      "ev_esc_baja_aj",
      "ev_esc_media_aj",
      "ev_esc_alta_aj"
    )
  ]
)

storage.mode(
  evidence_matrix
) <- "double"

if (
  anyNA(evidence_matrix) ||
  any(!is.finite(evidence_matrix))
) {
  stop(
    "Soft-runoff evidence contains invalid values."
  )
}

if (
  any(evidence_matrix <= 0)
) {
  stop(
    "Regularized runoff evidence must contain positive weights."
  )
}

evidence_sum_deviation <- max(
  abs(
    rowSums(evidence_matrix) - 1
  )
)

if (
  evidence_sum_deviation > 1e-10
) {
  stop(
    "Regularized runoff evidence does not sum to 1. ",
    "Maximum deviation = ",
    evidence_sum_deviation,
    "."
  )
}


# ------------------------------------------------------------
# 7. Match soft evidence to microcatchments by id
# ------------------------------------------------------------

evidence_index <- match(
  ureh$id,
  soft_runoff_evidence$id
)

if (anyNA(evidence_index)) {
  stop(
    "Soft-runoff evidence is missing for one or more microcatchments."
  )
}

ureh$ev_esc_baja_aj <-
  soft_runoff_evidence$ev_esc_baja_aj[
    evidence_index
  ]

ureh$ev_esc_media_aj <-
  soft_runoff_evidence$ev_esc_media_aj[
    evidence_index
  ]

ureh$ev_esc_alta_aj <-
  soft_runoff_evidence$ev_esc_alta_aj[
    evidence_index
  ]

if (anyNA(
  ureh[
    ,
    c(
      "ev_esc_baja_aj",
      "ev_esc_media_aj",
      "ev_esc_alta_aj"
    )
  ]
)) {
  stop(
    "Missing runoff evidence after matching by id."
  )
}

message(
  "Soft-runoff evidence matched to all microcatchments."
)


# ------------------------------------------------------------
# 8. Validate CPT dimension names
# ------------------------------------------------------------

validate_cpt_dimensions <- function(
    cpt,
    expected_dimensions,
    cpt_name
) {
  
  actual_dimensions <- names(
    dimnames(cpt)
  )
  
  if (!identical(
    actual_dimensions,
    expected_dimensions
  )) {
    stop(
      cpt_name,
      ": unexpected dimensions. Expected: ",
      paste(
        expected_dimensions,
        collapse = ", "
      ),
      "; found: ",
      paste(
        actual_dimensions,
        collapse = ", "
      ),
      "."
    )
  }
  
  invisible(TRUE)
}


validate_cpt_dimensions(
  cpt_infiltracion,
  c(
    "infiltracion",
    "cobertura",
    "pendiente",
    "precipitacion",
    "suelo"
  ),
  "CPT infiltration"
)

validate_cpt_dimensions(
  cpt_escorrentia,
  c(
    "escorrentia",
    "infiltracion",
    "precipitacion"
  ),
  "CPT runoff"
)

validate_cpt_dimensions(
  cpt_funcion_hidrologica,
  c(
    "funcion_hidrologica",
    "aporte",
    "escorrentia"
  ),
  "CPT hydrological function"
)

validate_cpt_dimensions(
  cpt_ureh,
  c(
    "ureh",
    "funcion_hidrologica",
    "cobertura"
  ),
  "CPT UREH"
)

message(
  "CPT dimensions validated for M0 inference."
)


# ------------------------------------------------------------
# 9. Probability-vector validator
# ------------------------------------------------------------

validate_probability_vector <- function(
    p,
    expected_names,
    object_name,
    tolerance = 1e-10
) {
  
  p <- as.numeric(p)
  
  if (length(p) != length(expected_names)) {
    stop(
      object_name,
      ": unexpected vector length."
    )
  }
  
  names(p) <- expected_names
  
  if (
    anyNA(p) ||
    any(!is.finite(p))
  ) {
    stop(
      object_name,
      ": missing or non-finite probabilities."
    )
  }
  
  if (
    any(p < -tolerance) ||
    any(p > 1 + tolerance)
  ) {
    stop(
      object_name,
      ": probabilities outside [0,1]."
    )
  }
  
  if (
    abs(sum(p) - 1) >
    tolerance
  ) {
    stop(
      object_name,
      ": probabilities do not sum to 1."
    )
  }
  
  p
}


# ------------------------------------------------------------
# 10. Normalized Shannon entropy
# ------------------------------------------------------------
#
# H = -sum(p_k log(p_k)) / log(K)
#
# K is fixed at the theoretical number of UREH states (6).
# Zero probabilities are removed only from the numerator because
# lim p*log(p) = 0 as p -> 0.
#
# ------------------------------------------------------------

normalized_shannon <- function(
    probabilities,
    n_states = 6L,
    tolerance = 1e-10
) {
  
  p <- as.numeric(
    probabilities
  )
  
  if (length(p) != n_states) {
    stop(
      "Entropy calculation received ",
      length(p),
      " probabilities; expected ",
      n_states,
      "."
    )
  }
  
  if (
    anyNA(p) ||
    any(!is.finite(p))
  ) {
    stop(
      "Entropy calculation received invalid probabilities."
    )
  }
  
  if (
    any(p < -tolerance) ||
    any(p > 1 + tolerance)
  ) {
    stop(
      "Entropy calculation received probabilities outside [0,1]."
    )
  }
  
  if (
    abs(sum(p) - 1) >
    tolerance
  ) {
    stop(
      "Entropy calculation received a distribution that does not sum to 1."
    )
  }
  
  p_positive <- p[
    p > 0
  ]
  
  entropy <-
    -sum(
      p_positive *
        log(p_positive)
    ) /
    log(n_states)
  
  if (
    entropy < -tolerance ||
    entropy > 1 + tolerance
  ) {
    stop(
      "Normalized entropy outside [0,1]."
    )
  }
  
  max(
    0,
    min(
      1,
      entropy
    )
  )
}


# ------------------------------------------------------------
# 11. M0 inference for one microcatchment
# ------------------------------------------------------------

infer_M0_unit <- function(
    suelo_i,
    cobertura_i,
    pendiente_i,
    precipitacion_i,
    aporte_i,
    ev_baja,
    ev_media,
    ev_alta
) {
  
  # ----------------------------------------------------------
  # 11.1 Normalize observed state labels
  # ----------------------------------------------------------
  
  suelo_i <-
    trimws(
      as.character(suelo_i)
    )
  
  cobertura_i <-
    trimws(
      as.character(cobertura_i)
    )
  
  pendiente_i <-
    trimws(
      as.character(pendiente_i)
    )
  
  precipitacion_i <-
    trimws(
      as.character(precipitacion_i)
    )
  
  aporte_i <-
    trimws(
      as.character(aporte_i)
    )
  
  
  # ----------------------------------------------------------
  # 11.2 Infiltration
  # ----------------------------------------------------------
  
  p_inf <- cpt_infiltracion[
    ,
    cobertura_i,
    pendiente_i,
    precipitacion_i,
    suelo_i,
    drop = TRUE
  ]
  
  p_inf <- validate_probability_vector(
    p_inf,
    estados$infiltracion,
    "Infiltration distribution"
  )
  
  
  # ----------------------------------------------------------
  # 11.3 Runoff implied by the Bayesian network
  # ----------------------------------------------------------
  
  p_esc_modelo <- setNames(
    rep(
      0,
      length(
        estados$escorrentia
      )
    ),
    estados$escorrentia
  )
  
  for (
    inf_i in estados$infiltracion
  ) {
    
    conditional_runoff <-
      cpt_escorrentia[
        ,
        inf_i,
        precipitacion_i,
        drop = TRUE
      ]
    
    conditional_runoff <-
      as.numeric(
        conditional_runoff
      )
    
    names(
      conditional_runoff
    ) <- estados$escorrentia
    
    p_esc_modelo <-
      p_esc_modelo +
      conditional_runoff *
      p_inf[inf_i]
  }
  
  p_esc_modelo <-
    validate_probability_vector(
      p_esc_modelo,
      estados$escorrentia,
      "Network runoff distribution"
    )
  
  
  # ----------------------------------------------------------
  # 11.4 Update runoff using regularized soft evidence
  # ----------------------------------------------------------
  
  likelihood_esc <- c(
    esc_baja =
      as.numeric(ev_baja),
    
    esc_media =
      as.numeric(ev_media),
    
    esc_alta =
      as.numeric(ev_alta)
  )
  
  if (
    anyNA(likelihood_esc) ||
    any(!is.finite(likelihood_esc)) ||
    any(likelihood_esc <= 0)
  ) {
    stop(
      "Invalid soft-runoff evidence supplied to M0."
    )
  }
  
  p_esc_final <-
    p_esc_modelo *
    likelihood_esc
  
  evidence_mass <-
    sum(
      p_esc_final
    )
  
  if (
    !is.finite(evidence_mass) ||
    evidence_mass <= 0
  ) {
    stop(
      "Network runoff distribution multiplied by soft evidence ",
      "produced zero or invalid probability mass."
    )
  }
  
  p_esc_final <-
    p_esc_final /
    evidence_mass
  
  p_esc_final <-
    validate_probability_vector(
      p_esc_final,
      estados$escorrentia,
      "Updated runoff distribution"
    )
  
  
  # ----------------------------------------------------------
  # 11.5 Hydrological function
  # ----------------------------------------------------------
  
  p_funcion <- setNames(
    rep(
      0,
      length(
        estados$funcion_hidrologica
      )
    ),
    estados$funcion_hidrologica
  )
  
  for (
    esc_i in estados$escorrentia
  ) {
    
    conditional_function <-
      cpt_funcion_hidrologica[
        ,
        aporte_i,
        esc_i,
        drop = TRUE
      ]
    
    conditional_function <-
      as.numeric(
        conditional_function
      )
    
    names(
      conditional_function
    ) <- estados$funcion_hidrologica
    
    p_funcion <-
      p_funcion +
      conditional_function *
      p_esc_final[esc_i]
  }
  
  p_funcion <-
    validate_probability_vector(
      p_funcion,
      estados$funcion_hidrologica,
      "Hydrological-function distribution"
    )
  
  
  # ----------------------------------------------------------
  # 11.6 UREH posterior
  # ----------------------------------------------------------
  
  p_ureh <- setNames(
    rep(
      0,
      length(
        estados$ureh
      )
    ),
    estados$ureh
  )
  
  for (
    function_i in
    estados$funcion_hidrologica
  ) {
    
    conditional_ureh <-
      cpt_ureh[
        ,
        function_i,
        cobertura_i,
        drop = TRUE
      ]
    
    conditional_ureh <-
      as.numeric(
        conditional_ureh
      )
    
    names(
      conditional_ureh
    ) <- estados$ureh
    
    p_ureh <-
      p_ureh +
      conditional_ureh *
      p_funcion[
        function_i
      ]
  }
  
  p_ureh <-
    validate_probability_vector(
      p_ureh,
      estados$ureh,
      "UREH posterior"
    )
  
  
  # ----------------------------------------------------------
  # 11.7 Return all internal distributions
  # ----------------------------------------------------------
  
  c(
    
    p_inf_baja =
      unname(
        p_inf["inf_baja"]
      ),
    
    p_inf_media =
      unname(
        p_inf["inf_media"]
      ),
    
    p_inf_alta =
      unname(
        p_inf["inf_alta"]
      ),
    
    p_esc_modelo_baja =
      unname(
        p_esc_modelo["esc_baja"]
      ),
    
    p_esc_modelo_media =
      unname(
        p_esc_modelo["esc_media"]
      ),
    
    p_esc_modelo_alta =
      unname(
        p_esc_modelo["esc_alta"]
      ),
    
    p_esc_final_baja =
      unname(
        p_esc_final["esc_baja"]
      ),
    
    p_esc_final_media =
      unname(
        p_esc_final["esc_media"]
      ),
    
    p_esc_final_alta =
      unname(
        p_esc_final["esc_alta"]
      ),
    
    p_func_aporte =
      unname(
        p_funcion["aporte"]
      ),
    
    p_func_transferencia =
      unname(
        p_funcion["transferencia"]
      ),
    
    p_func_recepcion =
      unname(
        p_funcion["recepcion"]
      ),
    
    aporte_conservada =
      unname(
        p_ureh[
          "aporte_conservada"
        ]
      ),
    
    aporte_degradada =
      unname(
        p_ureh[
          "aporte_degradada"
        ]
      ),
    
    transferencia_conservada =
      unname(
        p_ureh[
          "transferencia_conservada"
        ]
      ),
    
    transferencia_degradada =
      unname(
        p_ureh[
          "transferencia_degradada"
        ]
      ),
    
    receptora_conservada =
      unname(
        p_ureh[
          "receptora_conservada"
        ]
      ),
    
    receptora_degradada =
      unname(
        p_ureh[
          "receptora_degradada"
        ]
      )
  )
}


# ------------------------------------------------------------
# 12. Run M0 for all microcatchments
# ------------------------------------------------------------

message(
  "Running M0 inference for ",
  nrow(ureh),
  " microcatchments..."
)

M0_list <- vector(
  "list",
  nrow(ureh)
)

for (
  i in seq_len(
    nrow(ureh)
  )
) {
  
  M0_list[[i]] <- infer_M0_unit(
    
    suelo_i =
      ureh$suelo[i],
    
    cobertura_i =
      ureh$cobertura[i],
    
    pendiente_i =
      ureh$pendiente[i],
    
    precipitacion_i =
      ureh$precipitacion[i],
    
    aporte_i =
      ureh$aporte[i],
    
    ev_baja =
      ureh$ev_esc_baja_aj[i],
    
    ev_media =
      ureh$ev_esc_media_aj[i],
    
    ev_alta =
      ureh$ev_esc_alta_aj[i]
  )
}


# ------------------------------------------------------------
# 13. Assemble inference table
# ------------------------------------------------------------

M0_probabilities <- as.data.frame(
  do.call(
    rbind,
    M0_list
  )
)

M0_probabilities$id <-
  ureh$id

M0_probabilities <-
  M0_probabilities[
    ,
    c(
      "id",
      setdiff(
        names(M0_probabilities),
        "id"
      )
    )
  ]

if (
  nrow(M0_probabilities) !=
  n_expected_units
) {
  stop(
    "M0 inference produced an unexpected number of rows."
  )
}

if (
  anyDuplicated(
    M0_probabilities$id
  )
) {
  stop(
    "M0 inference produced duplicated ids."
  )
}

if (
  anyNA(
    M0_probabilities
  )
) {
  stop(
    "M0 inference produced missing values."
  )
}


# ------------------------------------------------------------
# 14. Validate all inferred distributions
# ------------------------------------------------------------

infiltration_columns <- c(
  "p_inf_baja",
  "p_inf_media",
  "p_inf_alta"
)

runoff_model_columns <- c(
  "p_esc_modelo_baja",
  "p_esc_modelo_media",
  "p_esc_modelo_alta"
)

runoff_final_columns <- c(
  "p_esc_final_baja",
  "p_esc_final_media",
  "p_esc_final_alta"
)

function_columns <- c(
  "p_func_aporte",
  "p_func_transferencia",
  "p_func_recepcion"
)

posterior_ureh_columns <-
  cols_ureh


validate_probability_rows <- function(
    x,
    columns,
    object_name,
    tolerance = 1e-10
) {
  
  m <- as.matrix(
    x[
      ,
      columns,
      drop = FALSE
    ]
  )
  
  storage.mode(m) <-
    "double"
  
  if (
    anyNA(m) ||
    any(!is.finite(m))
  ) {
    stop(
      object_name,
      " contains invalid probabilities."
    )
  }
  
  if (
    any(m < -tolerance) ||
    any(m > 1 + tolerance)
  ) {
    stop(
      object_name,
      " contains probabilities outside [0,1]."
    )
  }
  
  maximum_deviation <- max(
    abs(
      rowSums(m) - 1
    )
  )
  
  if (
    maximum_deviation >
    tolerance
  ) {
    stop(
      object_name,
      " does not sum to 1. Maximum deviation = ",
      maximum_deviation,
      "."
    )
  }
  
  invisible(
    maximum_deviation
  )
}


dev_inf <- validate_probability_rows(
  M0_probabilities,
  infiltration_columns,
  "M0 infiltration distributions"
)

dev_runoff_model <- validate_probability_rows(
  M0_probabilities,
  runoff_model_columns,
  "M0 network runoff distributions"
)

dev_runoff_final <- validate_probability_rows(
  M0_probabilities,
  runoff_final_columns,
  "M0 updated runoff distributions"
)

dev_function <- validate_probability_rows(
  M0_probabilities,
  function_columns,
  "M0 hydrological-function distributions"
)

dev_ureh <- validate_probability_rows(
  M0_probabilities,
  posterior_ureh_columns,
  "M0 UREH posterior distributions"
)

message(
  "All M0 probability distributions validated."
)


# ------------------------------------------------------------
# 15. Combine prepared data and M0 inference
# ------------------------------------------------------------

M0_complete <- ureh |>
  select(
    -ev_esc_baja_aj,
    -ev_esc_media_aj,
    -ev_esc_alta_aj
  ) |>
  left_join(
    M0_probabilities,
    by = "id"
  )

if (
  nrow(M0_complete) !=
  n_expected_units
) {
  stop(
    "M0 complete table has an unexpected number of rows."
  )
}

if (
  anyNA(
    M0_complete[
      ,
      posterior_ureh_columns
    ]
  )
) {
  stop(
    "UREH posterior missing after joining M0 results."
  )
}


# ------------------------------------------------------------
# 16. MAP classification
# ------------------------------------------------------------

posterior_matrix <- as.matrix(
  M0_complete[
    ,
    posterior_ureh_columns,
    drop = FALSE
  ]
)

storage.mode(
  posterior_matrix
) <- "double"


# Preserve historical tie-breaking:
# first state in the configured UREH order.

MAP_index <- max.col(
  posterior_matrix,
  ties.method = "first"
)

M0_complete$UREH_predicha <-
  posterior_ureh_columns[
    MAP_index
  ]

M0_complete$prob_UREH_predicha <-
  posterior_matrix[
    cbind(
      seq_len(
        nrow(posterior_matrix)
      ),
      MAP_index
    )
  ]


# ------------------------------------------------------------
# 17. Second-ranked UREH and classification margin
# ------------------------------------------------------------

probability_order <- t(
  apply(
    posterior_matrix,
    1,
    order,
    decreasing = TRUE
  )
)

second_index <-
  probability_order[
    ,
    2
  ]

M0_complete$UREH_segunda <-
  posterior_ureh_columns[
    second_index
  ]

M0_complete$prob_UREH_segunda <-
  posterior_matrix[
    cbind(
      seq_len(
        nrow(posterior_matrix)
      ),
      second_index
    )
  ]

M0_complete$margen_clasificacion <-
  M0_complete$prob_UREH_predicha -
  M0_complete$prob_UREH_segunda

if (
  any(
    M0_complete$margen_clasificacion <
    -1e-10
  )
) {
  stop(
    "Negative classification margins detected."
  )
}


# ------------------------------------------------------------
# 18. Normalized UREH entropy
# ------------------------------------------------------------

M0_complete$entropia_UREH <-
  apply(
    posterior_matrix,
    1,
    normalized_shannon,
    n_states =
      length(
        estados$ureh
      )
  )

if (
  anyNA(
    M0_complete$entropia_UREH
  )
) {
  stop(
    "Missing normalized entropy values detected."
  )
}

if (
  any(
    M0_complete$entropia_UREH <
    -1e-10
  ) ||
  any(
    M0_complete$entropia_UREH >
    1 + 1e-10
  )
) {
  stop(
    "Normalized UREH entropy outside [0,1]."
  )
}


# ------------------------------------------------------------
# 19. Classification output
# ------------------------------------------------------------

M0_classification <-
  M0_complete |>
  select(
    id,
    UREH_predicha,
    prob_UREH_predicha,
    UREH_segunda,
    prob_UREH_segunda,
    margen_clasificacion,
    entropia_UREH
  )


# ------------------------------------------------------------
# 20. Long-format UREH posterior output
# ------------------------------------------------------------

M0_posterior_long <-
  M0_complete |>
  select(
    id,
    all_of(
      posterior_ureh_columns
    )
  ) |>
  pivot_longer(
    cols =
      all_of(
        posterior_ureh_columns
      ),
    names_to =
      "UREH",
    values_to =
      "probabilidad"
  )

if (
  nrow(M0_posterior_long) !=
  n_expected_units *
  length(
    estados$ureh
  )
) {
  stop(
    "Unexpected number of rows in long-format UREH posterior."
  )
}


# ------------------------------------------------------------
# 21. Final classification validation
# ------------------------------------------------------------

if (
  anyNA(
    M0_classification
  )
) {
  stop(
    "M0 classification contains missing values."
  )
}

if (!all(
  M0_classification$UREH_predicha %in%
  estados$ureh
)) {
  stop(
    "Invalid MAP UREH classes detected."
  )
}

if (!all(
  M0_classification$UREH_segunda %in%
  estados$ureh
)) {
  stop(
    "Invalid second-ranked UREH classes detected."
  )
}

if (
  any(
    M0_classification$UREH_predicha ==
    M0_classification$UREH_segunda
  )
) {
  stop(
    "MAP and second-ranked UREH are identical for one or more units."
  )
}

message(
  "M0 MAP classifications and posterior-support metrics validated."
)


# ------------------------------------------------------------
# 22. Output paths
# ------------------------------------------------------------

complete_output_file <- here::here(
  "temp",
  "M0_UREH_resultados_completos.csv"
)

classification_output_file <- here::here(
  "results",
  "tables",
  "M0_UREH_clasificacion_dominante.csv"
)

posterior_output_file <- here::here(
  "results",
  "tables",
  "M0_UREH_probabilidades.csv"
)


# ------------------------------------------------------------
# 23. Write outputs
# ------------------------------------------------------------

write.csv(
  M0_complete,
  complete_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  M0_classification,
  classification_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  M0_posterior_long,
  posterior_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 24. Read-back validation
# ------------------------------------------------------------

complete_check <- read.csv(
  complete_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

classification_check <- read.csv(
  classification_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

posterior_check <- read.csv(
  posterior_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(complete_check) !=
  n_expected_units
) {
  stop(
    "M0 complete output failed read-back validation."
  )
}

if (
  nrow(classification_check) !=
  n_expected_units
) {
  stop(
    "M0 classification output failed read-back validation."
  )
}

if (
  nrow(posterior_check) !=
  n_expected_units *
  length(
    estados$ureh
  )
) {
  stop(
    "M0 posterior output failed read-back validation."
  )
}

if (
  anyDuplicated(
    classification_check$id
  )
) {
  stop(
    "Duplicate ids detected after classification read-back."
  )
}

posterior_sum_check <-
  posterior_check |>
  group_by(id) |>
  summarise(
    posterior_sum =
      sum(probabilidad),
    .groups = "drop"
  )

if (
  max(
    abs(
      posterior_sum_check$posterior_sum -
      1
    )
  ) >
  1e-10
) {
  stop(
    "Long-format UREH posterior failed read-back normalization check."
  )
}


# ------------------------------------------------------------
# 25. Console diagnostic summary
# ------------------------------------------------------------

class_counts <-
  M0_classification |>
  count(
    UREH_predicha,
    name = "n"
  ) |>
  mutate(
    percentage =
      100 *
      n /
      sum(n)
  ) |>
  arrange(
    match(
      UREH_predicha,
      estados$ureh
    )
  )

cat(
  "\n============================================\n"
)

cat(
  "M0 INFERENCE SUMMARY\n"
)

cat(
  "============================================\n"
)

cat(
  "Microcatchments:",
  nrow(M0_classification),
  "\n\n"
)

cat(
  "MAP UREH frequencies:\n"
)

print(
  class_counts,
  row.names = FALSE
)

cat(
  "\nPosterior support:\n"
)

cat(
  "Dominant probability - mean:",
  mean(
    M0_classification$prob_UREH_predicha
  ),
  "\n"
)

cat(
  "Dominant probability - median:",
  median(
    M0_classification$prob_UREH_predicha
  ),
  "\n"
)

cat(
  "Classification margin - mean:",
  mean(
    M0_classification$margen_clasificacion
  ),
  "\n"
)

cat(
  "Classification margin - median:",
  median(
    M0_classification$margen_clasificacion
  ),
  "\n"
)

cat(
  "Normalized entropy - mean:",
  mean(
    M0_classification$entropia_UREH
  ),
  "\n"
)

cat(
  "Normalized entropy - median:",
  median(
    M0_classification$entropia_UREH
  ),
  "\n"
)

cat(
  "Normalized entropy - range:",
  paste(
    range(
      M0_classification$entropia_UREH
    ),
    collapse = " - "
  ),
  "\n"
)

cat(
  "\nMaximum normalization deviations:\n"
)

cat(
  "Infiltration:",
  format(
    dev_inf,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Network runoff:",
  format(
    dev_runoff_model,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Updated runoff:",
  format(
    dev_runoff_final,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Hydrological function:",
  format(
    dev_function,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "UREH posterior:",
  format(
    dev_ureh,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 26. Completion
# ------------------------------------------------------------

message(
  "04_infer_M0.R completed successfully."
)

message(
  "Generated: temp/M0_UREH_resultados_completos.csv"
)

message(
  "Generated: results/tables/M0_UREH_clasificacion_dominante.csv"
)

message(
  "Generated: results/tables/M0_UREH_probabilidades.csv"
)