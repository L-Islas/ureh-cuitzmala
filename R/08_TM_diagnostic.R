# ============================================================
# 08. THORNTHWAITE-MATHER DIAGNOSTIC
# ============================================================
#
# Purpose
# -------
# Isolate the effect of the regularized Thornthwaite-Mather
# runoff evidence by comparing:
#
# A) Network only
#
#    P(runoff | infiltration, precipitation)
#
# B) Network + TM
#
#    P(runoff | infiltration, precipitation, E_TM)
#
# Both runoff distributions are then propagated through the
# same hydrological-function and UREH CPTs.
#
# The comparison is evaluated at two levels:
#
# 1. Runoff
# 2. Final UREH posterior
#
# Metrics
# -------
# - MAP stability
# - Jensen-Shannon distance
# - change in normalized Shannon entropy
#
# IMPORTANT
# ---------
# This diagnostic:
#
# - does not modify the DAG;
# - does not modify CPTs;
# - does not modify priors;
# - does not use hard evidence from Esc_sum;
# - does not define an additional main scenario;
# - must reproduce M0 when TM evidence is included.
#
# Inputs
# ------
# temp/ureh_prepared.csv
# temp/soft_runoff_evidence.csv
# temp/bayesian_network_model.rds
# temp/M0_UREH_resultados_completos.csv
#
# Outputs
# -------
# results/tables/tm_control_vs_M0.csv
# results/tables/tm_diagnostic_by_unit.csv
# results/tables/tm_effect_summary.csv
# results/tables/tm_runoff_MAP_transitions.csv
# results/tables/tm_UREH_MAP_transitions.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

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

M0_reference_file <- here::here(
  "temp",
  "M0_UREH_resultados_completos.csv"
)

required_files <- c(
  prepared_file,
  evidence_file,
  model_file,
  M0_reference_file
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
# 2. Read inputs
# ------------------------------------------------------------

ureh <- read.csv(
  prepared_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

soft_evidence <- read.csv(
  evidence_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

model <- readRDS(
  model_file
)

M0_reference <- read.csv(
  M0_reference_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Extract frozen model components
# ------------------------------------------------------------

required_model_objects <- c(
  "estados",
  "cpt_infiltracion",
  "cpt_escorrentia",
  "cpt_funcion_hidrologica",
  "cpt_ureh"
)

missing_model_objects <- setdiff(
  required_model_objects,
  names(model)
)

if (length(missing_model_objects) > 0L) {
  stop(
    "Incomplete Bayesian-network model object. Missing: ",
    paste(
      missing_model_objects,
      collapse = ", "
    )
  )
}

if (!identical(
  model$estados,
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

cols_ureh <-
  estados$ureh


# ------------------------------------------------------------
# 4. Validate unit counts and ids
# ------------------------------------------------------------

for (
  object_name in c(
    "ureh",
    "soft_evidence",
    "M0_reference"
  )
) {
  
  x <- get(
    object_name
  )
  
  if (
    nrow(x) !=
    n_expected_units
  ) {
    stop(
      object_name,
      " does not contain ",
      n_expected_units,
      " units."
    )
  }
  
  if (!("id" %in% names(x))) {
    stop(
      object_name,
      " does not contain an id column."
    )
  }
  
  if (
    anyNA(x$id) ||
    anyDuplicated(x$id)
  ) {
    stop(
      "Invalid ids in ",
      object_name,
      "."
    )
  }
}


# ------------------------------------------------------------
# 5. Match TM soft evidence by id
# ------------------------------------------------------------

required_evidence_columns <- c(
  "id",
  "ev_esc_baja_aj",
  "ev_esc_media_aj",
  "ev_esc_alta_aj"
)

missing_evidence_columns <- setdiff(
  required_evidence_columns,
  names(soft_evidence)
)

if (
  length(
    missing_evidence_columns
  ) > 0L
) {
  stop(
    "Soft-runoff evidence is missing columns: ",
    paste(
      missing_evidence_columns,
      collapse = ", "
    )
  )
}

evidence_index <- match(
  ureh$id,
  soft_evidence$id
)

if (
  anyNA(
    evidence_index
  )
) {
  stop(
    "TM soft evidence is missing for one or more microcatchments."
  )
}

ureh$ev_esc_baja_aj <-
  soft_evidence$ev_esc_baja_aj[
    evidence_index
  ]

ureh$ev_esc_media_aj <-
  soft_evidence$ev_esc_media_aj[
    evidence_index
  ]

ureh$ev_esc_alta_aj <-
  soft_evidence$ev_esc_alta_aj[
    evidence_index
  ]

evidence_columns <- c(
  "ev_esc_baja_aj",
  "ev_esc_media_aj",
  "ev_esc_alta_aj"
)

evidence_matrix <- as.matrix(
  ureh[
    ,
    evidence_columns,
    drop = FALSE
  ]
)

storage.mode(
  evidence_matrix
) <- "double"

if (
  anyNA(evidence_matrix) ||
  any(!is.finite(evidence_matrix)) ||
  any(evidence_matrix <= 0)
) {
  stop(
    "Invalid regularized TM evidence."
  )
}

if (
  max(
    abs(
      rowSums(
        evidence_matrix
      ) - 1
    )
  ) >
  1e-10
) {
  stop(
    "Regularized TM evidence does not sum to 1."
  )
}


# ------------------------------------------------------------
# 6. Validate observed states
# ------------------------------------------------------------

normalize_state <- function(x) {
  
  trimws(
    as.character(x)
  )
}


validate_state_variable <- function(
    x,
    allowed_states,
    variable_name
) {
  
  x <- normalize_state(x)
  
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
  
  if (
    length(
      invalid_states
    ) > 0L
  ) {
    stop(
      "Invalid states in ",
      variable_name,
      ": ",
      paste(
        invalid_states,
        collapse = ", "
      )
    )
  }
  
  invisible(TRUE)
}


validate_state_variable(
  ureh$suelo,
  estados$suelo,
  "suelo"
)

validate_state_variable(
  ureh$cobertura,
  estados$cobertura,
  "cobertura"
)

validate_state_variable(
  ureh$pendiente,
  estados$pendiente,
  "pendiente"
)

validate_state_variable(
  ureh$precipitacion,
  estados$precipitacion,
  "precipitacion"
)

validate_state_variable(
  ureh$aporte,
  estados$aporte,
  "aporte"
)


# ------------------------------------------------------------
# 7. Probability validator
# ------------------------------------------------------------

validate_probability_vector <- function(
    p,
    expected_states,
    object_name,
    tolerance = 1e-10
) {
  
  p <- as.numeric(p)
  
  if (
    length(p) !=
    length(expected_states)
  ) {
    stop(
      object_name,
      ": unexpected vector length."
    )
  }
  
  names(p) <-
    expected_states
  
  if (
    anyNA(p) ||
    any(!is.finite(p)) ||
    any(p < -tolerance) ||
    any(p > 1 + tolerance)
  ) {
    stop(
      object_name,
      ": invalid probabilities."
    )
  }
  
  if (
    abs(
      sum(p) - 1
    ) >
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
# 8. Normalized Shannon entropy
# ------------------------------------------------------------
#
# IMPORTANT:
#
# K is the theoretical number of states and remains fixed.
#
# Runoff:
#   K = 3
#
# UREH:
#   K = 6
#
# ------------------------------------------------------------

normalized_entropy <- function(
    p,
    n_states,
    tolerance = 1e-10
) {
  
  p <- as.numeric(p)
  
  if (
    length(p) !=
    n_states
  ) {
    stop(
      "Entropy received ",
      length(p),
      " probabilities but expected ",
      n_states,
      "."
    )
  }
  
  if (
    anyNA(p) ||
    any(!is.finite(p)) ||
    any(p < -tolerance) ||
    any(p > 1 + tolerance)
  ) {
    stop(
      "Invalid probabilities supplied to entropy."
    )
  }
  
  if (
    abs(
      sum(p) - 1
    ) >
    tolerance
  ) {
    stop(
      "Entropy received a distribution that does not sum to 1."
    )
  }
  
  positive <- p > 0
  
  H <-
    -sum(
      p[positive] *
        log(
          p[positive]
        )
    ) /
    log(
      n_states
    )
  
  max(
    0,
    min(
      1,
      H
    )
  )
}


# ------------------------------------------------------------
# 9. Jensen-Shannon distance
# ------------------------------------------------------------

jensen_shannon_distance <- function(
    p,
    q
) {
  
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  if (
    length(p) !=
    length(q)
  ) {
    stop(
      "JS distributions have different lengths."
    )
  }
  
  if (
    anyNA(p) ||
    anyNA(q) ||
    any(!is.finite(p)) ||
    any(!is.finite(q)) ||
    any(p < 0) ||
    any(q < 0) ||
    sum(p) <= 0 ||
    sum(q) <= 0
  ) {
    stop(
      "Invalid probabilities supplied to Jensen-Shannon distance."
    )
  }
  
  p <- p / sum(p)
  q <- q / sum(q)
  
  m <- 0.5 * (p + q)
  
  kl_base2 <- function(
    a,
    b
  ) {
    
    positive <- a > 0
    
    sum(
      a[positive] *
        log2(
          a[positive] /
            b[positive]
        )
    )
  }
  
  JS_divergence <-
    0.5 *
    kl_base2(
      p,
      m
    ) +
    0.5 *
    kl_base2(
      q,
      m
    )
  
  JS_divergence <- max(
    JS_divergence,
    0
  )
  
  distance <- sqrt(
    JS_divergence
  )
  
  if (
    !is.finite(distance) ||
    distance < -1e-12 ||
    distance > 1 + 1e-12
  ) {
    stop(
      "Jensen-Shannon distance outside [0,1]."
    )
  }
  
  distance
}


# ------------------------------------------------------------
# 10. Propagate a runoff distribution to function and UREH
# ------------------------------------------------------------

propagate_from_runoff <- function(
    p_runoff,
    contribution_state,
    cover_state
) {
  
  contribution_state <-
    normalize_state(
      contribution_state
    )
  
  cover_state <-
    normalize_state(
      cover_state
    )
  
  p_runoff <-
    validate_probability_vector(
      p_runoff,
      estados$escorrentia,
      "Runoff distribution"
    )
  
  
  # ----------------------------------------------------------
  # 10.1 Hydrological function
  # ----------------------------------------------------------
  
  p_function <- setNames(
    rep(
      0,
      length(
        estados$funcion_hidrologica
      )
    ),
    estados$funcion_hidrologica
  )
  
  for (
    runoff_state in
    estados$escorrentia
  ) {
    
    conditional_function <-
      cpt_funcion_hidrologica[
        ,
        contribution_state,
        runoff_state,
        drop = TRUE
      ]
    
    conditional_function <-
      setNames(
        as.numeric(
          conditional_function
        ),
        estados$funcion_hidrologica
      )
    
    p_function <-
      p_function +
      p_runoff[
        runoff_state
      ] *
      conditional_function
  }
  
  p_function <-
    validate_probability_vector(
      p_function,
      estados$funcion_hidrologica,
      "Hydrological-function distribution"
    )
  
  
  # ----------------------------------------------------------
  # 10.2 UREH
  # ----------------------------------------------------------
  
  p_UREH <- setNames(
    rep(
      0,
      length(
        estados$ureh
      )
    ),
    estados$ureh
  )
  
  for (
    function_state in
    estados$funcion_hidrologica
  ) {
    
    conditional_UREH <-
      cpt_ureh[
        ,
        function_state,
        cover_state,
        drop = TRUE
      ]
    
    conditional_UREH <-
      setNames(
        as.numeric(
          conditional_UREH
        ),
        estados$ureh
      )
    
    p_UREH <-
      p_UREH +
      p_function[
        function_state
      ] *
      conditional_UREH
  }
  
  p_UREH <-
    validate_probability_vector(
      p_UREH,
      estados$ureh,
      "UREH posterior"
    )
  
  
  list(
    function_distribution =
      p_function,
    
    UREH_distribution =
      p_UREH
  )
}


# ------------------------------------------------------------
# 11. Diagnostic for one microcatchment
# ------------------------------------------------------------

diagnose_unit <- function(i) {
  
  soil_state <-
    normalize_state(
      ureh$suelo[i]
    )
  
  cover_state <-
    normalize_state(
      ureh$cobertura[i]
    )
  
  slope_state <-
    normalize_state(
      ureh$pendiente[i]
    )
  
  precipitation_state <-
    normalize_state(
      ureh$precipitacion[i]
    )
  
  contribution_state <-
    normalize_state(
      ureh$aporte[i]
    )
  
  
  # ----------------------------------------------------------
  # 11.1 Infiltration
  # ----------------------------------------------------------
  
  p_infiltration <-
    cpt_infiltracion[
      ,
      cover_state,
      slope_state,
      precipitation_state,
      soil_state,
      drop = TRUE
    ]
  
  p_infiltration <-
    validate_probability_vector(
      p_infiltration,
      estados$infiltracion,
      "Infiltration distribution"
    )
  
  
  # ----------------------------------------------------------
  # 11.2 Network-only runoff
  # ----------------------------------------------------------
  
  p_runoff_network <- setNames(
    rep(
      0,
      length(
        estados$escorrentia
      )
    ),
    estados$escorrentia
  )
  
  for (
    infiltration_state in
    estados$infiltracion
  ) {
    
    conditional_runoff <-
      cpt_escorrentia[
        ,
        infiltration_state,
        precipitation_state,
        drop = TRUE
      ]
    
    conditional_runoff <-
      setNames(
        as.numeric(
          conditional_runoff
        ),
        estados$escorrentia
      )
    
    p_runoff_network <-
      p_runoff_network +
      p_infiltration[
        infiltration_state
      ] *
      conditional_runoff
  }
  
  p_runoff_network <-
    validate_probability_vector(
      p_runoff_network,
      estados$escorrentia,
      "Network-only runoff"
    )
  
  
  # ----------------------------------------------------------
  # 11.3 Network + TM runoff
  # ----------------------------------------------------------
  
  TM_likelihood <- setNames(
    c(
      ureh$ev_esc_baja_aj[i],
      ureh$ev_esc_media_aj[i],
      ureh$ev_esc_alta_aj[i]
    ),
    estados$escorrentia
  )
  
  p_runoff_TM <-
    p_runoff_network *
    TM_likelihood
  
  TM_mass <- sum(
    p_runoff_TM
  )
  
  if (
    !is.finite(TM_mass) ||
    TM_mass <= 0
  ) {
    stop(
      "Zero or invalid probability mass after TM update for id = ",
      ureh$id[i],
      "."
    )
  }
  
  p_runoff_TM <-
    p_runoff_TM /
    TM_mass
  
  p_runoff_TM <-
    validate_probability_vector(
      p_runoff_TM,
      estados$escorrentia,
      "TM-updated runoff"
    )
  
  
  # ----------------------------------------------------------
  # 11.4 Propagate both runoff distributions
  # ----------------------------------------------------------
  
  network_output <-
    propagate_from_runoff(
      p_runoff =
        p_runoff_network,
      
      contribution_state =
        contribution_state,
      
      cover_state =
        cover_state
    )
  
  TM_output <-
    propagate_from_runoff(
      p_runoff =
        p_runoff_TM,
      
      contribution_state =
        contribution_state,
      
      cover_state =
        cover_state
    )
  
  
  # ----------------------------------------------------------
  # 11.5 MAP states
  # ----------------------------------------------------------
  
  runoff_MAP_network <-
    names(
      p_runoff_network
    )[
      which.max(
        p_runoff_network
      )
    ]
  
  runoff_MAP_TM <-
    names(
      p_runoff_TM
    )[
      which.max(
        p_runoff_TM
      )
    ]
  
  UREH_MAP_network <-
    names(
      network_output$UREH_distribution
    )[
      which.max(
        network_output$UREH_distribution
      )
    ]
  
  UREH_MAP_TM <-
    names(
      TM_output$UREH_distribution
    )[
      which.max(
        TM_output$UREH_distribution
      )
    ]
  
  
  # ----------------------------------------------------------
  # 11.6 Entropy
  # ----------------------------------------------------------
  
  H_runoff_network <-
    normalized_entropy(
      p_runoff_network,
      n_states =
        length(
          estados$escorrentia
        )
    )
  
  H_runoff_TM <-
    normalized_entropy(
      p_runoff_TM,
      n_states =
        length(
          estados$escorrentia
        )
    )
  
  H_UREH_network <-
    normalized_entropy(
      network_output$UREH_distribution,
      n_states =
        length(
          estados$ureh
        )
    )
  
  H_UREH_TM <-
    normalized_entropy(
      TM_output$UREH_distribution,
      n_states =
        length(
          estados$ureh
        )
    )
  
  
  # ----------------------------------------------------------
  # 11.7 Return diagnostic information
  # ----------------------------------------------------------
  
  data.frame(
    
    id =
      ureh$id[i],
    
    
    # Network-only runoff
    
    runoff_network_low =
      unname(
        p_runoff_network[
          "esc_baja"
        ]
      ),
    
    runoff_network_medium =
      unname(
        p_runoff_network[
          "esc_media"
        ]
      ),
    
    runoff_network_high =
      unname(
        p_runoff_network[
          "esc_alta"
        ]
      ),
    
    
    # Network + TM runoff
    
    runoff_TM_low =
      unname(
        p_runoff_TM[
          "esc_baja"
        ]
      ),
    
    runoff_TM_medium =
      unname(
        p_runoff_TM[
          "esc_media"
        ]
      ),
    
    runoff_TM_high =
      unname(
        p_runoff_TM[
          "esc_alta"
        ]
      ),
    
    
    # Runoff diagnostic metrics
    
    MAP_runoff_network =
      runoff_MAP_network,
    
    MAP_runoff_TM =
      runoff_MAP_TM,
    
    JS_runoff =
      jensen_shannon_distance(
        p_runoff_network,
        p_runoff_TM
      ),
    
    H_runoff_network =
      H_runoff_network,
    
    H_runoff_TM =
      H_runoff_TM,
    
    delta_H_runoff =
      H_runoff_TM -
      H_runoff_network,
    
    
    # Network-only UREH posterior
    
    UREH_network_aporte_conservada =
      unname(
        network_output$UREH_distribution[
          "aporte_conservada"
        ]
      ),
    
    UREH_network_aporte_degradada =
      unname(
        network_output$UREH_distribution[
          "aporte_degradada"
        ]
      ),
    
    UREH_network_transferencia_conservada =
      unname(
        network_output$UREH_distribution[
          "transferencia_conservada"
        ]
      ),
    
    UREH_network_transferencia_degradada =
      unname(
        network_output$UREH_distribution[
          "transferencia_degradada"
        ]
      ),
    
    UREH_network_receptora_conservada =
      unname(
        network_output$UREH_distribution[
          "receptora_conservada"
        ]
      ),
    
    UREH_network_receptora_degradada =
      unname(
        network_output$UREH_distribution[
          "receptora_degradada"
        ]
      ),
    
    
    # Network + TM UREH posterior
    
    UREH_TM_aporte_conservada =
      unname(
        TM_output$UREH_distribution[
          "aporte_conservada"
        ]
      ),
    
    UREH_TM_aporte_degradada =
      unname(
        TM_output$UREH_distribution[
          "aporte_degradada"
        ]
      ),
    
    UREH_TM_transferencia_conservada =
      unname(
        TM_output$UREH_distribution[
          "transferencia_conservada"
        ]
      ),
    
    UREH_TM_transferencia_degradada =
      unname(
        TM_output$UREH_distribution[
          "transferencia_degradada"
        ]
      ),
    
    UREH_TM_receptora_conservada =
      unname(
        TM_output$UREH_distribution[
          "receptora_conservada"
        ]
      ),
    
    UREH_TM_receptora_degradada =
      unname(
        TM_output$UREH_distribution[
          "receptora_degradada"
        ]
      ),
    
    
    # UREH diagnostic metrics
    
    MAP_UREH_network =
      UREH_MAP_network,
    
    MAP_UREH_TM =
      UREH_MAP_TM,
    
    JS_UREH =
      jensen_shannon_distance(
        network_output$UREH_distribution,
        TM_output$UREH_distribution
      ),
    
    H_UREH_network =
      H_UREH_network,
    
    H_UREH_TM =
      H_UREH_TM,
    
    delta_H_UREH =
      H_UREH_TM -
      H_UREH_network,
    
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# 12. Run diagnostic for all microcatchments
# ------------------------------------------------------------

message("")
message(
  "============================================"
)
message(
  "THORNTHWAITE-MATHER DIAGNOSTIC"
)
message(
  "============================================"
)

diagnostic_list <- vector(
  "list",
  nrow(ureh)
)

for (
  i in seq_len(
    nrow(ureh)
  )
) {
  
  diagnostic_list[[i]] <-
    diagnose_unit(i)
  
  if (
    i %% 100L == 0L ||
    i == nrow(ureh)
  ) {
    message(
      "Processed: ",
      i,
      " / ",
      nrow(ureh)
    )
  }
}

TM_diagnostic <-
  bind_rows(
    diagnostic_list
  )


# ------------------------------------------------------------
# 13. Probability-distribution validation
# ------------------------------------------------------------

runoff_network_columns <- c(
  "runoff_network_low",
  "runoff_network_medium",
  "runoff_network_high"
)

runoff_TM_columns <- c(
  "runoff_TM_low",
  "runoff_TM_medium",
  "runoff_TM_high"
)

UREH_network_columns <- c(
  "UREH_network_aporte_conservada",
  "UREH_network_aporte_degradada",
  "UREH_network_transferencia_conservada",
  "UREH_network_transferencia_degradada",
  "UREH_network_receptora_conservada",
  "UREH_network_receptora_degradada"
)

UREH_TM_columns <- c(
  "UREH_TM_aporte_conservada",
  "UREH_TM_aporte_degradada",
  "UREH_TM_transferencia_conservada",
  "UREH_TM_transferencia_degradada",
  "UREH_TM_receptora_conservada",
  "UREH_TM_receptora_degradada"
)


validate_probability_table <- function(
    data,
    columns,
    object_name
) {
  
  m <- as.matrix(
    data[
      ,
      columns,
      drop = FALSE
    ]
  )
  
  storage.mode(m) <-
    "double"
  
  if (
    anyNA(m) ||
    any(!is.finite(m)) ||
    any(m < -1e-10) ||
    any(m > 1 + 1e-10)
  ) {
    stop(
      object_name,
      " contains invalid probabilities."
    )
  }
  
  maximum_deviation <- max(
    abs(
      rowSums(m) - 1
    )
  )
  
  if (
    maximum_deviation >
    1e-10
  ) {
    stop(
      object_name,
      " does not sum to 1."
    )
  }
  
  invisible(
    maximum_deviation
  )
}


dev_runoff_network <-
  validate_probability_table(
    TM_diagnostic,
    runoff_network_columns,
    "Network-only runoff"
  )

dev_runoff_TM <-
  validate_probability_table(
    TM_diagnostic,
    runoff_TM_columns,
    "Network + TM runoff"
  )

dev_UREH_network <-
  validate_probability_table(
    TM_diagnostic,
    UREH_network_columns,
    "Network-only UREH"
  )

dev_UREH_TM <-
  validate_probability_table(
    TM_diagnostic,
    UREH_TM_columns,
    "Network + TM UREH"
  )


# ------------------------------------------------------------
# 14. Mandatory control: Network + TM must reproduce M0
# ------------------------------------------------------------

M0_index <- match(
  TM_diagnostic$id,
  M0_reference$id
)

if (
  anyNA(
    M0_index
  )
) {
  stop(
    "Validated M0 reference does not contain all TM diagnostic ids."
  )
}

missing_M0_columns <- setdiff(
  c(
    cols_ureh,
    "UREH_predicha"
  ),
  names(
    M0_reference
  )
)

if (
  length(
    missing_M0_columns
  ) > 0L
) {
  stop(
    "M0 reference is missing columns: ",
    paste(
      missing_M0_columns,
      collapse = ", "
    )
  )
}


# Rename TM UREH matrix to canonical UREH names

TM_UREH_matrix <- as.matrix(
  TM_diagnostic[
    ,
    UREH_TM_columns,
    drop = FALSE
  ]
)

storage.mode(
  TM_UREH_matrix
) <- "double"

colnames(
  TM_UREH_matrix
) <- cols_ureh


M0_matrix <- as.matrix(
  M0_reference[
    M0_index,
    cols_ureh,
    drop = FALSE
  ]
)

storage.mode(
  M0_matrix
) <- "double"


maximum_M0_difference <- max(
  abs(
    TM_UREH_matrix -
      M0_matrix
  )
)

MAP_TM_equals_M0 <- mean(
  TM_diagnostic$MAP_UREH_TM ==
    M0_reference$UREH_predicha[
      M0_index
    ]
)

control_tolerance <- 1e-9

TM_control <- data.frame(
  
  n_units =
    nrow(
      TM_diagnostic
    ),
  
  MAP_agreement_TM_vs_M0 =
    MAP_TM_equals_M0,
  
  maximum_absolute_posterior_difference =
    maximum_M0_difference,
  
  tolerance =
    control_tolerance,
  
  control_passed =
    MAP_TM_equals_M0 == 1 &&
    maximum_M0_difference <
    control_tolerance,
  
  stringsAsFactors = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "CONTROL: NETWORK + TM vs M0\n"
)

cat(
  "============================================\n"
)

print(
  TM_control,
  row.names = FALSE
)

cat(
  "============================================\n\n"
)


if (
  !TM_control$control_passed
) {
  stop(
    paste0(
      "TM diagnostic control failed: ",
      "Network + TM does not reproduce M0."
    )
  )
}


# ------------------------------------------------------------
# 15. MAP changes
# ------------------------------------------------------------

TM_diagnostic <-
  TM_diagnostic |>
  mutate(
    
    runoff_MAP_changed =
      MAP_runoff_network !=
      MAP_runoff_TM,
    
    UREH_MAP_changed =
      MAP_UREH_network !=
      MAP_UREH_TM
  )


# ------------------------------------------------------------
# 16. Effect summary
# ------------------------------------------------------------

TM_effect_summary <- data.frame(
  
  level = c(
    "Runoff",
    "UREH"
  ),
  
  n_units = c(
    nrow(
      TM_diagnostic
    ),
    nrow(
      TM_diagnostic
    )
  ),
  
  n_MAP_stable = c(
    sum(
      !TM_diagnostic$runoff_MAP_changed
    ),
    sum(
      !TM_diagnostic$UREH_MAP_changed
    )
  ),
  
  percentage_MAP_stable = c(
    100 *
      mean(
        !TM_diagnostic$runoff_MAP_changed
      ),
    
    100 *
      mean(
        !TM_diagnostic$UREH_MAP_changed
      )
  ),
  
  n_MAP_changed = c(
    sum(
      TM_diagnostic$runoff_MAP_changed
    ),
    sum(
      TM_diagnostic$UREH_MAP_changed
    )
  ),
  
  percentage_MAP_changed = c(
    100 *
      mean(
        TM_diagnostic$runoff_MAP_changed
      ),
    
    100 *
      mean(
        TM_diagnostic$UREH_MAP_changed
      )
  ),
  
  median_JS = c(
    median(
      TM_diagnostic$JS_runoff
    ),
    
    median(
      TM_diagnostic$JS_UREH
    )
  ),
  
  mean_JS = c(
    mean(
      TM_diagnostic$JS_runoff
    ),
    
    mean(
      TM_diagnostic$JS_UREH
    )
  ),
  
  Q90_JS = c(
    unname(
      quantile(
        TM_diagnostic$JS_runoff,
        0.90
      )
    ),
    
    unname(
      quantile(
        TM_diagnostic$JS_UREH,
        0.90
      )
    )
  ),
  
  median_delta_H = c(
    median(
      TM_diagnostic$delta_H_runoff
    ),
    
    median(
      TM_diagnostic$delta_H_UREH
    )
  ),
  
  mean_delta_H = c(
    mean(
      TM_diagnostic$delta_H_runoff
    ),
    
    mean(
      TM_diagnostic$delta_H_UREH
    )
  ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 17. MAP transition tables
# ------------------------------------------------------------

TM_runoff_MAP_transitions <-
  TM_diagnostic |>
  count(
    MAP_runoff_network,
    MAP_runoff_TM,
    name = "n"
  ) |>
  mutate(
    
    percentage =
      100 *
      n /
      sum(n),
    
    changed =
      MAP_runoff_network !=
      MAP_runoff_TM
  ) |>
  arrange(
    desc(n)
  )


TM_UREH_MAP_transitions <-
  TM_diagnostic |>
  count(
    MAP_UREH_network,
    MAP_UREH_TM,
    name = "n"
  ) |>
  mutate(
    
    percentage =
      100 *
      n /
      sum(n),
    
    changed =
      MAP_UREH_network !=
      MAP_UREH_TM
  ) |>
  arrange(
    desc(n)
  )


# ------------------------------------------------------------
# 18. Output paths
# ------------------------------------------------------------

control_file <- here::here(
  "results",
  "tables",
  "tm_control_vs_M0.csv"
)

diagnostic_file <- here::here(
  "results",
  "tables",
  "tm_diagnostic_by_unit.csv"
)

summary_file <- here::here(
  "results",
  "tables",
  "tm_effect_summary.csv"
)

runoff_transitions_file <- here::here(
  "results",
  "tables",
  "tm_runoff_MAP_transitions.csv"
)

UREH_transitions_file <- here::here(
  "results",
  "tables",
  "tm_UREH_MAP_transitions.csv"
)


# ------------------------------------------------------------
# 19. Write outputs
# ------------------------------------------------------------

write.csv(
  TM_control,
  control_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  TM_diagnostic,
  diagnostic_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  TM_effect_summary,
  summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  TM_runoff_MAP_transitions,
  runoff_transitions_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  TM_UREH_MAP_transitions,
  UREH_transitions_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 20. Read-back validation
# ------------------------------------------------------------

control_check <- read.csv(
  control_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

diagnostic_check <- read.csv(
  diagnostic_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_check <- read.csv(
  summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

runoff_transitions_check <- read.csv(
  runoff_transitions_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

UREH_transitions_check <- read.csv(
  UREH_transitions_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(
    control_check
  ) != 1L
) {
  stop(
    "TM control failed read-back validation."
  )
}

if (
  nrow(
    diagnostic_check
  ) !=
  n_expected_units
) {
  stop(
    "TM diagnostic failed read-back validation."
  )
}

if (
  nrow(
    summary_check
  ) != 2L
) {
  stop(
    "TM effect summary failed read-back validation."
  )
}

if (
  sum(
    runoff_transitions_check$n
  ) !=
  n_expected_units
) {
  stop(
    "Runoff transition table failed read-back validation."
  )
}

if (
  sum(
    UREH_transitions_check$n
  ) !=
  n_expected_units
) {
  stop(
    "UREH transition table failed read-back validation."
  )
}


# ------------------------------------------------------------
# 21. Console results
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "THORNTHWAITE-MATHER EFFECT\n"
)

cat(
  "============================================\n"
)

print(
  TM_effect_summary,
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "RUNOFF MAP TRANSITIONS\n"
)

cat(
  "============================================\n"
)

print(
  TM_runoff_MAP_transitions,
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "UREH MAP TRANSITIONS\n"
)

cat(
  "============================================\n"
)

print(
  TM_UREH_MAP_transitions,
  row.names = FALSE
)


cat(
  "\nMaximum probability normalization deviations:\n"
)

cat(
  "Network runoff:",
  format(
    dev_runoff_network,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Network + TM runoff:",
  format(
    dev_runoff_TM,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Network UREH:",
  format(
    dev_UREH_network,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "Network + TM UREH:",
  format(
    dev_UREH_TM,
    scientific = TRUE
  ),
  "\n\n"
)


# ------------------------------------------------------------
# 22. Completion
# ------------------------------------------------------------

message(
  "08_TM_diagnostic.R completed successfully."
)

message(
  "Generated: results/tables/tm_control_vs_M0.csv"
)

message(
  "Generated: results/tables/tm_diagnostic_by_unit.csv"
)

message(
  "Generated: results/tables/tm_effect_summary.csv"
)

message(
  "Generated: results/tables/tm_runoff_MAP_transitions.csv"
)

message(
  "Generated: results/tables/tm_UREH_MAP_transitions.csv"
)