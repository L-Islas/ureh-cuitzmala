# ============================================================
# 06. SENSITIVITY ANALYSIS
# ============================================================
#
# Purpose
# -------
# Evaluate the sensitivity of the UREH regionalization to the
# removal of selected evidence while keeping fixed:
#
# - the DAG;
# - node definitions and states;
# - conditional probability tables;
# - empirical root-node priors;
# - the inference logic;
# - all evidence not explicitly removed.
#
# IMPORTANT
# ---------
# Removing evidence does NOT remove nodes and does NOT assign
# artificial replacement states.
#
# Unobserved root nodes are marginalized over all their states,
# weighted by their empirical priors.
#
# Main scenarios
# --------------
# M0: complete model (control).
#
# M1: remove land cover + soil + slope.
#
# M2: remove precipitation + soft runoff evidence.
#
# M3: remove upstream contribution.
#
# Diagnostic scenarios for M1
# ---------------------------
# "Sin cobertura":
#   remove land cover only.
#
# "Sin suelo + pendiente":
#   remove soil and slope while retaining land cover.
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
# results/tables/sensitivity_engine_control.csv
# results/tables/sensitivity_scenarios.csv
# results/tables/sensitivity_posteriors.csv
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
    ),
    "\nRun scripts 01-04 first."
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

M0_reference <- read.csv(
  M0_reference_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Validate number and identity of units
# ------------------------------------------------------------

if (nrow(ureh) != n_expected_units) {
  stop(
    "Prepared dataset must contain ",
    n_expected_units,
    " microcatchments."
  )
}

if (nrow(soft_runoff_evidence) != n_expected_units) {
  stop(
    "Soft-runoff evidence must contain ",
    n_expected_units,
    " microcatchments."
  )
}

if (nrow(M0_reference) != n_expected_units) {
  stop(
    "M0 reference must contain ",
    n_expected_units,
    " microcatchments."
  )
}

for (
  x_name in c(
    "ureh",
    "soft_runoff_evidence",
    "M0_reference"
  )
) {
  
  x <- get(x_name)
  
  if (!("id" %in% names(x))) {
    stop(
      x_name,
      " does not contain an id field."
    )
  }
  
  if (
    anyNA(x$id) ||
    anyDuplicated(x$id)
  ) {
    stop(
      "Invalid ids in ",
      x_name,
      "."
    )
  }
  
  if (!identical(
    sort(as.integer(x$id)),
    seq_len(n_expected_units)
  )) {
    stop(
      "Unexpected ids in ",
      x_name,
      "."
    )
  }
}


# ------------------------------------------------------------
# 4. Validate model object
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

if (length(cols_ureh) != 6L) {
  stop(
    "Sensitivity analysis expects six UREH states."
  )
}


# ------------------------------------------------------------
# 5. Prepare root-node priors
# ------------------------------------------------------------

prepare_prior <- function(
    x,
    expected_states,
    node_name
) {
  
  values <- as.numeric(x)
  
  state_names <- names(x)
  
  if (
    is.null(state_names) &&
    !is.null(dimnames(x))
  ) {
    state_names <- dimnames(x)[[1]]
  }
  
  names(values) <-
    state_names
  
  if (is.null(
    names(values)
  )) {
    stop(
      "Prior for ",
      node_name,
      " has no state names."
    )
  }
  
  missing_states <- setdiff(
    expected_states,
    names(values)
  )
  
  if (length(missing_states) > 0L) {
    stop(
      "Prior for ",
      node_name,
      " is missing states: ",
      paste(
        missing_states,
        collapse = ", "
      )
    )
  }
  
  values <-
    values[
      expected_states
    ]
  
  if (
    anyNA(values) ||
    any(!is.finite(values)) ||
    any(values < 0) ||
    sum(values) <= 0
  ) {
    stop(
      "Invalid prior for ",
      node_name,
      "."
    )
  }
  
  values <-
    values /
    sum(values)
  
  if (
    abs(
      sum(values) - 1
    ) >
    1e-10
  ) {
    stop(
      "Prior for ",
      node_name,
      " does not sum to 1."
    )
  }
  
  values
}


priors <- list(
  
  suelo = prepare_prior(
    model$priors$suelo,
    estados$suelo,
    "suelo"
  ),
  
  cobertura = prepare_prior(
    model$priors$cobertura,
    estados$cobertura,
    "cobertura"
  ),
  
  pendiente = prepare_prior(
    model$priors$pendiente,
    estados$pendiente,
    "pendiente"
  ),
  
  precipitacion = prepare_prior(
    model$priors$precipitacion,
    estados$precipitacion,
    "precipitacion"
  ),
  
  aporte = prepare_prior(
    model$priors$aporte,
    estados$aporte,
    "aporte"
  )
)

message(
  "Root-node priors loaded and validated."
)


# ------------------------------------------------------------
# 6. Match soft runoff evidence by id
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

evidence_index <- match(
  ureh$id,
  soft_runoff_evidence$id
)

if (anyNA(evidence_index)) {
  stop(
    "Soft-runoff evidence is missing for one or more units."
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
    "Invalid regularized runoff evidence."
  )
}

if (
  max(
    abs(
      rowSums(evidence_matrix) - 1
    )
  ) >
  1e-10
) {
  stop(
    "Regularized runoff evidence does not sum to 1."
  )
}

message(
  "Soft-runoff evidence matched and validated."
)


# ------------------------------------------------------------
# 7. Validate observed root-node states
# ------------------------------------------------------------

normalize_state <- function(x) {
  trimws(
    as.character(x)
  )
}


validate_observed_variable <- function(
    x,
    allowed_states,
    variable_name
) {
  
  x <- normalize_state(x)
  
  if (anyNA(x)) {
    stop(
      "Variable ",
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
      "Unrecognized states in ",
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


validate_observed_variable(
  ureh$suelo,
  estados$suelo,
  "suelo"
)

validate_observed_variable(
  ureh$cobertura,
  estados$cobertura,
  "cobertura"
)

validate_observed_variable(
  ureh$pendiente,
  estados$pendiente,
  "pendiente"
)

validate_observed_variable(
  ureh$precipitacion,
  estados$precipitacion,
  "precipitacion"
)

validate_observed_variable(
  ureh$aporte,
  estados$aporte,
  "aporte"
)

message(
  "Observed root-node states validated."
)


# ------------------------------------------------------------
# 8. Validate CPT dimension structure
# ------------------------------------------------------------

validate_cpt_dimensions <- function(
    cpt,
    expected_dimensions,
    cpt_name
) {
  
  actual_dimensions <-
    names(
      dimnames(cpt)
    )
  
  if (!identical(
    actual_dimensions,
    expected_dimensions
  )) {
    stop(
      cpt_name,
      ": unexpected dimensions."
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


# ------------------------------------------------------------
# 9. Root-node support function
# ------------------------------------------------------------
#
# Observed root:
#   weight 1 for the observed state.
#
# Unobserved root:
#   all states retained and weighted by their empirical prior.
#
# ------------------------------------------------------------

root_support <- function(
    observed,
    node_states,
    node_prior
) {
  
  if (!is.null(observed)) {
    
    observed <-
      normalize_state(
        observed
      )
    
    if (
      length(observed) != 1L ||
      is.na(observed) ||
      !(observed %in% node_states)
    ) {
      stop(
        "Invalid observed root-node state: ",
        observed
      )
    }
    
    result <- setNames(
      1,
      observed
    )
    
    return(result)
  }
  
  result <-
    node_prior[
      node_states
    ]
  
  if (
    anyNA(result) ||
    abs(sum(result) - 1) >
    1e-10
  ) {
    stop(
      "Invalid marginal support for an unobserved root node."
    )
  }
  
  result
}


# ------------------------------------------------------------
# 10. Normalized Shannon entropy
# ------------------------------------------------------------
#
# K remains fixed at six UREH states.
#
# ------------------------------------------------------------

normalized_shannon <- function(
    probabilities,
    n_states =
      length(estados$ureh)
) {
  
  p <-
    as.numeric(
      probabilities
    )
  
  if (
    length(p) != n_states ||
    anyNA(p) ||
    any(!is.finite(p)) ||
    any(p < -1e-10) ||
    any(p > 1 + 1e-10) ||
    abs(sum(p) - 1) > 1e-10
  ) {
    stop(
      "Invalid probability distribution supplied to entropy calculation."
    )
  }
  
  p_positive <-
    p[
      p > 0
    ]
  
  entropy <-
    -sum(
      p_positive *
        log(p_positive)
    ) /
    log(n_states)
  
  max(
    0,
    min(
      1,
      entropy
    )
  )
}


# ------------------------------------------------------------
# 11. Partial-evidence inference engine
# ------------------------------------------------------------
#
# For each unobserved root node:
#
#   P(root state) = empirical prior.
#
# All root-state combinations are propagated through:
#
# roots
#   ↓
# infiltration
#   ↓
# runoff
#   ↓
# hydrological function
#   ↓
# UREH
#
# If runoff soft evidence is retained, it acts as a likelihood
# on the runoff state.
#
# Normalization occurs after all configurations have been
# accumulated.
#
# ------------------------------------------------------------

infer_UREH_partial <- function(
    suelo_obs = NULL,
    cobertura_obs = NULL,
    pendiente_obs = NULL,
    precipitacion_obs = NULL,
    aporte_obs = NULL,
    runoff_evidence = NULL
) {
  
  # ----------------------------------------------------------
  # 11.1 Root supports
  # ----------------------------------------------------------
  
  support_soil <- root_support(
    suelo_obs,
    estados$suelo,
    priors$suelo
  )
  
  support_cover <- root_support(
    cobertura_obs,
    estados$cobertura,
    priors$cobertura
  )
  
  support_slope <- root_support(
    pendiente_obs,
    estados$pendiente,
    priors$pendiente
  )
  
  support_precip <- root_support(
    precipitacion_obs,
    estados$precipitacion,
    priors$precipitacion
  )
  
  support_contribution <- root_support(
    aporte_obs,
    estados$aporte,
    priors$aporte
  )
  
  
  # ----------------------------------------------------------
  # 11.2 Runoff likelihood
  # ----------------------------------------------------------
  
  if (is.null(
    runoff_evidence
  )) {
    
    runoff_likelihood <- setNames(
      rep(
        1,
        length(
          estados$escorrentia
        )
      ),
      estados$escorrentia
    )
    
  } else {
    
    runoff_evidence <-
      as.numeric(
        runoff_evidence
      )
    
    if (
      length(runoff_evidence) !=
      length(estados$escorrentia) ||
      anyNA(runoff_evidence) ||
      any(!is.finite(runoff_evidence)) ||
      any(runoff_evidence < 0) ||
      sum(runoff_evidence) <= 0
    ) {
      stop(
        "Invalid soft-runoff evidence."
      )
    }
    
    runoff_likelihood <- setNames(
      runoff_evidence,
      estados$escorrentia
    )
  }
  
  
  # ----------------------------------------------------------
  # 11.3 All compatible root-state configurations
  # ----------------------------------------------------------
  
  root_combinations <- expand.grid(
    
    suelo =
      names(
        support_soil
      ),
    
    cobertura =
      names(
        support_cover
      ),
    
    pendiente =
      names(
        support_slope
      ),
    
    precipitacion =
      names(
        support_precip
      ),
    
    aporte =
      names(
        support_contribution
      ),
    
    stringsAsFactors = FALSE
  )
  
  root_combinations$weight <-
    
    support_soil[
      root_combinations$suelo
    ] *
    
    support_cover[
      root_combinations$cobertura
    ] *
    
    support_slope[
      root_combinations$pendiente
    ] *
    
    support_precip[
      root_combinations$precipitacion
    ] *
    
    support_contribution[
      root_combinations$aporte
    ]
  
  
  # ----------------------------------------------------------
  # 11.4 UREH posterior numerator
  # ----------------------------------------------------------
  
  UREH_numerator <- setNames(
    rep(
      0,
      length(
        estados$ureh
      )
    ),
    estados$ureh
  )
  
  
  # ----------------------------------------------------------
  # 11.5 Exact propagation
  # ----------------------------------------------------------
  
  for (
    j in seq_len(
      nrow(
        root_combinations
      )
    )
  ) {
    
    cfg <-
      root_combinations[
        j,
        ,
        drop = FALSE
      ]
    
    root_weight <-
      as.numeric(
        cfg$weight
      )
    
    if (root_weight == 0) {
      next
    }
    
    
    # Infiltration distribution
    
    p_infiltration <-
      cpt_infiltracion[
        ,
        cfg$cobertura,
        cfg$pendiente,
        cfg$precipitacion,
        cfg$suelo,
        drop = TRUE
      ]
    
    p_infiltration <- setNames(
      as.numeric(
        p_infiltration
      ),
      estados$infiltracion
    )
    
    
    # Marginalize infiltration
    
    for (
      infiltration_state in
      estados$infiltracion
    ) {
      
      infiltration_weight <-
        root_weight *
        p_infiltration[
          infiltration_state
        ]
      
      if (
        infiltration_weight == 0
      ) {
        next
      }
      
      
      # Runoff distribution
      
      p_runoff <-
        cpt_escorrentia[
          ,
          infiltration_state,
          cfg$precipitacion,
          drop = TRUE
        ]
      
      p_runoff <- setNames(
        as.numeric(
          p_runoff
        ),
        estados$escorrentia
      )
      
      
      # Marginalize runoff
      
      for (
        runoff_state in
        estados$escorrentia
      ) {
        
        runoff_weight <-
          
          infiltration_weight *
          
          p_runoff[
            runoff_state
          ] *
          
          runoff_likelihood[
            runoff_state
          ]
        
        if (
          runoff_weight == 0
        ) {
          next
        }
        
        
        # Hydrological-function distribution
        
        p_function <-
          cpt_funcion_hidrologica[
            ,
            cfg$aporte,
            runoff_state,
            drop = TRUE
          ]
        
        p_function <- setNames(
          as.numeric(
            p_function
          ),
          estados$funcion_hidrologica
        )
        
        
        # Marginalize hydrological function
        
        for (
          function_state in
          estados$funcion_hidrologica
        ) {
          
          function_weight <-
            
            runoff_weight *
            
            p_function[
              function_state
            ]
          
          if (
            function_weight == 0
          ) {
            next
          }
          
          
          # Conditional UREH distribution
          
          p_UREH <-
            cpt_ureh[
              ,
              function_state,
              cfg$cobertura,
              drop = TRUE
            ]
          
          p_UREH <- setNames(
            as.numeric(
              p_UREH
            ),
            estados$ureh
          )
          
          
          UREH_numerator <-
            UREH_numerator +
            function_weight *
            p_UREH
        }
      }
    }
  }
  
  
  # ----------------------------------------------------------
  # 11.6 Final normalization
  # ----------------------------------------------------------
  
  total_mass <-
    sum(
      UREH_numerator
    )
  
  if (
    !is.finite(total_mass) ||
    total_mass <= 0
  ) {
    stop(
      "Marginalization produced zero or non-finite probability mass."
    )
  }
  
  posterior <-
    UREH_numerator /
    total_mass
  
  if (
    anyNA(posterior) ||
    any(!is.finite(posterior)) ||
    any(posterior < -1e-10) ||
    any(posterior > 1 + 1e-10) ||
    abs(sum(posterior) - 1) >
    1e-10
  ) {
    stop(
      "Invalid UREH posterior after marginalization."
    )
  }
  
  posterior
}


# ------------------------------------------------------------
# 12. Define evidence-removal scenarios
# ------------------------------------------------------------
#
# TRUE:
#   evidence remains observed.
#
# FALSE:
#   evidence is removed and the root node is marginalized.
#
# use_runoff_evidence = FALSE:
#   the Thornthwaite-Mather soft runoff evidence is not used.
#
# ------------------------------------------------------------

scenarios <- list(
  
  M0 = list(
    
    scenario_type =
      "control",
    
    evidence_removed =
      "Ninguna",
    
    question =
      "Referencia completa",
    
    use_soil =
      TRUE,
    
    use_cover =
      TRUE,
    
    use_slope =
      TRUE,
    
    use_precipitation =
      TRUE,
    
    use_contribution =
      TRUE,
    
    use_runoff_evidence =
      TRUE
  ),
  
  
  M1 = list(
    
    scenario_type =
      "principal",
    
    evidence_removed =
      "Cobertura + suelo + pendiente",
    
    question =
      "Dependencia de informacion estructural local",
    
    use_soil =
      FALSE,
    
    use_cover =
      FALSE,
    
    use_slope =
      FALSE,
    
    use_precipitation =
      TRUE,
    
    use_contribution =
      TRUE,
    
    use_runoff_evidence =
      TRUE
  ),
  
  
  M2 = list(
    
    scenario_type =
      "principal",
    
    evidence_removed =
      "Precipitacion + evidencia blanda de escorrentia",
    
    question =
      "Dependencia de informacion climatica/hidrologica local",
    
    use_soil =
      TRUE,
    
    use_cover =
      TRUE,
    
    use_slope =
      TRUE,
    
    use_precipitation =
      FALSE,
    
    use_contribution =
      TRUE,
    
    use_runoff_evidence =
      FALSE
  ),
  
  
  M3 = list(
    
    scenario_type =
      "principal",
    
    evidence_removed =
      "Aporte acumulado de aguas arriba",
    
    question =
      "Dependencia de informacion asociada con el aporte aguas arriba",
    
    use_soil =
      TRUE,
    
    use_cover =
      TRUE,
    
    use_slope =
      TRUE,
    
    use_precipitation =
      TRUE,
    
    use_contribution =
      FALSE,
    
    use_runoff_evidence =
      TRUE
  ),
  
  
  "Sin cobertura" = list(
    
    scenario_type =
      "diagnostico",
    
    evidence_removed =
      "Cobertura",
    
    question =
      "Distinguir el efecto de cobertura dentro de M1",
    
    use_soil =
      TRUE,
    
    use_cover =
      FALSE,
    
    use_slope =
      TRUE,
    
    use_precipitation =
      TRUE,
    
    use_contribution =
      TRUE,
    
    use_runoff_evidence =
      TRUE
  ),
  
  
  "Sin suelo + pendiente" = list(
    
    scenario_type =
      "diagnostico",
    
    evidence_removed =
      "Suelo + pendiente",
    
    question =
      "Distinguir el efecto indirecto de suelo y pendiente dentro de M1",
    
    use_soil =
      FALSE,
    
    use_cover =
      TRUE,
    
    use_slope =
      FALSE,
    
    use_precipitation =
      TRUE,
    
    use_contribution =
      TRUE,
    
    use_runoff_evidence =
      TRUE
  )
)


# ------------------------------------------------------------
# 13. Run one scenario
# ------------------------------------------------------------

run_scenario <- function(
    scenario_name,
    config
) {
  
  message("")
  message(
    "============================================"
  )
  message(
    "SCENARIO: ",
    scenario_name
  )
  message(
    "Evidence removed: ",
    config$evidence_removed
  )
  message(
    "============================================"
  )
  
  posterior_matrix <- matrix(
    NA_real_,
    nrow =
      nrow(ureh),
    ncol =
      length(cols_ureh),
    dimnames =
      list(
        NULL,
        cols_ureh
      )
  )
  
  
  for (
    i in seq_len(
      nrow(ureh)
    )
  ) {
    
    runoff_evidence_i <- NULL
    
    if (
      config$use_runoff_evidence
    ) {
      
      runoff_evidence_i <- c(
        
        ureh$ev_esc_baja_aj[i],
        
        ureh$ev_esc_media_aj[i],
        
        ureh$ev_esc_alta_aj[i]
      )
    }
    
    
    posterior_matrix[
      i,
    ] <- infer_UREH_partial(
      
      suelo_obs =
        if (
          config$use_soil
        ) {
          ureh$suelo[i]
        } else {
          NULL
        },
      
      cobertura_obs =
        if (
          config$use_cover
        ) {
          ureh$cobertura[i]
        } else {
          NULL
        },
      
      pendiente_obs =
        if (
          config$use_slope
        ) {
          ureh$pendiente[i]
        } else {
          NULL
        },
      
      precipitacion_obs =
        if (
          config$use_precipitation
        ) {
          ureh$precipitacion[i]
        } else {
          NULL
        },
      
      aporte_obs =
        if (
          config$use_contribution
        ) {
          ureh$aporte[i]
        } else {
          NULL
        },
      
      runoff_evidence =
        runoff_evidence_i
    )
    
    
    if (
      i %% 100L == 0L ||
      i == nrow(ureh)
    ) {
      message(
        "  Processed: ",
        i,
        " / ",
        nrow(ureh)
      )
    }
  }
  
  
  # ----------------------------------------------------------
  # 13.1 Validate posterior matrix
  # ----------------------------------------------------------
  
  if (
    anyNA(posterior_matrix) ||
    any(!is.finite(posterior_matrix))
  ) {
    stop(
      "Invalid posterior probabilities in scenario ",
      scenario_name,
      "."
    )
  }
  
  if (
    any(
      posterior_matrix <
      -1e-10
    ) ||
    any(
      posterior_matrix >
      1 + 1e-10
    )
  ) {
    stop(
      "Posterior probabilities outside [0,1] in scenario ",
      scenario_name,
      "."
    )
  }
  
  maximum_deviation <-
    max(
      abs(
        rowSums(
          posterior_matrix
        ) - 1
      )
    )
  
  if (
    maximum_deviation >
    1e-10
  ) {
    stop(
      "Posterior normalization failed in scenario ",
      scenario_name,
      "."
    )
  }
  
  
  # ----------------------------------------------------------
  # 13.2 Create scenario output
  # ----------------------------------------------------------
  
  output <-
    as.data.frame(
      posterior_matrix,
      stringsAsFactors = FALSE
    )
  
  output$id <-
    ureh$id
  
  output$escenario <-
    scenario_name
  
  output$tipo_escenario <-
    config$scenario_type
  
  output$evidencia_retirada <-
    config$evidence_removed
  
  output$pregunta <-
    config$question
  
  
  # MAP UREH
  
  MAP_index <- max.col(
    posterior_matrix,
    ties.method = "first"
  )
  
  output$UREH_predicha <-
    cols_ureh[
      MAP_index
    ]
  
  output$prob_UREH_predicha <-
    posterior_matrix[
      cbind(
        seq_len(
          nrow(
            posterior_matrix
          )
        ),
        MAP_index
      )
    ]
  
  
  # Normalized entropy
  
  output$entropia_UREH <-
    apply(
      posterior_matrix,
      1,
      normalized_shannon
    )
  
  
  output <-
    output[
      ,
      c(
        "id",
        "escenario",
        "tipo_escenario",
        "evidencia_retirada",
        "pregunta",
        cols_ureh,
        "UREH_predicha",
        "prob_UREH_predicha",
        "entropia_UREH"
      )
    ]
  
  
  output
}


# ------------------------------------------------------------
# 14. Mandatory engine control: M0 vs validated M0
# ------------------------------------------------------------
#
# The partial-evidence engine must reproduce the validated M0
# posterior before sensitivity scenarios are allowed to run.
#
# ------------------------------------------------------------

M0_from_engine <- run_scenario(
  "M0",
  scenarios[["M0"]]
)

reference_index <- match(
  M0_from_engine$id,
  M0_reference$id
)

if (anyNA(
  reference_index
)) {
  stop(
    "Validated M0 reference does not contain all sensitivity-engine ids."
  )
}

missing_M0_columns <- setdiff(
  c(
    cols_ureh,
    "UREH_predicha"
  ),
  names(M0_reference)
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

engine_matrix <-
  as.matrix(
    M0_from_engine[
      ,
      cols_ureh,
      drop = FALSE
    ]
  )

reference_matrix <-
  as.matrix(
    M0_reference[
      reference_index,
      cols_ureh,
      drop = FALSE
    ]
  )

storage.mode(
  engine_matrix
) <- "double"

storage.mode(
  reference_matrix
) <- "double"

maximum_engine_difference <-
  max(
    abs(
      engine_matrix -
        reference_matrix
    )
  )

MAP_agreement <-
  mean(
    M0_from_engine$UREH_predicha ==
      M0_reference$UREH_predicha[
        reference_index
      ]
  )

engine_tolerance <-
  1e-9

engine_control <- data.frame(
  
  n_units =
    nrow(
      M0_from_engine
    ),
  
  MAP_agreement =
    MAP_agreement,
  
  maximum_absolute_posterior_difference =
    maximum_engine_difference,
  
  tolerance =
    engine_tolerance,
  
  stringsAsFactors = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "SENSITIVITY ENGINE CONTROL\n"
)

cat(
  "============================================\n"
)

print(
  engine_control,
  row.names = FALSE
)

cat(
  "============================================\n\n"
)


if (
  MAP_agreement != 1 ||
  maximum_engine_difference >=
  engine_tolerance
) {
  
  stop(
    "Sensitivity-engine control failed. ",
    "M1-M3 and diagnostic scenarios will not be executed."
  )
}

message(
  "ENGINE CONTROL PASSED: partial-evidence inference reproduces M0."
)


# ------------------------------------------------------------
# 15. Run remaining scenarios
# ------------------------------------------------------------

scenario_results <- list(
  M0 =
    M0_from_engine
)

for (
  scenario_name in
  setdiff(
    names(scenarios),
    "M0"
  )
) {
  
  scenario_results[[scenario_name]] <-
    run_scenario(
      scenario_name,
      scenarios[[scenario_name]]
    )
}


sensitivity_posteriors <-
  bind_rows(
    scenario_results
  )


# ------------------------------------------------------------
# 16. Final posterior validations
# ------------------------------------------------------------

n_scenarios <-
  length(
    scenarios
  )

expected_rows <-
  n_expected_units *
  n_scenarios

if (
  nrow(
    sensitivity_posteriors
  ) != expected_rows
) {
  stop(
    "Unexpected number of rows in sensitivity output."
  )
}

if (
  anyDuplicated(
    sensitivity_posteriors[
      ,
      c(
        "id",
        "escenario"
      )
    ]
  )
) {
  stop(
    "Duplicated id-scenario combinations detected."
  )
}

if (
  anyNA(
    sensitivity_posteriors[
      ,
      cols_ureh,
      drop = FALSE
    ]
  )
) {
  stop(
    "Missing UREH posterior probabilities in sensitivity output."
  )
}

all_posterior_sums <-
  rowSums(
    sensitivity_posteriors[
      ,
      cols_ureh,
      drop = FALSE
    ]
  )

maximum_posterior_deviation <-
  max(
    abs(
      all_posterior_sums - 1
    )
  )

if (
  maximum_posterior_deviation >
  1e-10
) {
  stop(
    "One or more sensitivity posteriors do not sum to 1."
  )
}


# ------------------------------------------------------------
# 17. Scenario metadata
# ------------------------------------------------------------

scenario_metadata <-
  bind_rows(
    
    lapply(
      names(scenarios),
      function(
    scenario_name
      ) {
        
        config <-
          scenarios[[scenario_name]]
        
        data.frame(
          
          escenario =
            scenario_name,
          
          tipo_escenario =
            config$scenario_type,
          
          evidencia_retirada =
            config$evidence_removed,
          
          pregunta =
            config$question,
          
          stringsAsFactors = FALSE
        )
      }
    )
  )


# ------------------------------------------------------------
# 18. Preliminary MAP-change diagnostic
# ------------------------------------------------------------
#
# This is only an execution diagnostic.
# Formal sensitivity metrics are produced in script 07.
#
# ------------------------------------------------------------

M0_MAP <-
  sensitivity_posteriors |>
  filter(
    escenario == "M0"
  ) |>
  select(
    id,
    UREH_M0 =
      UREH_predicha
  )


MAP_change_preview <-
  bind_rows(
    
    lapply(
      setdiff(
        names(scenarios),
        "M0"
      ),
      function(
    scenario_name
      ) {
        
        scenario_data <-
          sensitivity_posteriors |>
          filter(
            escenario ==
              scenario_name
          ) |>
          select(
            id,
            UREH_scenario =
              UREH_predicha
          ) |>
          left_join(
            M0_MAP,
            by = "id"
          )
        
        data.frame(
          
          escenario =
            scenario_name,
          
          n_cambios_MAP =
            sum(
              scenario_data$UREH_scenario !=
                scenario_data$UREH_M0
            ),
          
          porcentaje_MAP_estable =
            100 *
            mean(
              scenario_data$UREH_scenario ==
                scenario_data$UREH_M0
            ),
          
          stringsAsFactors = FALSE
        )
      }
    )
  )


# ------------------------------------------------------------
# 19. Output paths
# ------------------------------------------------------------

engine_control_file <- here::here(
  "results",
  "tables",
  "sensitivity_engine_control.csv"
)

scenario_metadata_file <- here::here(
  "results",
  "tables",
  "sensitivity_scenarios.csv"
)

sensitivity_posteriors_file <- here::here(
  "results",
  "tables",
  "sensitivity_posteriors.csv"
)


# ------------------------------------------------------------
# 20. Write outputs
# ------------------------------------------------------------

write.csv(
  engine_control,
  engine_control_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  scenario_metadata,
  scenario_metadata_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_posteriors,
  sensitivity_posteriors_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 21. Read-back validation
# ------------------------------------------------------------

engine_control_check <- read.csv(
  engine_control_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

scenario_metadata_check <- read.csv(
  scenario_metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sensitivity_check <- read.csv(
  sensitivity_posteriors_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(
    engine_control_check
  ) != 1L
) {
  stop(
    "Sensitivity engine-control file failed read-back validation."
  )
}

if (
  nrow(
    scenario_metadata_check
  ) != n_scenarios
) {
  stop(
    "Sensitivity scenario metadata failed read-back validation."
  )
}

if (
  nrow(
    sensitivity_check
  ) != expected_rows
) {
  stop(
    "Sensitivity posterior file failed read-back validation."
  )
}

if (
  anyDuplicated(
    sensitivity_check[
      ,
      c(
        "id",
        "escenario"
      )
    ]
  )
) {
  stop(
    "Duplicated id-scenario pairs after read-back."
  )
}

posterior_check_matrix <-
  as.matrix(
    sensitivity_check[
      ,
      cols_ureh,
      drop = FALSE
    ]
  )

storage.mode(
  posterior_check_matrix
) <- "double"

if (
  max(
    abs(
      rowSums(
        posterior_check_matrix
      ) - 1
    )
  ) >
  1e-10
) {
  stop(
    "Sensitivity posterior normalization failed after read-back."
  )
}


# ------------------------------------------------------------
# 22. Console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "SENSITIVITY ANALYSIS COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "Scenarios:",
  paste(
    names(scenarios),
    collapse = " | "
  ),
  "\n"
)

cat(
  "Units per scenario:",
  n_expected_units,
  "\n"
)

cat(
  "Total posterior rows:",
  nrow(
    sensitivity_posteriors
  ),
  "\n"
)

cat(
  "Maximum posterior normalization deviation:",
  format(
    maximum_posterior_deviation,
    scientific = TRUE
  ),
  "\n\n"
)

cat(
  "Preliminary MAP changes relative to M0:\n"
)

print(
  MAP_change_preview,
  row.names = FALSE
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 23. Completion
# ------------------------------------------------------------

message(
  "06_sensitivity_analysis.R completed successfully."
)

message(
  "Generated: results/tables/sensitivity_engine_control.csv"
)

message(
  "Generated: results/tables/sensitivity_scenarios.csv"
)

message(
  "Generated: results/tables/sensitivity_posteriors.csv"
)

message(
  "Next: R/07_summarize_sensitivity.R"
)