# ============================================================
# 03. BUILD BAYESIAN NETWORK
# ============================================================
#
# Purpose
# -------
# Construct and parameterize the Bayesian network used for UREH
# inference.
#
# This script:
#
# 1. reads the prepared microcatchment dataset;
# 2. reads the four operational CPT files;
# 3. validates model states and CPT structure;
# 4. constructs the directed acyclic graph (DAG);
# 5. calculates empirical priors for the five root nodes;
# 6. converts the operational CPT tables to bnlearn arrays;
# 7. constructs the parameterized Bayesian network;
# 8. writes reproducible intermediate products to temp/.
#
# Inputs
# ------
# temp/ureh_prepared.csv
#
# data/model/model_states.csv
# data/model/cpt_infiltration.csv
# data/model/cpt_runoff.csv
# data/model/cpt_hydrological_function.csv
# data/model/cpt_ureh.csv
#
# Outputs
# -------
# temp/bayesian_network_model.rds
# temp/root_priors.csv
# temp/dag_arcs.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(bnlearn)
library(dplyr)
library(tidyr)

temp_dir <- here::here("temp")

dir.create(
  temp_dir,
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

states_file <- here::here(
  "data",
  "model",
  "model_states.csv"
)

cpt_infiltration_file <- here::here(
  "data",
  "model",
  "cpt_infiltration.csv"
)

cpt_runoff_file <- here::here(
  "data",
  "model",
  "cpt_runoff.csv"
)

cpt_function_file <- here::here(
  "data",
  "model",
  "cpt_hydrological_function.csv"
)

cpt_ureh_file <- here::here(
  "data",
  "model",
  "cpt_ureh.csv"
)

required_files <- c(
  prepared_file,
  states_file,
  cpt_infiltration_file,
  cpt_runoff_file,
  cpt_function_file,
  cpt_ureh_file
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
# 2. Read prepared analysis data
# ------------------------------------------------------------

ureh <- read.csv(
  prepared_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

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
    "Missing fields in temp/ureh_prepared.csv: ",
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
    " prepared microcatchments, but found ",
    nrow(ureh),
    "."
  )
}

if (anyNA(ureh$id)) {
  stop(
    "Prepared data contain missing ids."
  )
}

if (anyDuplicated(ureh$id)) {
  stop(
    "Prepared data contain duplicated ids."
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
# 3. Validate model_states.csv against R/00_config.R
# ------------------------------------------------------------

model_states <- read.csv(
  states_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_state_columns <- c(
  "node",
  "state",
  "state_order"
)

if (!identical(
  names(model_states),
  required_state_columns
)) {
  stop(
    "Unexpected columns in model_states.csv."
  )
}

if (anyNA(model_states)) {
  stop(
    "model_states.csv contains missing values."
  )
}

if (anyDuplicated(
  model_states[
    ,
    c("node", "state")
  ]
)) {
  stop(
    "Duplicated node-state combinations in model_states.csv."
  )
}

states_from_config <- do.call(
  rbind,
  lapply(
    names(estados),
    function(node_name) {
      
      data.frame(
        node = node_name,
        state = estados[[node_name]],
        state_order =
          seq_along(
            estados[[node_name]]
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(states_from_config) <- NULL

if (!identical(
  model_states,
  states_from_config
)) {
  stop(
    "model_states.csv does not exactly match ",
    "the state definitions in R/00_config.R."
  )
}

if (nrow(model_states) != 34L) {
  stop(
    "Expected 34 total model states, found ",
    nrow(model_states),
    "."
  )
}

message(
  "Model state definitions validated."
)


# ------------------------------------------------------------
# 4. Validate observed root-node states
# ------------------------------------------------------------

validate_observed_states <- function(
    x,
    allowed,
    variable_name
) {
  
  if (anyNA(x)) {
    stop(
      variable_name,
      " contains missing values."
    )
  }
  
  invalid <- setdiff(
    unique(x),
    allowed
  )
  
  if (length(invalid) > 0L) {
    stop(
      variable_name,
      " contains unrecognized states: ",
      paste(
        invalid,
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
  "Observed root-node states validated."
)


# ------------------------------------------------------------
# 5. Read operational CPT files
# ------------------------------------------------------------

cpt_infiltration_df <- read.csv(
  cpt_infiltration_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cpt_runoff_df <- read.csv(
  cpt_runoff_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cpt_function_df <- read.csv(
  cpt_function_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cpt_ureh_df <- read.csv(
  cpt_ureh_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 6. Generic operational-CPT validator
# ------------------------------------------------------------

validate_operational_cpt <- function(
    x,
    expected_columns,
    parent_columns,
    probability_columns,
    expected_rows,
    cpt_name,
    tolerance = 1e-10
) {
  
  if (!identical(
    names(x),
    expected_columns
  )) {
    stop(
      cpt_name,
      ": unexpected column structure.\n",
      "Expected: ",
      paste(
        expected_columns,
        collapse = ", "
      ),
      "\nFound: ",
      paste(
        names(x),
        collapse = ", "
      )
    )
  }
  
  if (nrow(x) != expected_rows) {
    stop(
      cpt_name,
      ": expected ",
      expected_rows,
      " rows, found ",
      nrow(x),
      "."
    )
  }
  
  if (anyNA(x)) {
    stop(
      cpt_name,
      ": missing values detected."
    )
  }
  
  if (anyDuplicated(
    x[
      ,
      parent_columns,
      drop = FALSE
    ]
  )) {
    stop(
      cpt_name,
      ": duplicated parent-state combinations."
    )
  }
  
  probabilities <- as.matrix(
    x[
      ,
      probability_columns,
      drop = FALSE
    ]
  )
  
  storage.mode(probabilities) <- "double"
  
  if (
    any(!is.finite(probabilities))
  ) {
    stop(
      cpt_name,
      ": non-finite probabilities detected."
    )
  }
  
  if (
    any(
      probabilities <
      -tolerance
    ) ||
    any(
      probabilities >
      1 + tolerance
    )
  ) {
    stop(
      cpt_name,
      ": probabilities outside [0,1]."
    )
  }
  
  row_deviation <- max(
    abs(
      rowSums(probabilities) - 1
    )
  )
  
  if (row_deviation > tolerance) {
    stop(
      cpt_name,
      ": probabilities do not sum to 1. ",
      "Maximum deviation = ",
      row_deviation,
      "."
    )
  }
  
  invisible(TRUE)
}


# ------------------------------------------------------------
# 7. Validate the four operational CPTs
# ------------------------------------------------------------

validate_operational_cpt(
  cpt_infiltration_df,
  expected_columns = c(
    "suelo",
    "cobertura",
    "pendiente",
    "precipitacion",
    "inf_baja",
    "inf_media",
    "inf_alta"
  ),
  parent_columns = c(
    "suelo",
    "cobertura",
    "pendiente",
    "precipitacion"
  ),
  probability_columns = c(
    "inf_baja",
    "inf_media",
    "inf_alta"
  ),
  expected_rows = 180L,
  cpt_name = "CPT infiltration"
)

validate_operational_cpt(
  cpt_runoff_df,
  expected_columns = c(
    "infiltracion",
    "precipitacion",
    "esc_baja",
    "esc_media",
    "esc_alta"
  ),
  parent_columns = c(
    "infiltracion",
    "precipitacion"
  ),
  probability_columns = c(
    "esc_baja",
    "esc_media",
    "esc_alta"
  ),
  expected_rows = 9L,
  cpt_name = "CPT runoff"
)

validate_operational_cpt(
  cpt_function_df,
  expected_columns = c(
    "aporte",
    "escorrentia",
    "aporte_prob",
    "transferencia_prob",
    "recepcion_prob"
  ),
  parent_columns = c(
    "aporte",
    "escorrentia"
  ),
  probability_columns = c(
    "aporte_prob",
    "transferencia_prob",
    "recepcion_prob"
  ),
  expected_rows = 12L,
  cpt_name = "CPT hydrological function"
)

validate_operational_cpt(
  cpt_ureh_df,
  expected_columns = c(
    "funcion_hidrologica",
    "cobertura",
    "aporte_conservada",
    "aporte_degradada",
    "transferencia_conservada",
    "transferencia_degradada",
    "receptora_conservada",
    "receptora_degradada"
  ),
  parent_columns = c(
    "funcion_hidrologica",
    "cobertura"
  ),
  probability_columns =
    estados$ureh,
  expected_rows = 15L,
  cpt_name = "CPT UREH"
)

message(
  "Operational CPT probability structure validated."
)


# ------------------------------------------------------------
# 8. Validate CPT parent states
# ------------------------------------------------------------

if (!setequal(
  cpt_infiltration_df$suelo,
  estados$suelo
)) {
  stop(
    "CPT infiltration contains unexpected soil states."
  )
}

if (!setequal(
  cpt_infiltration_df$cobertura,
  estados$cobertura
)) {
  stop(
    "CPT infiltration contains unexpected land-cover states."
  )
}

if (!setequal(
  cpt_infiltration_df$pendiente,
  estados$pendiente
)) {
  stop(
    "CPT infiltration contains unexpected slope states."
  )
}

if (!setequal(
  cpt_infiltration_df$precipitacion,
  estados$precipitacion
)) {
  stop(
    "CPT infiltration contains unexpected precipitation states."
  )
}

if (!setequal(
  cpt_runoff_df$infiltracion,
  estados$infiltracion
)) {
  stop(
    "CPT runoff contains unexpected infiltration states."
  )
}

if (!setequal(
  cpt_runoff_df$precipitacion,
  estados$precipitacion
)) {
  stop(
    "CPT runoff contains unexpected precipitation states."
  )
}

if (!setequal(
  cpt_function_df$aporte,
  estados$aporte
)) {
  stop(
    "CPT hydrological function contains unexpected ",
    "upstream-contribution states."
  )
}

if (!setequal(
  cpt_function_df$escorrentia,
  estados$escorrentia
)) {
  stop(
    "CPT hydrological function contains unexpected runoff states."
  )
}

if (!setequal(
  cpt_ureh_df$funcion_hidrologica,
  estados$funcion_hidrologica
)) {
  stop(
    "CPT UREH contains unexpected hydrological-function states."
  )
}

if (!setequal(
  cpt_ureh_df$cobertura,
  estados$cobertura
)) {
  stop(
    "CPT UREH contains unexpected land-cover states."
  )
}

message(
  "Operational CPT state definitions validated."
)


# ------------------------------------------------------------
# 9. Verify complete parent-state combinations
# ------------------------------------------------------------

expected_inf_combinations <- expand.grid(
  suelo = estados$suelo,
  cobertura = estados$cobertura,
  pendiente = estados$pendiente,
  precipitacion =
    estados$precipitacion,
  stringsAsFactors = FALSE
)

expected_runoff_combinations <- expand.grid(
  infiltracion =
    estados$infiltracion,
  precipitacion =
    estados$precipitacion,
  stringsAsFactors = FALSE
)

expected_function_combinations <- expand.grid(
  aporte =
    estados$aporte,
  escorrentia =
    estados$escorrentia,
  stringsAsFactors = FALSE
)

expected_ureh_combinations <- expand.grid(
  funcion_hidrologica =
    estados$funcion_hidrologica,
  cobertura =
    estados$cobertura,
  stringsAsFactors = FALSE
)

combination_key <- function(
    x,
    columns
) {
  apply(
    x[
      ,
      columns,
      drop = FALSE
    ],
    1,
    paste,
    collapse = "||"
  )
}

if (!setequal(
  combination_key(
    cpt_infiltration_df,
    c(
      "suelo",
      "cobertura",
      "pendiente",
      "precipitacion"
    )
  ),
  combination_key(
    expected_inf_combinations,
    c(
      "suelo",
      "cobertura",
      "pendiente",
      "precipitacion"
    )
  )
)) {
  stop(
    "CPT infiltration does not contain the complete ",
    "parent-state Cartesian product."
  )
}

if (!setequal(
  combination_key(
    cpt_runoff_df,
    c(
      "infiltracion",
      "precipitacion"
    )
  ),
  combination_key(
    expected_runoff_combinations,
    c(
      "infiltracion",
      "precipitacion"
    )
  )
)) {
  stop(
    "CPT runoff does not contain the complete ",
    "parent-state Cartesian product."
  )
}

if (!setequal(
  combination_key(
    cpt_function_df,
    c(
      "aporte",
      "escorrentia"
    )
  ),
  combination_key(
    expected_function_combinations,
    c(
      "aporte",
      "escorrentia"
    )
  )
)) {
  stop(
    "CPT hydrological function does not contain the complete ",
    "parent-state Cartesian product."
  )
}

if (!setequal(
  combination_key(
    cpt_ureh_df,
    c(
      "funcion_hidrologica",
      "cobertura"
    )
  ),
  combination_key(
    expected_ureh_combinations,
    c(
      "funcion_hidrologica",
      "cobertura"
    )
  )
)) {
  stop(
    "CPT UREH does not contain the complete ",
    "parent-state Cartesian product."
  )
}

message(
  "All parent-state combinations are complete."
)


# ------------------------------------------------------------
# 10. Construct the DAG
# ------------------------------------------------------------

dag <- model2network(
  paste0(
    "[suelo]",
    "[cobertura]",
    "[pendiente]",
    "[precipitacion]",
    "[aporte]",
    "[infiltracion|",
    "suelo:cobertura:pendiente:precipitacion]",
    "[escorrentia|",
    "infiltracion:precipitacion]",
    "[funcion_hidrologica|",
    "aporte:escorrentia]",
    "[ureh|",
    "funcion_hidrologica:cobertura]"
  )
)


# ------------------------------------------------------------
# 11. Validate DAG nodes and arcs
# ------------------------------------------------------------

expected_nodes <- names(
  estados
)

if (!setequal(
  nodes(dag),
  expected_nodes
)) {
  stop(
    "The DAG nodes do not match the model-state definitions."
  )
}

expected_arcs <- data.frame(
  from = c(
    "suelo",
    "cobertura",
    "pendiente",
    "precipitacion",
    "infiltracion",
    "precipitacion",
    "aporte",
    "escorrentia",
    "funcion_hidrologica",
    "cobertura"
  ),
  to = c(
    "infiltracion",
    "infiltracion",
    "infiltracion",
    "infiltracion",
    "escorrentia",
    "escorrentia",
    "funcion_hidrologica",
    "funcion_hidrologica",
    "ureh",
    "ureh"
  ),
  stringsAsFactors = FALSE
)

dag_arcs <- as.data.frame(
  arcs(dag),
  stringsAsFactors = FALSE
)

names(dag_arcs) <- c(
  "from",
  "to"
)

arc_key <- function(x) {
  paste(
    x$from,
    x$to,
    sep = "->"
  )
}

if (!setequal(
  arc_key(dag_arcs),
  arc_key(expected_arcs)
)) {
  stop(
    "The constructed DAG does not contain the expected arcs."
  )
}

if (!bnlearn::acyclic(dag)) {
  stop(
    "The Bayesian-network graph is not acyclic."
  )
}

message(
  "Bayesian-network DAG validated."
)


# ------------------------------------------------------------
# 12. Convert infiltration CPT to bnlearn array
# ------------------------------------------------------------

inf_long <- cpt_infiltration_df |>
  pivot_longer(
    cols = all_of(
      estados$infiltracion
    ),
    names_to = "infiltracion",
    values_to = "probabilidad"
  ) |>
  mutate(
    infiltracion = factor(
      infiltracion,
      levels =
        estados$infiltracion
    ),
    cobertura = factor(
      cobertura,
      levels =
        estados$cobertura
    ),
    pendiente = factor(
      pendiente,
      levels =
        estados$pendiente
    ),
    precipitacion = factor(
      precipitacion,
      levels =
        estados$precipitacion
    ),
    suelo = factor(
      suelo,
      levels =
        estados$suelo
    )
  )

if (anyNA(inf_long)) {
  stop(
    "CPT infiltration contains invalid states after factorization."
  )
}

cpt_infiltracion <- xtabs(
  probabilidad ~
    infiltracion +
    cobertura +
    pendiente +
    precipitacion +
    suelo,
  data = inf_long
)


# ------------------------------------------------------------
# 13. Convert runoff CPT to bnlearn array
# ------------------------------------------------------------

runoff_long <- cpt_runoff_df |>
  pivot_longer(
    cols = all_of(
      estados$escorrentia
    ),
    names_to = "escorrentia",
    values_to = "probabilidad"
  ) |>
  mutate(
    escorrentia = factor(
      escorrentia,
      levels =
        estados$escorrentia
    ),
    infiltracion = factor(
      infiltracion,
      levels =
        estados$infiltracion
    ),
    precipitacion = factor(
      precipitacion,
      levels =
        estados$precipitacion
    )
  )

if (anyNA(runoff_long)) {
  stop(
    "CPT runoff contains invalid states after factorization."
  )
}

cpt_escorrentia <- xtabs(
  probabilidad ~
    escorrentia +
    infiltracion +
    precipitacion,
  data = runoff_long
)


# ------------------------------------------------------------
# 14. Convert hydrological-function CPT to bnlearn array
# ------------------------------------------------------------

function_long <- cpt_function_df |>
  pivot_longer(
    cols = c(
      "aporte_prob",
      "transferencia_prob",
      "recepcion_prob"
    ),
    names_to = "funcion_hidrologica",
    values_to = "probabilidad"
  ) |>
  mutate(
    funcion_hidrologica = recode(
      funcion_hidrologica,
      "aporte_prob" =
        "aporte",
      "transferencia_prob" =
        "transferencia",
      "recepcion_prob" =
        "recepcion"
    ),
    funcion_hidrologica = factor(
      funcion_hidrologica,
      levels =
        estados$funcion_hidrologica
    ),
    aporte = factor(
      aporte,
      levels =
        estados$aporte
    ),
    escorrentia = factor(
      escorrentia,
      levels =
        estados$escorrentia
    )
  )

if (anyNA(function_long)) {
  stop(
    "CPT hydrological function contains invalid states ",
    "after factorization."
  )
}

cpt_funcion_hidrologica <- xtabs(
  probabilidad ~
    funcion_hidrologica +
    aporte +
    escorrentia,
  data = function_long
)


# ------------------------------------------------------------
# 15. Convert UREH CPT to bnlearn array
# ------------------------------------------------------------

ureh_long <- cpt_ureh_df |>
  pivot_longer(
    cols = all_of(
      estados$ureh
    ),
    names_to = "ureh",
    values_to = "probabilidad"
  ) |>
  mutate(
    ureh = factor(
      ureh,
      levels =
        estados$ureh
    ),
    funcion_hidrologica =
      factor(
        funcion_hidrologica,
        levels =
          estados$funcion_hidrologica
      ),
    cobertura = factor(
      cobertura,
      levels =
        estados$cobertura
    )
  )

if (anyNA(ureh_long)) {
  stop(
    "CPT UREH contains invalid states after factorization."
  )
}

cpt_ureh <- xtabs(
  probabilidad ~
    ureh +
    funcion_hidrologica +
    cobertura,
  data = ureh_long
)


# ------------------------------------------------------------
# 16. Validate CPT array dimensions
# ------------------------------------------------------------

validate_cpt_dimensions <- function(
    cpt,
    expected_dimensions,
    expected_states,
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
      ": unexpected dimension order.\n",
      "Expected: ",
      paste(
        expected_dimensions,
        collapse = ", "
      ),
      "\nFound: ",
      paste(
        actual_dimensions,
        collapse = ", "
      )
    )
  }
  
  for (
    i in seq_along(
      expected_dimensions
    )
  ) {
    
    actual_states <-
      dimnames(cpt)[[i]]
    
    expected_i <-
      expected_states[[i]]
    
    if (!identical(
      actual_states,
      expected_i
    )) {
      stop(
        cpt_name,
        ": unexpected states in dimension ",
        expected_dimensions[i],
        "."
      )
    }
  }
  
  invisible(TRUE)
}

validate_cpt_dimensions(
  cpt_infiltracion,
  expected_dimensions = c(
    "infiltracion",
    "cobertura",
    "pendiente",
    "precipitacion",
    "suelo"
  ),
  expected_states = list(
    estados$infiltracion,
    estados$cobertura,
    estados$pendiente,
    estados$precipitacion,
    estados$suelo
  ),
  cpt_name = "CPT infiltration"
)

validate_cpt_dimensions(
  cpt_escorrentia,
  expected_dimensions = c(
    "escorrentia",
    "infiltracion",
    "precipitacion"
  ),
  expected_states = list(
    estados$escorrentia,
    estados$infiltracion,
    estados$precipitacion
  ),
  cpt_name = "CPT runoff"
)

validate_cpt_dimensions(
  cpt_funcion_hidrologica,
  expected_dimensions = c(
    "funcion_hidrologica",
    "aporte",
    "escorrentia"
  ),
  expected_states = list(
    estados$funcion_hidrologica,
    estados$aporte,
    estados$escorrentia
  ),
  cpt_name =
    "CPT hydrological function"
)

validate_cpt_dimensions(
  cpt_ureh,
  expected_dimensions = c(
    "ureh",
    "funcion_hidrologica",
    "cobertura"
  ),
  expected_states = list(
    estados$ureh,
    estados$funcion_hidrologica,
    estados$cobertura
  ),
  cpt_name = "CPT UREH"
)


# ------------------------------------------------------------
# 17. Validate CPT array probability sums
# ------------------------------------------------------------

sum_infiltration <- apply(
  cpt_infiltracion,
  c(2, 3, 4, 5),
  sum
)

sum_runoff <- apply(
  cpt_escorrentia,
  c(2, 3),
  sum
)

sum_function <- apply(
  cpt_funcion_hidrologica,
  c(2, 3),
  sum
)

sum_ureh <- apply(
  cpt_ureh,
  c(2, 3),
  sum
)

validar_suma_uno(
  sum_infiltration,
  nombre =
    "CPT infiltration"
)

validar_suma_uno(
  sum_runoff,
  nombre =
    "CPT runoff"
)

validar_suma_uno(
  sum_function,
  nombre =
    "CPT hydrological function"
)

validar_suma_uno(
  sum_ureh,
  nombre =
    "CPT UREH"
)

message(
  "bnlearn CPT arrays validated."
)


# ------------------------------------------------------------
# 18. Empirical priors for root nodes
# ------------------------------------------------------------
#
# These probabilities are calculated from the prepared
# microcatchment dataset and are therefore not fixed CPT files.
# ------------------------------------------------------------

prior_suelo <- prop.table(
  table(
    factor(
      ureh$suelo,
      levels =
        estados$suelo
    )
  )
)

prior_cobertura <- prop.table(
  table(
    factor(
      ureh$cobertura,
      levels =
        estados$cobertura
    )
  )
)

prior_pendiente <- prop.table(
  table(
    factor(
      ureh$pendiente,
      levels =
        estados$pendiente
    )
  )
)

prior_precipitacion <- prop.table(
  table(
    factor(
      ureh$precipitacion,
      levels =
        estados$precipitacion
    )
  )
)

prior_aporte <- prop.table(
  table(
    factor(
      ureh$aporte,
      levels =
        estados$aporte
    )
  )
)


# ------------------------------------------------------------
# 19. Validate root priors
# ------------------------------------------------------------

priors <- list(
  suelo =
    prior_suelo,
  cobertura =
    prior_cobertura,
  pendiente =
    prior_pendiente,
  precipitacion =
    prior_precipitacion,
  aporte =
    prior_aporte
)

for (
  node_name in names(priors)
) {
  
  p <- priors[[node_name]]
  
  if (anyNA(p)) {
    stop(
      "Prior ",
      node_name,
      " contains missing values."
    )
  }
  
  if (
    any(!is.finite(p))
  ) {
    stop(
      "Prior ",
      node_name,
      " contains non-finite values."
    )
  }
  
  if (
    any(p < 0) ||
    any(p > 1)
  ) {
    stop(
      "Prior ",
      node_name,
      " contains probabilities outside [0,1]."
    )
  }
  
  if (
    abs(
      sum(p) - 1
    ) >
    1e-10
  ) {
    stop(
      "Prior ",
      node_name,
      " does not sum to 1."
    )
  }
  
  if (!identical(
    names(p),
    estados[[node_name]]
  )) {
    stop(
      "Prior ",
      node_name,
      " has unexpected state names/order."
    )
  }
}

message(
  "Empirical root-node priors validated."
)


# ------------------------------------------------------------
# 20. Create tidy root-prior table
# ------------------------------------------------------------

root_priors <- bind_rows(
  lapply(
    names(priors),
    function(node_name) {
      
      data.frame(
        node =
          node_name,
        state =
          names(
            priors[[node_name]]
          ),
        probability =
          as.numeric(
            priors[[node_name]]
          ),
        stringsAsFactors = FALSE
      )
    }
  )
)

if (
  abs(
    sum(
      root_priors$probability[
        root_priors$node ==
        "suelo"
      ]
    ) - 1
  ) >
  1e-10
) {
  stop(
    "Root-prior output validation failed."
  )
}


# ------------------------------------------------------------
# 21. Parameterize Bayesian network
# ------------------------------------------------------------

fit_ureh <- custom.fit(
  dag,
  dist = list(
    
    suelo =
      prior_suelo,
    
    cobertura =
      prior_cobertura,
    
    pendiente =
      prior_pendiente,
    
    precipitacion =
      prior_precipitacion,
    
    aporte =
      prior_aporte,
    
    infiltracion =
      cpt_infiltracion,
    
    escorrentia =
      cpt_escorrentia,
    
    funcion_hidrologica =
      cpt_funcion_hidrologica,
    
    ureh =
      cpt_ureh
  )
)


# ------------------------------------------------------------
# 22. Validate fitted network
# ------------------------------------------------------------

if (!inherits(
  fit_ureh,
  "bn.fit"
)) {
  stop(
    "custom.fit() did not return a bn.fit object."
  )
}

if (!setequal(
  names(fit_ureh),
  expected_nodes
)) {
  stop(
    "The fitted network does not contain the expected nodes."
  )
}

for (
  node_name in expected_nodes
) {
  
  if (is.null(
    fit_ureh[[node_name]]$prob
  )) {
    stop(
      "Node ",
      node_name,
      " has no probability distribution in the fitted network."
    )
  }
}

message(
  "Parameterized Bayesian network validated."
)


# ------------------------------------------------------------
# 23. Assemble model object
# ------------------------------------------------------------

bayesian_network_model <- list(
  
  dag =
    dag,
  
  fit_ureh =
    fit_ureh,
  
  estados =
    estados,
  
  cpt_infiltracion =
    cpt_infiltracion,
  
  cpt_escorrentia =
    cpt_escorrentia,
  
  cpt_funcion_hidrologica =
    cpt_funcion_hidrologica,
  
  cpt_ureh =
    cpt_ureh,
  
  priors =
    priors
)


# ------------------------------------------------------------
# 24. Output paths
# ------------------------------------------------------------

model_output_file <- here::here(
  "temp",
  "bayesian_network_model.rds"
)

priors_output_file <- here::here(
  "temp",
  "root_priors.csv"
)

arcs_output_file <- here::here(
  "temp",
  "dag_arcs.csv"
)


# ------------------------------------------------------------
# 25. Write outputs
# ------------------------------------------------------------

saveRDS(
  bayesian_network_model,
  model_output_file
)

write.csv(
  root_priors,
  priors_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  dag_arcs,
  arcs_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 26. Read-back validation
# ------------------------------------------------------------

model_check <- readRDS(
  model_output_file
)

priors_check <- read.csv(
  priors_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

arcs_check <- read.csv(
  arcs_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!inherits(
  model_check$fit_ureh,
  "bn.fit"
)) {
  stop(
    "Model RDS failed read-back validation."
  )
}

if (!identical(
  model_check$estados,
  estados
)) {
  stop(
    "State definitions changed after model RDS read-back."
  )
}

if (nrow(priors_check) !=
    sum(
      lengths(
        estados[
          c(
            "suelo",
            "cobertura",
            "pendiente",
            "precipitacion",
            "aporte"
          )
        ]
      )
    )) {
  stop(
    "Unexpected number of rows in root_priors.csv."
  )
}

if (!setequal(
  arc_key(arcs_check),
  arc_key(expected_arcs)
)) {
  stop(
    "DAG arcs failed CSV read-back validation."
  )
}


# ------------------------------------------------------------
# 27. Diagnostic summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "BAYESIAN NETWORK SUMMARY\n"
)

cat(
  "============================================\n"
)

cat(
  "Nodes:",
  length(
    nodes(dag)
  ),
  "\n"
)

cat(
  "Arcs:",
  nrow(
    dag_arcs
  ),
  "\n"
)

cat(
  "Root nodes:",
  paste(
    names(priors),
    collapse = ", "
  ),
  "\n"
)

cat(
  "\nEmpirical root priors:\n"
)

print(
  root_priors
)

cat(
  "\nCPT dimensions:\n"
)

cat(
  "Infiltration:",
  paste(
    dim(
      cpt_infiltracion
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Runoff:",
  paste(
    dim(
      cpt_escorrentia
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Hydrological function:",
  paste(
    dim(
      cpt_funcion_hidrologica
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "UREH:",
  paste(
    dim(
      cpt_ureh
    ),
    collapse = " x "
  ),
  "\n"
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 28. Completion
# ------------------------------------------------------------

message(
  "03_build_bayesian_network.R completed successfully."
)

message(
  "Generated: temp/bayesian_network_model.rds"
)

message(
  "Generated: temp/root_priors.csv"
)

message(
  "Generated: temp/dag_arcs.csv"
)