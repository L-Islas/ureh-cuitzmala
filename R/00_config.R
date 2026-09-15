# ============================================================
# UREH-CUITZMALA
# Central configuration
# ============================================================
#
# This file defines:
#   - repository paths;
#   - model constants;
#   - Bayesian-network states;
#   - general probability checks;
#   - normalized Shannon entropy.
#
# It contains no statistical analysis.
# ============================================================


# ------------------------------------------------------------
# 1. Repository root
# ------------------------------------------------------------

if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install it with install.packages('here')."
  )
}

repo_root <- here::here()


# ------------------------------------------------------------
# 2. Directory structure
# ------------------------------------------------------------

dir_data_input <- file.path(
  repo_root,
  "data",
  "input"
)

dir_model <- file.path(
  repo_root,
  "data",
  "model"
)

dir_results <- file.path(
  repo_root,
  "results"
)

dir_tables <- file.path(
  dir_results,
  "tables"
)

dir_figures <- file.path(
  dir_results,
  "figures"
)

dir_spatial <- file.path(
  dir_results,
  "spatial"
)

dir_temp <- file.path(
  repo_root,
  "temp"
)

dir_checks <- file.path(
  repo_root,
  "checks"
)

for (d in c(
  dir_results,
  dir_tables,
  dir_figures,
  dir_spatial,
  dir_temp,
  dir_checks
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 3. Main input files
# ------------------------------------------------------------

file_ureh_units <- file.path(
  dir_data_input,
  "ureh_units.gpkg"
)

file_data_dictionary <- file.path(
  dir_data_input,
  "data_dictionary.csv"
)

file_model_states <- file.path(
  dir_model,
  "model_states.csv"
)

cpt_files <- c(
  infiltration = file.path(
    dir_model,
    "cpt_infiltration.csv"
  ),
  runoff = file.path(
    dir_model,
    "cpt_runoff.csv"
  ),
  hydrological_function = file.path(
    dir_model,
    "cpt_hydrological_function.csv"
  ),
  ureh = file.path(
    dir_model,
    "cpt_ureh.csv"
  )
)


# ------------------------------------------------------------
# 4. General constants
# ------------------------------------------------------------

n_expected_units <- 899L
expected_epsg <- 32613L

probability_tolerance <- 1e-10


# ------------------------------------------------------------
# 5. Hydrological discretization
# ------------------------------------------------------------

n_precipitacion <- 3L
n_escorrentia   <- 3L
n_aporte        <- 4L


# ------------------------------------------------------------
# 6. Soft runoff evidence
# ------------------------------------------------------------

# Thresholds used to construct the soft runoff evidence.
# Units correspond to the runoff variable used in the analysis.

b1_esc <- 272431
b2_esc <- 455102
m1_esc <- 675599
m2_esc <- 1021751

# Regularization parameter for soft evidence.

epsilon_esc <- 0.05


# ------------------------------------------------------------
# 7. Bayesian-network states
# ------------------------------------------------------------

estados <- list(
  
  suelo = c(
    "suelo_infiltracion_alta",
    "suelo_infiltracion_media",
    "suelo_infiltracion_baja",
    "suelo_salino"
  ),
  
  cobertura = c(
    "forestal_conservada",
    "forestal_secundaria",
    "hidrofila_humedal",
    "agropecuaria",
    "antropica"
  ),
  
  pendiente = c(
    "pend_baja",
    "pend_media",
    "pend_alta"
  ),
  
  precipitacion = c(
    "prec_baja",
    "prec_media",
    "prec_alta"
  ),
  
  aporte = c(
    "aporte_muy_bajo",
    "aporte_bajo",
    "aporte_medio",
    "aporte_alto"
  ),
  
  infiltracion = c(
    "inf_baja",
    "inf_media",
    "inf_alta"
  ),
  
  escorrentia = c(
    "esc_baja",
    "esc_media",
    "esc_alta"
  ),
  
  funcion_hidrologica = c(
    "aporte",
    "transferencia",
    "recepcion"
  ),
  
  ureh = c(
    "aporte_conservada",
    "aporte_degradada",
    "transferencia_conservada",
    "transferencia_degradada",
    "receptora_conservada",
    "receptora_degradada"
  )
)

cols_ureh <- estados$ureh


# ------------------------------------------------------------
# 8. Probability validation
# ------------------------------------------------------------

validar_suma_uno <- function(
    x,
    tol = probability_tolerance,
    nombre = "objeto"
) {
  
  x <- as.numeric(x)
  
  if (anyNA(x) || any(!is.finite(x))) {
    stop(
      nombre,
      ": contiene valores NA o no finitos."
    )
  }
  
  if (!all(abs(x - 1) < tol)) {
    stop(
      nombre,
      ": hay probabilidades que no suman 1 dentro de la tolerancia."
    )
  }
  
  invisible(TRUE)
}


# ------------------------------------------------------------
# 9. Normalized Shannon entropy
# ------------------------------------------------------------

calcular_entropia_normalizada <- function(
    probabilidades,
    n_clases = length(probabilidades),
    tol = probability_tolerance
) {
  
  p <- as.numeric(probabilidades)
  
  if (
    anyNA(p) ||
    any(!is.finite(p)) ||
    any(p < 0)
  ) {
    stop(
      "Vector de probabilidades inválido."
    )
  }
  
  if (n_clases < 2) {
    stop(
      "n_clases debe ser >= 2."
    )
  }
  
  suma_p <- sum(p)
  
  if (abs(suma_p - 1) > tol) {
    stop(
      "Las probabilidades no suman 1 dentro de la tolerancia."
    )
  }
  
  # Remove zero terms only after fixing the theoretical
  # number of possible states used for normalization.
  p <- p / suma_p
  p_pos <- p[p > 0]
  
  -sum(
    p_pos * log(p_pos)
  ) / log(n_clases)
}


# ------------------------------------------------------------
# 10. Configuration message
# ------------------------------------------------------------

message(
  "UREH configuration loaded."
)