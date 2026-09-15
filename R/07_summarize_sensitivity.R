# ============================================================
# 07. SUMMARIZE SENSITIVITY
# ============================================================
#
# Purpose
# -------
# Compare evidence-removal scenarios against M0 using the three
# metrics defined for the robustness analysis:
#
# 1. MAP stability;
# 2. Jensen-Shannon distance between complete UREH posteriors;
# 3. change in normalized Shannon entropy.
#
# Main scenarios
# --------------
# M1
# M2
# M3
#
# Structural diagnostics
# ----------------------
# Sin cobertura
# Sin suelo + pendiente
#
# IMPORTANT
# ---------
# Diagnostic scenarios are not additional main scenarios and
# must not be interpreted as a ranking of variable importance.
#
# Input
# -----
# results/tables/sensitivity_posteriors.csv
#
# Outputs
# -------
# results/tables/sensitivity_by_unit.csv
# results/tables/sensitivity_summary_principal.csv
# results/tables/sensitivity_summary_diagnostic.csv
# results/tables/sensitivity_summary_all.csv
# results/tables/sensitivity_diagnostic_transitions.csv
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
# 1. Input
# ------------------------------------------------------------

posterior_file <- here::here(
  "results",
  "tables",
  "sensitivity_posteriors.csv"
)

if (!file.exists(posterior_file)) {
  stop(
    "Sensitivity posterior file not found:\n",
    posterior_file,
    "\nRun R/06_sensitivity_analysis.R first."
  )
}

post <- read.csv(
  posterior_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 2. UREH probability columns
# ------------------------------------------------------------

probability_columns <- estados$ureh

if (length(probability_columns) != 6L) {
  stop(
    "Sensitivity summary expects exactly six UREH states."
  )
}


# ------------------------------------------------------------
# 3. Required columns
# ------------------------------------------------------------

required_columns <- c(
  "id",
  "escenario",
  "tipo_escenario",
  "evidencia_retirada",
  "pregunta",
  probability_columns,
  "UREH_predicha",
  "prob_UREH_predicha",
  "entropia_UREH"
)

missing_columns <- setdiff(
  required_columns,
  names(post)
)

if (length(missing_columns) > 0L) {
  stop(
    "Sensitivity posterior file is missing columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 4. Expected scenarios
# ------------------------------------------------------------

expected_scenarios <- c(
  "M0",
  "M1",
  "M2",
  "M3",
  "Sin cobertura",
  "Sin suelo + pendiente"
)

missing_scenarios <- setdiff(
  expected_scenarios,
  unique(post$escenario)
)

if (length(missing_scenarios) > 0L) {
  stop(
    "Missing sensitivity scenarios: ",
    paste(
      missing_scenarios,
      collapse = ", "
    )
  )
}

unexpected_scenarios <- setdiff(
  unique(post$escenario),
  expected_scenarios
)

if (length(unexpected_scenarios) > 0L) {
  stop(
    "Unexpected sensitivity scenarios: ",
    paste(
      unexpected_scenarios,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 5. Validate number of units per scenario
# ------------------------------------------------------------

scenario_counts <- table(
  post$escenario
)

if (!all(
  scenario_counts[expected_scenarios] ==
  n_expected_units
)) {
  stop(
    "One or more scenarios do not contain ",
    n_expected_units,
    " microcatchments:\n",
    paste(
      names(scenario_counts),
      scenario_counts,
      collapse = "; "
    )
  )
}


# ------------------------------------------------------------
# 6. Validate id-scenario combinations
# ------------------------------------------------------------

if (anyNA(post$id)) {
  stop(
    "Sensitivity posterior file contains missing ids."
  )
}

if (
  anyDuplicated(
    post[
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


# ------------------------------------------------------------
# 7. Validate posterior probabilities
# ------------------------------------------------------------

posterior_matrix <- as.matrix(
  post[
    ,
    probability_columns,
    drop = FALSE
  ]
)

storage.mode(
  posterior_matrix
) <- "double"

if (
  anyNA(posterior_matrix) ||
  any(!is.finite(posterior_matrix))
) {
  stop(
    "Sensitivity posterior probabilities contain invalid values."
  )
}

if (
  any(posterior_matrix < -1e-10) ||
  any(posterior_matrix > 1 + 1e-10)
) {
  stop(
    "Sensitivity posterior probabilities fall outside [0,1]."
  )
}

posterior_sum_deviation <- max(
  abs(
    rowSums(
      posterior_matrix
    ) - 1
  )
)

if (
  posterior_sum_deviation >
  1e-9
) {
  stop(
    "UREH posteriors do not sum to 1 within tolerance."
  )
}


# ------------------------------------------------------------
# 8. Validate entropy
# ------------------------------------------------------------

if (
  anyNA(post$entropia_UREH) ||
  any(!is.finite(post$entropia_UREH)) ||
  any(post$entropia_UREH < -1e-10) ||
  any(post$entropia_UREH > 1 + 1e-10)
) {
  stop(
    "Invalid normalized entropy values detected."
  )
}


# ------------------------------------------------------------
# 9. Jensen-Shannon distance
# ------------------------------------------------------------
#
# d_JS(P,Q) = sqrt(JS(P,Q))
#
# JS divergence is calculated using log base 2.
#
# Interpretation:
#
# 0 = identical posterior distributions
# 1 = theoretical maximum distance
#
# ------------------------------------------------------------

jensen_shannon_distance <- function(
    p,
    q
) {
  
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  if (length(p) != length(q)) {
    stop(
      "Jensen-Shannon distributions have different lengths."
    )
  }
  
  if (
    anyNA(p) ||
    anyNA(q) ||
    any(!is.finite(p)) ||
    any(!is.finite(q)) ||
    any(p < 0) ||
    any(q < 0)
  ) {
    stop(
      "Invalid probabilities supplied to Jensen-Shannon distance."
    )
  }
  
  if (
    sum(p) <= 0 ||
    sum(q) <= 0
  ) {
    stop(
      "Zero probability mass supplied to Jensen-Shannon distance."
    )
  }
  
  
  # Normalize defensively
  
  p <- p / sum(p)
  q <- q / sum(q)
  
  
  # Mixture distribution
  
  m <- 0.5 * (p + q)
  
  
  # Kullback-Leibler divergence, base 2
  
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
  
  
  # Numerical protection
  
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
      "Calculated Jensen-Shannon distance outside [0,1]."
    )
  }
  
  distance
}


# ------------------------------------------------------------
# 10. Extract M0 reference
# ------------------------------------------------------------

M0 <- post |>
  filter(
    escenario == "M0"
  ) |>
  arrange(
    id
  )

if (
  nrow(M0) !=
  n_expected_units
) {
  stop(
    "M0 does not contain the expected number of units."
  )
}

if (
  anyDuplicated(
    M0$id
  )
) {
  stop(
    "M0 contains duplicated ids."
  )
}


# ------------------------------------------------------------
# 11. Scenarios compared against M0
# ------------------------------------------------------------

comparison_scenarios <- c(
  "M1",
  "M2",
  "M3",
  "Sin cobertura",
  "Sin suelo + pendiente"
)


# ------------------------------------------------------------
# 12. Compare each scenario with M0 by microcatchment
# ------------------------------------------------------------

comparison_list <- vector(
  mode = "list",
  length =
    length(
      comparison_scenarios
    )
)

names(
  comparison_list
) <- comparison_scenarios


for (
  scenario_name in
  comparison_scenarios
) {
  
  scenario_data <- post |>
    filter(
      escenario ==
        scenario_name
    ) |>
    arrange(
      id
    )
  
  
  # ----------------------------------------------------------
  # 12.1 Exact id correspondence
  # ----------------------------------------------------------
  
  if (!identical(
    as.integer(M0$id),
    as.integer(scenario_data$id)
  )) {
    stop(
      "Ids in scenario '",
      scenario_name,
      "' do not exactly match M0."
    )
  }
  
  
  # ----------------------------------------------------------
  # 12.2 Posterior matrices
  # ----------------------------------------------------------
  
  posterior_M0 <- as.matrix(
    M0[
      ,
      probability_columns,
      drop = FALSE
    ]
  )
  
  posterior_scenario <- as.matrix(
    scenario_data[
      ,
      probability_columns,
      drop = FALSE
    ]
  )
  
  storage.mode(
    posterior_M0
  ) <- "double"
  
  storage.mode(
    posterior_scenario
  ) <- "double"
  
  
  # ----------------------------------------------------------
  # 12.3 Jensen-Shannon distance by unit
  # ----------------------------------------------------------
  
  JS_distance <- vapply(
    seq_len(
      nrow(M0)
    ),
    FUN = function(i) {
      
      jensen_shannon_distance(
        posterior_M0[
          i,
        ],
        posterior_scenario[
          i,
        ]
      )
    },
    FUN.VALUE = numeric(1)
  )
  
  
  # ----------------------------------------------------------
  # 12.4 Assemble unit-level comparison
  # ----------------------------------------------------------
  
  comparison_list[[scenario_name]] <-
    data.frame(
      
      id =
        M0$id,
      
      escenario =
        scenario_name,
      
      tipo_escenario =
        scenario_data$tipo_escenario,
      
      evidencia_retirada =
        scenario_data$evidencia_retirada,
      
      pregunta =
        scenario_data$pregunta,
      
      UREH_M0 =
        M0$UREH_predicha,
      
      UREH_escenario =
        scenario_data$UREH_predicha,
      
      MAP_estable =
        M0$UREH_predicha ==
        scenario_data$UREH_predicha,
      
      distancia_JS =
        JS_distance,
      
      entropia_M0 =
        M0$entropia_UREH,
      
      entropia_escenario =
        scenario_data$entropia_UREH,
      
      delta_entropia =
        scenario_data$entropia_UREH -
        M0$entropia_UREH,
      
      stringsAsFactors = FALSE
    )
}


sensitivity_by_unit <- bind_rows(
  comparison_list
)


# ------------------------------------------------------------
# 13. Validate sensitivity metrics
# ------------------------------------------------------------

metric_columns <- c(
  "MAP_estable",
  "distancia_JS",
  "delta_entropia"
)

if (
  anyNA(
    sensitivity_by_unit[
      ,
      metric_columns,
      drop = FALSE
    ]
  )
) {
  stop(
    "Sensitivity metrics contain missing values."
  )
}

if (
  any(
    !is.finite(
      sensitivity_by_unit$distancia_JS
    )
  )
) {
  stop(
    "Non-finite Jensen-Shannon distances detected."
  )
}

if (
  any(
    sensitivity_by_unit$distancia_JS <
    -1e-12
  ) ||
  any(
    sensitivity_by_unit$distancia_JS >
    1 + 1e-12
  )
) {
  stop(
    "Jensen-Shannon distance outside [0,1]."
  )
}

if (
  any(
    !is.finite(
      sensitivity_by_unit$delta_entropia
    )
  )
) {
  stop(
    "Non-finite entropy changes detected."
  )
}

expected_comparison_rows <-
  n_expected_units *
  length(
    comparison_scenarios
  )

if (
  nrow(
    sensitivity_by_unit
  ) !=
  expected_comparison_rows
) {
  stop(
    "Unexpected number of unit-level sensitivity comparisons."
  )
}

if (
  anyDuplicated(
    sensitivity_by_unit[
      ,
      c(
        "id",
        "escenario"
      )
    ]
  )
) {
  stop(
    "Duplicated unit-scenario comparisons detected."
  )
}


# ------------------------------------------------------------
# 14. Sensitivity-summary function
# ------------------------------------------------------------

summarize_sensitivity <- function(
    data
) {
  
  data |>
    group_by(
      escenario,
      tipo_escenario,
      evidencia_retirada
    ) |>
    summarise(
      
      n_unidades =
        n(),
      
      n_MAP_estable =
        sum(
          MAP_estable
        ),
      
      n_MAP_cambio =
        sum(
          !MAP_estable
        ),
      
      porcentaje_MAP_estable =
        100 *
        mean(
          MAP_estable
        ),
      
      mediana_distancia_JS =
        median(
          distancia_JS
        ),
      
      mediana_delta_entropia =
        median(
          delta_entropia
        ),
      
      .groups = "drop"
    )
}


# ------------------------------------------------------------
# 15. Complete summary
# ------------------------------------------------------------

sensitivity_summary_all <-
  summarize_sensitivity(
    sensitivity_by_unit
  )


# ------------------------------------------------------------
# 16. Main sensitivity analysis: M1-M3
# ------------------------------------------------------------

main_scenario_order <- c(
  "M1",
  "M2",
  "M3"
)

sensitivity_summary_principal <-
  sensitivity_summary_all |>
  filter(
    tipo_escenario ==
      "principal"
  ) |>
  mutate(
    escenario = factor(
      escenario,
      levels =
        main_scenario_order
    )
  ) |>
  arrange(
    escenario
  ) |>
  mutate(
    escenario =
      as.character(
        escenario
      )
  )

if (
  nrow(
    sensitivity_summary_principal
  ) != 3L
) {
  stop(
    "Principal sensitivity summary must contain M1, M2 and M3."
  )
}


# ------------------------------------------------------------
# 17. Structural diagnostic
# ------------------------------------------------------------

diagnostic_scenario_order <- c(
  "Sin cobertura",
  "Sin suelo + pendiente"
)

sensitivity_summary_diagnostic <-
  sensitivity_summary_all |>
  filter(
    tipo_escenario ==
      "diagnostico"
  ) |>
  mutate(
    escenario = factor(
      escenario,
      levels =
        diagnostic_scenario_order
    )
  ) |>
  arrange(
    escenario
  ) |>
  mutate(
    escenario =
      as.character(
        escenario
      )
  )

if (
  nrow(
    sensitivity_summary_diagnostic
  ) != 2L
) {
  stop(
    "Structural diagnostic summary must contain two scenarios."
  )
}


# ------------------------------------------------------------
# 18. Order complete summary
# ------------------------------------------------------------

complete_order <- c(
  "M1",
  "M2",
  "M3",
  "Sin cobertura",
  "Sin suelo + pendiente"
)

sensitivity_summary_all <-
  sensitivity_summary_all |>
  mutate(
    escenario = factor(
      escenario,
      levels =
        complete_order
    )
  ) |>
  arrange(
    escenario
  ) |>
  mutate(
    escenario =
      as.character(
        escenario
      )
  )


# ------------------------------------------------------------
# 19. MAP transitions for structural diagnostic
# ------------------------------------------------------------
#
# This table helps diagnose the mechanism behind M1.
#
# In particular, it allows examination of whether changes caused
# by removal of land-cover evidence occur mainly along the
# conserved/degraded component while hydrological function is
# retained.
#
# ------------------------------------------------------------

sensitivity_diagnostic_transitions <-
  sensitivity_by_unit |>
  filter(
    tipo_escenario ==
      "diagnostico"
  ) |>
  count(
    escenario,
    UREH_M0,
    UREH_escenario,
    name = "n"
  ) |>
  group_by(
    escenario
  ) |>
  mutate(
    
    porcentaje =
      100 *
      n /
      sum(n),
    
    cambio_MAP =
      UREH_M0 !=
      UREH_escenario
  ) |>
  ungroup() |>
  mutate(
    escenario = factor(
      escenario,
      levels =
        diagnostic_scenario_order
    )
  ) |>
  arrange(
    escenario,
    desc(n)
  ) |>
  mutate(
    escenario =
      as.character(
        escenario
      )
  )


# ------------------------------------------------------------
# 20. Output paths
# ------------------------------------------------------------

by_unit_file <- here::here(
  "results",
  "tables",
  "sensitivity_by_unit.csv"
)

principal_summary_file <- here::here(
  "results",
  "tables",
  "sensitivity_summary_principal.csv"
)

diagnostic_summary_file <- here::here(
  "results",
  "tables",
  "sensitivity_summary_diagnostic.csv"
)

complete_summary_file <- here::here(
  "results",
  "tables",
  "sensitivity_summary_all.csv"
)

diagnostic_transitions_file <- here::here(
  "results",
  "tables",
  "sensitivity_diagnostic_transitions.csv"
)


# ------------------------------------------------------------
# 21. Write outputs
# ------------------------------------------------------------

write.csv(
  sensitivity_by_unit,
  by_unit_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_summary_principal,
  principal_summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_summary_diagnostic,
  diagnostic_summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_summary_all,
  complete_summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  sensitivity_diagnostic_transitions,
  diagnostic_transitions_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 22. Read-back validation
# ------------------------------------------------------------

by_unit_check <- read.csv(
  by_unit_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

principal_check <- read.csv(
  principal_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

diagnostic_check <- read.csv(
  diagnostic_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

complete_check <- read.csv(
  complete_summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

transitions_check <- read.csv(
  diagnostic_transitions_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (
  nrow(
    by_unit_check
  ) !=
  expected_comparison_rows
) {
  stop(
    "Unit-level sensitivity output failed read-back validation."
  )
}

if (
  nrow(
    principal_check
  ) != 3L
) {
  stop(
    "Principal sensitivity summary failed read-back validation."
  )
}

if (
  nrow(
    diagnostic_check
  ) != 2L
) {
  stop(
    "Diagnostic sensitivity summary failed read-back validation."
  )
}

if (
  nrow(
    complete_check
  ) != 5L
) {
  stop(
    "Complete sensitivity summary failed read-back validation."
  )
}

if (
  nrow(
    transitions_check
  ) < 2L
) {
  stop(
    "Diagnostic transition table failed read-back validation."
  )
}


# ------------------------------------------------------------
# 23. Console output
# ------------------------------------------------------------
#
# Convert to ordinary data.frame before printing so R displays
# all columns instead of truncating tibble output.
#
# ------------------------------------------------------------

cat(
  "\n====================================================\n"
)

cat(
  "MAIN SENSITIVITY: M1-M3 vs M0\n"
)

cat(
  "====================================================\n\n"
)

print(
  as.data.frame(
    sensitivity_summary_principal
  ),
  row.names = FALSE
)


cat(
  "\n====================================================\n"
)

cat(
  "LOCAL STRUCTURAL-BLOCK DIAGNOSTIC\n"
)

cat(
  "====================================================\n\n"
)

print(
  as.data.frame(
    sensitivity_summary_diagnostic
  ),
  row.names = FALSE
)


cat(
  "\nMetric interpretation:\n"
)

cat(
  "- porcentaje_MAP_estable: percentage of units retaining the M0 MAP class.\n"
)

cat(
  "- distancia_JS = 0: posterior identical to M0.\n"
)

cat(
  "- larger distancia_JS: greater reorganization of the complete posterior.\n"
)

cat(
  "- delta_entropia > 0: normalized entropy increases relative to M0.\n"
)

cat(
  "- delta_entropia < 0: normalized entropy decreases relative to M0.\n\n"
)


# ------------------------------------------------------------
# 24. Completion
# ------------------------------------------------------------

message(
  "07_summarize_sensitivity.R completed successfully."
)

message(
  "Generated: results/tables/sensitivity_by_unit.csv"
)

message(
  "Generated: results/tables/sensitivity_summary_principal.csv"
)

message(
  "Generated: results/tables/sensitivity_summary_diagnostic.csv"
)

message(
  "Generated: results/tables/sensitivity_summary_all.csv"
)

message(
  "Generated: results/tables/sensitivity_diagnostic_transitions.csv"
)