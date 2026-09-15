# ============================================================
# 12. DESCRIPTIVE COMPARISON: C2 CLUSTERS vs M0 UREH
# ============================================================
#
# Purpose
# -------
# Describe the functional differentiation present within the
# multivariate-similarity groups obtained with C2.
#
# C2 represents:
#
#   land cover
#   soil
#   slope
#   precipitation
#   runoff
#   upstream contribution
#
# The Bayesian regionalization represents:
#
#   six UREH categories
#
# This script asks:
#
# - How are the six UREH distributed within each C2 cluster?
# - How are the three hydrological functions distributed within
#   each C2 cluster?
# - How are conserved/degraded conditions distributed within
#   each C2 cluster?
# - Do multivariately similar groups contain more than one
#   inferred hydrological function?
#
# IMPORTANT
# ---------
# This analysis is descriptive.
#
# It is NOT:
#
# - an external validation of the Bayesian network;
# - an accuracy assessment;
# - a cluster-purity assessment;
# - an equivalence test between clustering and UREH;
# - an ARI comparison between clusters and UREH.
#
# The two approaches represent different conceptual objects:
#
# clustering -> multivariate similarity
# UREH       -> probabilistically inferred ecohydrological function
#
# Inputs
# ------
# results/tables/clustering_C1_C2_assignments.csv
# results/tables/M0_UREH_clasificacion_dominante.csv
#
# Outputs
# -------
# results/tables/clustering_C2_UREH_by_unit.csv
# results/tables/clustering_C2_UREH_counts.csv
# results/tables/clustering_C2_UREH_percentages.csv
# results/tables/clustering_C2_function_counts.csv
# results/tables/clustering_C2_function_percentages.csv
# results/tables/clustering_C2_condition_counts.csv
# results/tables/clustering_C2_condition_percentages.csv
# results/tables/clustering_C2_UREH_summary.csv
#
# Figures are generated later in:
#
# R/13_make_figures.R
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(dplyr)
library(tidyr)

results_tables_dir <- here::here(
  "results",
  "tables"
)

dir.create(
  results_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

k_final <- 12L


# ------------------------------------------------------------
# 1. Input paths
# ------------------------------------------------------------

clusters_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_assignments.csv"
)

UREH_file <- here::here(
  "results",
  "tables",
  "M0_UREH_clasificacion_dominante.csv"
)

required_files <- c(
  clusters_file,
  UREH_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (
  length(
    missing_files
  ) > 0L
) {
  stop(
    "Required files not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\nRun R/04_infer_M0.R and R/11_compare_C1_C2.R first."
  )
}


# ------------------------------------------------------------
# 2. Read inputs
# ------------------------------------------------------------

clusters <- read.csv(
  clusters_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

UREH <- read.csv(
  UREH_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Required fields
# ------------------------------------------------------------

required_cluster_fields <- c(
  "id",
  "cluster_C2"
)

required_UREH_fields <- c(
  "id",
  "UREH_predicha"
)

missing_cluster_fields <- setdiff(
  required_cluster_fields,
  names(clusters)
)

missing_UREH_fields <- setdiff(
  required_UREH_fields,
  names(UREH)
)

if (
  length(
    missing_cluster_fields
  ) > 0L
) {
  stop(
    "Clustering assignments are missing fields: ",
    paste(
      missing_cluster_fields,
      collapse = ", "
    )
  )
}

if (
  length(
    missing_UREH_fields
  ) > 0L
) {
  stop(
    "M0 classification is missing fields: ",
    paste(
      missing_UREH_fields,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------
# 4. General validation
# ------------------------------------------------------------

if (
  nrow(clusters) != n_expected_units ||
  nrow(UREH) != n_expected_units
) {
  stop(
    "Clusters and M0 UREH must each contain ",
    n_expected_units,
    " units."
  )
}


if (
  anyNA(clusters$id) ||
  anyNA(UREH$id) ||
  anyDuplicated(clusters$id) ||
  anyDuplicated(UREH$id)
) {
  stop(
    "Invalid ids in clustering or M0 classification."
  )
}


if (!identical(
  sort(
    as.integer(
      clusters$id
    )
  ),
  seq_len(
    n_expected_units
  )
)) {
  stop(
    "Unexpected ids in clustering assignments."
  )
}


if (!identical(
  sort(
    as.integer(
      UREH$id
    )
  ),
  seq_len(
    n_expected_units
  )
)) {
  stop(
    "Unexpected ids in M0 classification."
  )
}


# ------------------------------------------------------------
# 5. Validate C2 clusters
# ------------------------------------------------------------

clusters$cluster_C2 <- as.integer(
  clusters$cluster_C2
)

if (
  anyNA(
    clusters$cluster_C2
  )
) {
  stop(
    "cluster_C2 contains missing or non-integer values."
  )
}


expected_clusters <- seq_len(
  k_final
)

observed_clusters <- sort(
  unique(
    clusters$cluster_C2
  )
)

if (!identical(
  observed_clusters,
  expected_clusters
)) {
  stop(
    "C2 must contain clusters 1 through ",
    k_final,
    "."
  )
}


# ------------------------------------------------------------
# 6. Validate six UREH states
# ------------------------------------------------------------

UREH_levels <- estados$ureh

if (
  length(
    UREH_levels
  ) != 6L
) {
  stop(
    "Expected six UREH states."
  )
}


invalid_UREH <- setdiff(
  unique(
    UREH$UREH_predicha
  ),
  UREH_levels
)

if (
  length(
    invalid_UREH
  ) > 0L
) {
  stop(
    "Unrecognized UREH states: ",
    paste(
      invalid_UREH,
      collapse = ", "
    )
  )
}


if (
  anyNA(
    UREH$UREH_predicha
  )
) {
  stop(
    "M0 classification contains missing UREH states."
  )
}


# ------------------------------------------------------------
# 7. Join C2 and UREH
# ------------------------------------------------------------

C2_UREH_by_unit <-
  clusters |>
  select(
    id,
    cluster_C2
  ) |>
  left_join(
    UREH |>
      select(
        id,
        UREH_predicha
      ),
    by = "id"
  ) |>
  arrange(
    id
  )


if (
  nrow(
    C2_UREH_by_unit
  ) !=
  n_expected_units
) {
  stop(
    "C2-UREH join changed the number of units."
  )
}


if (
  anyNA(
    C2_UREH_by_unit$UREH_predicha
  )
) {
  stop(
    "One or more C2 units could not be matched with M0 UREH."
  )
}


# ------------------------------------------------------------
# 8. Derive UREH semantic axes
# ------------------------------------------------------------
#
# The six UREH combine:
#
# Function:
#   aporte
#   transferencia
#   recepcion
#
# Relative biophysical condition:
#   conservada
#   degradada
#
# These are derived from the final UREH MAP label.
#
# They are NOT an independently calculated MAP classification
# of the hydrological-function node.
#
# ------------------------------------------------------------

C2_UREH_by_unit <-
  C2_UREH_by_unit |>
  mutate(
    
    funcion =
      case_when(
        
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
    
    condicion =
      case_when(
        
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
  anyNA(
    C2_UREH_by_unit$funcion
  ) ||
  anyNA(
    C2_UREH_by_unit$condicion
  )
) {
  stop(
    "Could not derive function or condition from all UREH labels."
  )
}


function_levels <- c(
  "aporte",
  "transferencia",
  "recepcion"
)

condition_levels <- c(
  "conservada",
  "degradada"
)


if (!setequal(
  unique(
    C2_UREH_by_unit$funcion
  ),
  function_levels
)) {
  stop(
    "Unexpected derived hydrological-function states."
  )
}


if (!setequal(
  unique(
    C2_UREH_by_unit$condicion
  ),
  condition_levels
)) {
  stop(
    "Unexpected derived condition states."
  )
}


# ------------------------------------------------------------
# 9. C2 x UREH counts
# ------------------------------------------------------------

C2_UREH_counts <-
  C2_UREH_by_unit |>
  count(
    cluster_C2,
    UREH_predicha,
    name = "n"
  ) |>
  complete(
    cluster_C2 =
      expected_clusters,
    
    UREH_predicha =
      UREH_levels,
    
    fill =
      list(
        n = 0L
      )
  ) |>
  arrange(
    cluster_C2,
    match(
      UREH_predicha,
      UREH_levels
    )
  )


if (
  sum(
    C2_UREH_counts$n
  ) !=
  n_expected_units
) {
  stop(
    "C2 x UREH counts do not sum to the total number of units."
  )
}


# ------------------------------------------------------------
# 10. C2 x UREH percentages
# ------------------------------------------------------------
#
# Percentages are calculated WITHIN each C2 cluster.
#
# Each cluster must sum to 100%.
#
# ------------------------------------------------------------

C2_UREH_percentages <-
  C2_UREH_counts |>
  group_by(
    cluster_C2
  ) |>
  mutate(
    
    n_cluster =
      sum(n),
    
    percentage =
      100 *
      n /
      n_cluster
  ) |>
  ungroup()


UREH_percentage_check <-
  C2_UREH_percentages |>
  group_by(
    cluster_C2
  ) |>
  summarise(
    
    total_percentage =
      sum(
        percentage
      ),
    
    .groups =
      "drop"
  )


if (
  any(
    abs(
      UREH_percentage_check$total_percentage -
      100
    ) >
    1e-9
  )
) {
  stop(
    "UREH percentages do not sum to 100% within every C2 cluster."
  )
}


# ------------------------------------------------------------
# 11. C2 x hydrological-function counts
# ------------------------------------------------------------

C2_function_counts <-
  C2_UREH_by_unit |>
  count(
    cluster_C2,
    funcion,
    name = "n"
  ) |>
  complete(
    cluster_C2 =
      expected_clusters,
    
    funcion =
      function_levels,
    
    fill =
      list(
        n = 0L
      )
  ) |>
  arrange(
    cluster_C2,
    match(
      funcion,
      function_levels
    )
  )


if (
  sum(
    C2_function_counts$n
  ) !=
  n_expected_units
) {
  stop(
    "C2 x function counts do not sum to the total number of units."
  )
}


# ------------------------------------------------------------
# 12. C2 x hydrological-function percentages
# ------------------------------------------------------------

C2_function_percentages <-
  C2_function_counts |>
  group_by(
    cluster_C2
  ) |>
  mutate(
    
    n_cluster =
      sum(n),
    
    percentage =
      100 *
      n /
      n_cluster
  ) |>
  ungroup()


function_percentage_check <-
  C2_function_percentages |>
  group_by(
    cluster_C2
  ) |>
  summarise(
    
    total_percentage =
      sum(
        percentage
      ),
    
    .groups =
      "drop"
  )


if (
  any(
    abs(
      function_percentage_check$total_percentage -
      100
    ) >
    1e-9
  )
) {
  stop(
    "Function percentages do not sum to 100% within every C2 cluster."
  )
}


# ------------------------------------------------------------
# 13. C2 x condition counts
# ------------------------------------------------------------

C2_condition_counts <-
  C2_UREH_by_unit |>
  count(
    cluster_C2,
    condicion,
    name = "n"
  ) |>
  complete(
    cluster_C2 =
      expected_clusters,
    
    condicion =
      condition_levels,
    
    fill =
      list(
        n = 0L
      )
  ) |>
  arrange(
    cluster_C2,
    match(
      condicion,
      condition_levels
    )
  )


if (
  sum(
    C2_condition_counts$n
  ) !=
  n_expected_units
) {
  stop(
    "C2 x condition counts do not sum to the total number of units."
  )
}


# ------------------------------------------------------------
# 14. C2 x condition percentages
# ------------------------------------------------------------

C2_condition_percentages <-
  C2_condition_counts |>
  group_by(
    cluster_C2
  ) |>
  mutate(
    
    n_cluster =
      sum(n),
    
    percentage =
      100 *
      n /
      n_cluster
  ) |>
  ungroup()


condition_percentage_check <-
  C2_condition_percentages |>
  group_by(
    cluster_C2
  ) |>
  summarise(
    
    total_percentage =
      sum(
        percentage
      ),
    
    .groups =
      "drop"
  )


if (
  any(
    abs(
      condition_percentage_check$total_percentage -
      100
    ) >
    1e-9
  )
) {
  stop(
    "Condition percentages do not sum to 100% within every C2 cluster."
  )
}


# ------------------------------------------------------------
# 15. Cross-check against basin-wide M0 classification
# ------------------------------------------------------------
#
# The C2 x UREH contingency table must reproduce the original
# basin-wide M0 class counts after summing over clusters.
#
# ------------------------------------------------------------

M0_counts <- table(
  factor(
    UREH$UREH_predicha,
    levels =
      UREH_levels
  )
)


C2_UREH_total_counts <-
  C2_UREH_counts |>
  group_by(
    UREH_predicha
  ) |>
  summarise(
    
    n =
      sum(n),
    
    .groups =
      "drop"
  ) |>
  arrange(
    match(
      UREH_predicha,
      UREH_levels
    )
  )


if (!identical(
  as.integer(
    C2_UREH_total_counts$n
  ),
  as.integer(
    M0_counts
  )
)) {
  stop(
    "C2 x UREH counts do not reproduce basin-wide M0 frequencies."
  )
}


# ------------------------------------------------------------
# 16. Most frequent UREH in each C2 cluster
# ------------------------------------------------------------
#
# This is a descriptive mode only.
#
# It is NOT interpreted as:
#
# - accuracy;
# - purity;
# - agreement score;
# - validation.
#
# ------------------------------------------------------------

dominant_UREH <-
  C2_UREH_percentages |>
  group_by(
    cluster_C2
  ) |>
  slice_max(
    order_by =
      percentage,
    n =
      1,
    with_ties =
      FALSE
  ) |>
  transmute(
    
    cluster_C2,
    
    n_cluster,
    
    most_frequent_UREH =
      UREH_predicha,
    
    percentage_most_frequent_UREH =
      percentage
  ) |>
  ungroup()


# ------------------------------------------------------------
# 17. Most frequent function in each C2 cluster
# ------------------------------------------------------------

dominant_function <-
  C2_function_percentages |>
  group_by(
    cluster_C2
  ) |>
  slice_max(
    order_by =
      percentage,
    n =
      1,
    with_ties =
      FALSE
  ) |>
  transmute(
    
    cluster_C2,
    
    most_frequent_function =
      funcion,
    
    percentage_most_frequent_function =
      percentage
  ) |>
  ungroup()


# ------------------------------------------------------------
# 18. Most frequent condition in each C2 cluster
# ------------------------------------------------------------

dominant_condition <-
  C2_condition_percentages |>
  group_by(
    cluster_C2
  ) |>
  slice_max(
    order_by =
      percentage,
    n =
      1,
    with_ties =
      FALSE
  ) |>
  transmute(
    
    cluster_C2,
    
    most_frequent_condition =
      condicion,
    
    percentage_most_frequent_condition =
      percentage
  ) |>
  ungroup()


# ------------------------------------------------------------
# 19. Functional diversity within C2 clusters
# ------------------------------------------------------------
#
# Count how many functions and UREH categories are actually
# represented in each C2 cluster.
#
# This is descriptive and does not constitute a diversity index.
#
# ------------------------------------------------------------

represented_categories <-
  C2_UREH_by_unit |>
  group_by(
    cluster_C2
  ) |>
  summarise(
    
    n_UREH_present =
      n_distinct(
        UREH_predicha
      ),
    
    n_functions_present =
      n_distinct(
        funcion
      ),
    
    n_conditions_present =
      n_distinct(
        condicion
      ),
    
    .groups =
      "drop"
  )


# ------------------------------------------------------------
# 20. Cluster-level descriptive summary
# ------------------------------------------------------------

C2_UREH_summary <-
  dominant_UREH |>
  left_join(
    dominant_function,
    by =
      "cluster_C2"
  ) |>
  left_join(
    dominant_condition,
    by =
      "cluster_C2"
  ) |>
  left_join(
    represented_categories,
    by =
      "cluster_C2"
  ) |>
  arrange(
    cluster_C2
  )


if (
  nrow(
    C2_UREH_summary
  ) !=
  k_final
) {
  stop(
    "C2-UREH summary must contain exactly 12 clusters."
  )
}


# ------------------------------------------------------------
# 21. Output paths
# ------------------------------------------------------------

by_unit_file <- here::here(
  "results",
  "tables",
  "clustering_C2_UREH_by_unit.csv"
)

UREH_counts_file <- here::here(
  "results",
  "tables",
  "clustering_C2_UREH_counts.csv"
)

UREH_percentages_file <- here::here(
  "results",
  "tables",
  "clustering_C2_UREH_percentages.csv"
)

function_counts_file <- here::here(
  "results",
  "tables",
  "clustering_C2_function_counts.csv"
)

function_percentages_file <- here::here(
  "results",
  "tables",
  "clustering_C2_function_percentages.csv"
)

condition_counts_file <- here::here(
  "results",
  "tables",
  "clustering_C2_condition_counts.csv"
)

condition_percentages_file <- here::here(
  "results",
  "tables",
  "clustering_C2_condition_percentages.csv"
)

summary_file <- here::here(
  "results",
  "tables",
  "clustering_C2_UREH_summary.csv"
)


# ------------------------------------------------------------
# 22. Write outputs
# ------------------------------------------------------------

write.csv(
  C2_UREH_by_unit,
  by_unit_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_UREH_counts,
  UREH_counts_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_UREH_percentages,
  UREH_percentages_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_function_counts,
  function_counts_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_function_percentages,
  function_percentages_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_condition_counts,
  condition_counts_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_condition_percentages,
  condition_percentages_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  C2_UREH_summary,
  summary_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 23. Read-back validation
# ------------------------------------------------------------

by_unit_check <- read.csv(
  by_unit_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

UREH_counts_check <- read.csv(
  UREH_counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

UREH_percentages_check <- read.csv(
  UREH_percentages_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

function_counts_check <- read.csv(
  function_counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

function_percentages_check <- read.csv(
  function_percentages_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

condition_counts_check <- read.csv(
  condition_counts_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

condition_percentages_check <- read.csv(
  condition_percentages_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_check <- read.csv(
  summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (
  nrow(
    by_unit_check
  ) !=
  n_expected_units
) {
  stop(
    "Unit-level C2-UREH table failed read-back validation."
  )
}


if (
  nrow(
    UREH_counts_check
  ) !=
  k_final *
  length(
    UREH_levels
  )
) {
  stop(
    "C2 x UREH counts failed read-back validation."
  )
}


if (
  nrow(
    UREH_percentages_check
  ) !=
  k_final *
  length(
    UREH_levels
  )
) {
  stop(
    "C2 x UREH percentages failed read-back validation."
  )
}


if (
  nrow(
    function_counts_check
  ) !=
  k_final *
  length(
    function_levels
  )
) {
  stop(
    "C2 x function counts failed read-back validation."
  )
}


if (
  nrow(
    function_percentages_check
  ) !=
  k_final *
  length(
    function_levels
  )
) {
  stop(
    "C2 x function percentages failed read-back validation."
  )
}


if (
  nrow(
    condition_counts_check
  ) !=
  k_final *
  length(
    condition_levels
  )
) {
  stop(
    "C2 x condition counts failed read-back validation."
  )
}


if (
  nrow(
    condition_percentages_check
  ) !=
  k_final *
  length(
    condition_levels
  )
) {
  stop(
    "C2 x condition percentages failed read-back validation."
  )
}


if (
  nrow(
    summary_check
  ) !=
  k_final
) {
  stop(
    "C2-UREH cluster summary failed read-back validation."
  )
}


# ------------------------------------------------------------
# 24. Console output
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "C2 vs UREH — DESCRIPTIVE CLUSTER SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  as.data.frame(
    C2_UREH_summary
  ),
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "C2 vs HYDROLOGICAL FUNCTION — PERCENTAGES\n"
)

cat(
  "============================================\n"
)

print(
  as.data.frame(
    C2_function_percentages
  ),
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "C2 vs RELATIVE CONDITION — PERCENTAGES\n"
)

cat(
  "============================================\n"
)

print(
  as.data.frame(
    C2_condition_percentages
  ),
  row.names = FALSE
)


cat(
  "\n============================================\n\n"
)


# ------------------------------------------------------------
# 25. Completion
# ------------------------------------------------------------

message(
  "12_compare_C2_UREH.R completed successfully."
)

message(
  "Generated: results/tables/clustering_C2_UREH_by_unit.csv"
)

message(
  "Generated: results/tables/clustering_C2_UREH_counts.csv"
)

message(
  "Generated: results/tables/clustering_C2_UREH_percentages.csv"
)

message(
  "Generated: results/tables/clustering_C2_function_counts.csv"
)

message(
  "Generated: results/tables/clustering_C2_function_percentages.csv"
)

message(
  "Generated: results/tables/clustering_C2_condition_counts.csv"
)

message(
  "Generated: results/tables/clustering_C2_condition_percentages.csv"
)

message(
  "Generated: results/tables/clustering_C2_UREH_summary.csv"
)

message(
  "Next: R/13_make_figures.R"
)