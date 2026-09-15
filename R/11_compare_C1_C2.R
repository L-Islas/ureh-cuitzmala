# ============================================================
# 11. COMPARE C1 AND C2 CLUSTERING SOLUTIONS
# ============================================================
#
# Purpose
# -------
# Compare the finalized k = 12 PAM solutions:
#
# C1:
#   local bio-physical similarity
#
# C2:
#   local bio-physical similarity + upstream contribution
#
# The script:
#
# 1. extracts k = 12 assignments for C1 and C2;
# 2. calculates the Adjusted Rand Index (ARI);
# 3. aligns C2 labels with C1 labels for interpretation;
# 4. identifies units reassigned between C1 and C2;
# 5. summarizes real C1 -> C2 transitions;
# 6. generates compact cluster profiles.
#
# IMPORTANT
# ---------
# Cluster numbers have no intrinsic meaning.
#
# Label alignment is performed only to make the two partitions
# directly readable. ARI itself is invariant to cluster labels.
#
# This comparison is an independent structural reference.
# It is NOT interpreted as external validation of the Bayesian
# network or of UREH.
#
# Inputs
# ------
# temp/clustering_base.csv
# results/tables/clustering_assignments_k2_20.csv
# results/tables/clustering_k_evaluation.csv
#
# Outputs
# -------
# results/tables/clustering_C1_C2_assignments.csv
# results/tables/clustering_C1_C2_summary.csv
# results/tables/clustering_C2_to_C1_labels.csv
# results/tables/clustering_C1_C2_transition_matrix.csv
# results/tables/clustering_C1_C2_transitions.csv
# results/tables/clustering_C1_C2_reassigned_units.csv
# results/tables/clustering_C1_profiles.csv
# results/tables/clustering_C2_profiles.csv
# results/tables/clustering_profiles_comparable.csv
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

k_final <- 12L


# ------------------------------------------------------------
# 1. Input paths
# ------------------------------------------------------------

base_file <- here::here(
  "temp",
  "clustering_base.csv"
)

assignments_file <- here::here(
  "results",
  "tables",
  "clustering_assignments_k2_20.csv"
)

evaluation_file <- here::here(
  "results",
  "tables",
  "clustering_k_evaluation.csv"
)

required_files <- c(
  base_file,
  assignments_file,
  evaluation_file
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
    "Required clustering files not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\nRun R/09_prepare_clustering.R and R/10_gower_PAM.R first."
  )
}


# ------------------------------------------------------------
# 2. Read inputs
# ------------------------------------------------------------

clustering_base <- read.csv(
  base_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

assignments <- read.csv(
  assignments_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

k_evaluation <- read.csv(
  evaluation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 3. Validate clustering base
# ------------------------------------------------------------

required_base_fields <- c(
  "id",
  "USyV",
  "Suelos",
  "Grad_Pend",
  "P_sum",
  "Esc_sum",
  "Agua_D",
  "Agua_D_log"
)

missing_base_fields <- setdiff(
  required_base_fields,
  names(
    clustering_base
  )
)

if (
  length(
    missing_base_fields
  ) > 0L
) {
  stop(
    "Clustering base is missing fields: ",
    paste(
      missing_base_fields,
      collapse = ", "
    )
  )
}


if (
  nrow(
    clustering_base
  ) !=
  n_expected_units
) {
  stop(
    "Clustering base must contain ",
    n_expected_units,
    " units."
  )
}


if (
  anyNA(
    clustering_base$id
  ) ||
  anyDuplicated(
    clustering_base$id
  )
) {
  stop(
    "Invalid ids in clustering base."
  )
}


# ------------------------------------------------------------
# 4. Validate assignment table
# ------------------------------------------------------------

required_assignment_fields <- c(
  "id",
  "escenario",
  "k",
  "cluster",
  "silhouette"
)

missing_assignment_fields <- setdiff(
  required_assignment_fields,
  names(
    assignments
  )
)

if (
  length(
    missing_assignment_fields
  ) > 0L
) {
  stop(
    "Clustering assignments are missing fields: ",
    paste(
      missing_assignment_fields,
      collapse = ", "
    )
  )
}


if (
  anyNA(
    assignments[
      ,
      required_assignment_fields,
      drop = FALSE
    ]
  )
) {
  stop(
    "Clustering assignments contain missing required values."
  )
}


if (
  anyDuplicated(
    assignments[
      ,
      c(
        "id",
        "escenario",
        "k"
      )
    ]
  )
) {
  stop(
    "Duplicated id-scenario-k assignments detected."
  )
}


# ------------------------------------------------------------
# 5. Validate k = 12 evaluation
# ------------------------------------------------------------

k12_evaluation <-
  k_evaluation |>
  filter(
    k == k_final,
    escenario %in%
      c(
        "C1",
        "C2"
      )
  ) |>
  arrange(
    escenario
  )


if (
  nrow(
    k12_evaluation
  ) != 2L
) {
  stop(
    "Expected exactly one k = 12 evaluation for C1 and C2."
  )
}


# ------------------------------------------------------------
# 6. Extract finalized k = 12 solutions
# ------------------------------------------------------------

solution_C1 <-
  assignments |>
  filter(
    escenario == "C1",
    k == k_final
  ) |>
  transmute(
    
    id,
    
    cluster_C1 =
      as.integer(
        cluster
      ),
    
    silhouette_C1 =
      silhouette
  ) |>
  arrange(
    id
  )


solution_C2 <-
  assignments |>
  filter(
    escenario == "C2",
    k == k_final
  ) |>
  transmute(
    
    id,
    
    cluster_C2_original =
      as.integer(
        cluster
      ),
    
    silhouette_C2 =
      silhouette
  ) |>
  arrange(
    id
  )


if (
  nrow(
    solution_C1
  ) !=
  n_expected_units ||
  nrow(
    solution_C2
  ) !=
  n_expected_units
) {
  stop(
    "C1 and C2 must each contain ",
    n_expected_units,
    " assignments at k = 12."
  )
}


if (!identical(
  as.integer(
    solution_C1$id
  ),
  as.integer(
    solution_C2$id
  )
)) {
  stop(
    "C1 and C2 ids do not match at k = 12."
  )
}


if (
  length(
    unique(
      solution_C1$cluster_C1
    )
  ) !=
  k_final ||
  length(
    unique(
      solution_C2$cluster_C2_original
    )
  ) !=
  k_final
) {
  stop(
    "C1 and C2 must each contain exactly 12 clusters."
  )
}


# ------------------------------------------------------------
# 7. Adjusted Rand Index
# ------------------------------------------------------------
#
# Direct implementation avoids an additional package dependency.
#
# ARI compares two partitions and is invariant to the arbitrary
# numerical labels assigned to clusters.
#
# ------------------------------------------------------------

comb2 <- function(n) {
  
  n <- as.numeric(
    n
  )
  
  ifelse(
    n < 2,
    0,
    n *
      (n - 1) /
      2
  )
}


adjusted_rand_index <- function(
    x,
    y
) {
  
  if (
    length(x) !=
    length(y)
  ) {
    stop(
      "ARI partitions have different lengths."
    )
  }
  
  if (
    anyNA(x) ||
    anyNA(y)
  ) {
    stop(
      "ARI partitions contain missing values."
    )
  }
  
  
  contingency <- table(
    x,
    y
  )
  
  n <- sum(
    contingency
  )
  
  sum_pairs_cells <- sum(
    comb2(
      contingency
    )
  )
  
  sum_pairs_rows <- sum(
    comb2(
      rowSums(
        contingency
      )
    )
  )
  
  sum_pairs_columns <- sum(
    comb2(
      colSums(
        contingency
      )
    )
  )
  
  total_pairs <- comb2(
    n
  )
  
  expected_index <-
    (
      sum_pairs_rows *
        sum_pairs_columns
    ) /
    total_pairs
  
  maximum_index <-
    0.5 *
    (
      sum_pairs_rows +
        sum_pairs_columns
    )
  
  
  if (
    maximum_index ==
    expected_index
  ) {
    return(
      1
    )
  }
  
  
  (
    sum_pairs_cells -
      expected_index
  ) /
    (
      maximum_index -
        expected_index
    )
}


ARI_C1_C2 <- adjusted_rand_index(
  solution_C1$cluster_C1,
  solution_C2$cluster_C2_original
)


if (
  !is.finite(
    ARI_C1_C2
  ) ||
  ARI_C1_C2 < -1 ||
  ARI_C1_C2 > 1
) {
  stop(
    "Invalid Adjusted Rand Index."
  )
}


# ------------------------------------------------------------
# 8. Raw C1 x C2 correspondence
# ------------------------------------------------------------

raw_correspondence <- table(
  
  C1 =
    solution_C1$cluster_C1,
  
  C2 =
    solution_C2$cluster_C2_original
)


# ------------------------------------------------------------
# 9. Align C2 labels with C1
# ------------------------------------------------------------
#
# Cluster labels are arbitrary.
#
# For each C2 cluster, identify the C1 cluster containing the
# largest number of the same units.
#
# The resulting correspondence must be:
#
# - unique for every C2 cluster;
# - one-to-one;
# - complete across all 12 clusters.
#
# If these conditions fail, the script stops rather than silently
# forcing an ambiguous relabeling.
#
# ------------------------------------------------------------

C2_clusters <- sort(
  unique(
    solution_C2$cluster_C2_original
  )
)


label_correspondence_list <- lapply(
  C2_clusters,
  function(
    C2_cluster
  ) {
    
    column_counts <-
      raw_correspondence[
        ,
        as.character(
          C2_cluster
        )
      ]
    
    
    maximum_overlap <- max(
      column_counts
    )
    
    
    candidates <- as.integer(
      names(
        column_counts
      )[
        column_counts ==
          maximum_overlap
      ]
    )
    
    
    if (
      length(
        candidates
      ) != 1L
    ) {
      stop(
        "Ambiguous C2 -> C1 label correspondence for C2 cluster ",
        C2_cluster,
        "."
      )
    }
    
    
    data.frame(
      
      cluster_C2_original =
        C2_cluster,
      
      cluster_C1_equivalent =
        candidates,
      
      n_overlap =
        maximum_overlap,
      
      stringsAsFactors = FALSE
    )
  }
)


label_correspondence <-
  bind_rows(
    label_correspondence_list
  )


# One-to-one correspondence

if (
  anyDuplicated(
    label_correspondence$
    cluster_C1_equivalent
  )
) {
  stop(
    paste0(
      "Simple C2 -> C1 label alignment is not one-to-one. ",
      "An explicit optimal assignment would be required."
    )
  )
}


if (
  nrow(
    label_correspondence
  ) !=
  k_final ||
  length(
    unique(
      label_correspondence$
      cluster_C1_equivalent
    )
  ) !=
  k_final
) {
  stop(
    "Incomplete C2 -> C1 label correspondence."
  )
}


# ------------------------------------------------------------
# 10. Construct final aligned assignment table
# ------------------------------------------------------------

clustering_C1_C2_assignments <-
  solution_C1 |>
  left_join(
    solution_C2,
    by = "id"
  ) |>
  left_join(
    label_correspondence,
    by =
      "cluster_C2_original"
  ) |>
  rename(
    cluster_C2 =
      cluster_C1_equivalent
  ) |>
  mutate(
    
    changed_cluster =
      cluster_C1 !=
      cluster_C2
  ) |>
  left_join(
    clustering_base,
    by = "id"
  ) |>
  arrange(
    id
  )


if (
  nrow(
    clustering_C1_C2_assignments
  ) !=
  n_expected_units
) {
  stop(
    "Final C1/C2 assignment table does not contain ",
    n_expected_units,
    " units."
  )
}


if (
  anyNA(
    clustering_C1_C2_assignments$cluster_C2
  )
) {
  stop(
    "One or more C2 labels could not be aligned."
  )
}


# ------------------------------------------------------------
# 11. C1/C2 stability
# ------------------------------------------------------------

n_stable <- sum(
  !clustering_C1_C2_assignments$
    changed_cluster
)

n_reassigned <- sum(
  clustering_C1_C2_assignments$
    changed_cluster
)

percentage_stable <-
  100 *
  n_stable /
  n_expected_units

percentage_reassigned <-
  100 *
  n_reassigned /
  n_expected_units


clustering_C1_C2_summary <- data.frame(
  
  k =
    k_final,
  
  n_units =
    n_expected_units,
  
  ARI =
    ARI_C1_C2,
  
  n_stable =
    n_stable,
  
  percentage_stable =
    percentage_stable,
  
  n_reassigned =
    n_reassigned,
  
  percentage_reassigned =
    percentage_reassigned,
  
  silhouette_mean_C1 =
    k12_evaluation$
    silhouette_media[
      k12_evaluation$escenario ==
        "C1"
    ],
  
  silhouette_mean_C2 =
    k12_evaluation$
    silhouette_media[
      k12_evaluation$escenario ==
        "C2"
    ],
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 12. Aligned transition matrix
# ------------------------------------------------------------

aligned_transition_table <- table(
  
  C1 =
    clustering_C1_C2_assignments$
    cluster_C1,
  
  C2 =
    clustering_C1_C2_assignments$
    cluster_C2
)


clustering_C1_C2_transition_matrix <-
  as.data.frame.matrix(
    aligned_transition_table
  )


clustering_C1_C2_transition_matrix <-
  cbind(
    
    cluster_C1 =
      as.integer(
        rownames(
          clustering_C1_C2_transition_matrix
        )
      ),
    
    clustering_C1_C2_transition_matrix
  )


rownames(
  clustering_C1_C2_transition_matrix
) <- NULL


# ------------------------------------------------------------
# 13. Real C1 -> C2 transitions
# ------------------------------------------------------------

clustering_C1_C2_transitions <-
  clustering_C1_C2_assignments |>
  filter(
    changed_cluster
  ) |>
  count(
    cluster_C1,
    cluster_C2,
    name = "n"
  ) |>
  arrange(
    desc(n),
    cluster_C1,
    cluster_C2
  )


if (
  sum(
    clustering_C1_C2_transitions$n
  ) !=
  n_reassigned
) {
  stop(
    "Transition counts do not equal the number of reassigned units."
  )
}


# ------------------------------------------------------------
# 14. Reassigned units
# ------------------------------------------------------------

clustering_C1_C2_reassigned_units <-
  clustering_C1_C2_assignments |>
  filter(
    changed_cluster
  ) |>
  select(
    id,
    cluster_C1,
    cluster_C2_original,
    cluster_C2,
    silhouette_C1,
    silhouette_C2,
    USyV,
    Suelos,
    Grad_Pend,
    P_sum,
    Esc_sum,
    Agua_D,
    Agua_D_log
  ) |>
  arrange(
    cluster_C1,
    cluster_C2,
    id
  )


if (
  nrow(
    clustering_C1_C2_reassigned_units
  ) !=
  n_reassigned
) {
  stop(
    "Reassigned-unit table has an unexpected number of rows."
  )
}


# ------------------------------------------------------------
# 15. Modal-category helper
# ------------------------------------------------------------

mode_character <- function(x) {
  
  x <- x[
    !is.na(x)
  ]
  
  
  if (
    length(x) == 0L
  ) {
    return(
      NA_character_
    )
  }
  
  
  tab <- table(
    x
  )
  
  maximum_n <- max(
    tab
  )
  
  
  candidates <- names(
    tab
  )[
    tab ==
      maximum_n
  ]
  
  
  # If modal categories are tied, preserve all tied categories.
  
  paste(
    candidates,
    collapse = " / "
  )
}


# ------------------------------------------------------------
# 16. C1 cluster profiles
# ------------------------------------------------------------

clustering_C1_profiles <-
  clustering_C1_C2_assignments |>
  group_by(
    cluster =
      cluster_C1
  ) |>
  summarise(
    
    n =
      n(),
    
    USyV_mode =
      mode_character(
        USyV
      ),
    
    soil_mode =
      mode_character(
        Suelos
      ),
    
    slope_mode =
      mode_character(
        Grad_Pend
      ),
    
    n_P_sum_available =
      sum(
        !is.na(
          P_sum
        )
      ),
    
    P_sum_median =
      median(
        P_sum,
        na.rm = TRUE
      ),
    
    Esc_sum_median =
      median(
        Esc_sum,
        na.rm = TRUE
      ),
    
    silhouette_median =
      median(
        silhouette_C1,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) |>
  arrange(
    cluster
  )


# ------------------------------------------------------------
# 17. C2 cluster profiles
# ------------------------------------------------------------
#
# Agua_D is reported in the original scale for interpretation,
# although clustering uses log1p(Agua_D).
#
# ------------------------------------------------------------

clustering_C2_profiles <-
  clustering_C1_C2_assignments |>
  group_by(
    cluster =
      cluster_C2
  ) |>
  summarise(
    
    n =
      n(),
    
    USyV_mode =
      mode_character(
        USyV
      ),
    
    soil_mode =
      mode_character(
        Suelos
      ),
    
    slope_mode =
      mode_character(
        Grad_Pend
      ),
    
    n_P_sum_available =
      sum(
        !is.na(
          P_sum
        )
      ),
    
    P_sum_median =
      median(
        P_sum,
        na.rm = TRUE
      ),
    
    Esc_sum_median =
      median(
        Esc_sum,
        na.rm = TRUE
      ),
    
    Agua_D_median =
      median(
        Agua_D,
        na.rm = TRUE
      ),
    
    Agua_D_log_median =
      median(
        Agua_D_log,
        na.rm = TRUE
      ),
    
    silhouette_median =
      median(
        silhouette_C2,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) |>
  arrange(
    cluster
  )


# ------------------------------------------------------------
# 18. Validate cluster profiles
# ------------------------------------------------------------

if (
  nrow(
    clustering_C1_profiles
  ) !=
  k_final ||
  nrow(
    clustering_C2_profiles
  ) !=
  k_final
) {
  stop(
    "Expected 12 cluster profiles for C1 and C2."
  )
}


if (
  sum(
    clustering_C1_profiles$n
  ) !=
  n_expected_units ||
  sum(
    clustering_C2_profiles$n
  ) !=
  n_expected_units
) {
  stop(
    "Cluster-profile sizes do not sum to the total number of units."
  )
}


# ------------------------------------------------------------
# 19. Comparable C1/C2 profile table
# ------------------------------------------------------------

C1_profile_long <-
  clustering_C1_profiles |>
  mutate(
    
    escenario =
      "C1",
    
    Agua_D_median =
      NA_real_,
    
    Agua_D_log_median =
      NA_real_
  ) |>
  select(
    escenario,
    cluster,
    n,
    USyV_mode,
    soil_mode,
    slope_mode,
    n_P_sum_available,
    P_sum_median,
    Esc_sum_median,
    Agua_D_median,
    Agua_D_log_median,
    silhouette_median
  )


C2_profile_long <-
  clustering_C2_profiles |>
  mutate(
    escenario =
      "C2"
  ) |>
  select(
    escenario,
    cluster,
    n,
    USyV_mode,
    soil_mode,
    slope_mode,
    n_P_sum_available,
    P_sum_median,
    Esc_sum_median,
    Agua_D_median,
    Agua_D_log_median,
    silhouette_median
  )


clustering_profiles_comparable <-
  bind_rows(
    C1_profile_long,
    C2_profile_long
  ) |>
  arrange(
    cluster,
    escenario
  )


# ------------------------------------------------------------
# 20. Output paths
# ------------------------------------------------------------

assignments_output_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_assignments.csv"
)

summary_output_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_summary.csv"
)

label_output_file <- here::here(
  "results",
  "tables",
  "clustering_C2_to_C1_labels.csv"
)

transition_matrix_output_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_transition_matrix.csv"
)

transitions_output_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_transitions.csv"
)

reassigned_output_file <- here::here(
  "results",
  "tables",
  "clustering_C1_C2_reassigned_units.csv"
)

C1_profiles_output_file <- here::here(
  "results",
  "tables",
  "clustering_C1_profiles.csv"
)

C2_profiles_output_file <- here::here(
  "results",
  "tables",
  "clustering_C2_profiles.csv"
)

comparable_profiles_output_file <- here::here(
  "results",
  "tables",
  "clustering_profiles_comparable.csv"
)


# ------------------------------------------------------------
# 21. Write outputs
# ------------------------------------------------------------

write.csv(
  clustering_C1_C2_assignments |>
    select(
      id,
      cluster_C1,
      cluster_C2_original,
      cluster_C2,
      changed_cluster,
      silhouette_C1,
      silhouette_C2
    ),
  assignments_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_C1_C2_summary,
  summary_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  label_correspondence,
  label_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_C1_C2_transition_matrix,
  transition_matrix_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_C1_C2_transitions,
  transitions_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_C1_C2_reassigned_units,
  reassigned_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_C1_profiles,
  C1_profiles_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_C2_profiles,
  C2_profiles_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


write.csv(
  clustering_profiles_comparable,
  comparable_profiles_output_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 22. Read-back validation
# ------------------------------------------------------------

assignments_check <- read.csv(
  assignments_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_check <- read.csv(
  summary_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

labels_check <- read.csv(
  label_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

reassigned_check <- read.csv(
  reassigned_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

C1_profiles_check <- read.csv(
  C1_profiles_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

C2_profiles_check <- read.csv(
  C2_profiles_output_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (
  nrow(
    assignments_check
  ) !=
  n_expected_units
) {
  stop(
    "Final clustering assignments failed read-back validation."
  )
}


if (
  nrow(
    summary_check
  ) != 1L
) {
  stop(
    "C1/C2 summary failed read-back validation."
  )
}


if (
  nrow(
    labels_check
  ) !=
  k_final
) {
  stop(
    "Label correspondence failed read-back validation."
  )
}


if (
  nrow(
    reassigned_check
  ) !=
  n_reassigned
) {
  stop(
    "Reassigned-unit output failed read-back validation."
  )
}


if (
  nrow(
    C1_profiles_check
  ) !=
  k_final ||
  nrow(
    C2_profiles_check
  ) !=
  k_final
) {
  stop(
    "Cluster profiles failed read-back validation."
  )
}


# ------------------------------------------------------------
# 23. Console results
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "FINAL C1 vs C2 COMPARISON — k = 12\n"
)

cat(
  "============================================\n"
)

print(
  clustering_C1_C2_summary,
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "C2 -> C1 LABEL CORRESPONDENCE\n"
)

cat(
  "============================================\n"
)

print(
  label_correspondence,
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "REAL C1 -> C2 TRANSITIONS\n"
)

cat(
  "============================================\n"
)


if (
  nrow(
    clustering_C1_C2_transitions
  ) == 0L
) {
  
  cat(
    "No units were reassigned.\n"
  )
  
} else {
  
  print(
    clustering_C1_C2_transitions,
    row.names = FALSE
  )
}


cat(
  "\n============================================\n"
)

cat(
  "C1 CLUSTER PROFILES\n"
)

cat(
  "============================================\n"
)

print(
  as.data.frame(
    clustering_C1_profiles
  ),
  row.names = FALSE
)


cat(
  "\n============================================\n"
)

cat(
  "C2 CLUSTER PROFILES\n"
)

cat(
  "============================================\n"
)

print(
  as.data.frame(
    clustering_C2_profiles
  ),
  row.names = FALSE
)


cat(
  "\n============================================\n\n"
)


# ------------------------------------------------------------
# 24. Completion
# ------------------------------------------------------------

message(
  "11_compare_C1_C2.R completed successfully."
)

message(
  "k used: ",
  k_final,
  " for both configurations."
)

message(
  "ARI = ",
  round(
    ARI_C1_C2,
    6
  )
)

message(
  "Stable units = ",
  n_stable,
  " / ",
  n_expected_units,
  " (",
  round(
    percentage_stable,
    2
  ),
  "%)"
)

message(
  "Reassigned units = ",
  n_reassigned,
  " / ",
  n_expected_units,
  " (",
  round(
    percentage_reassigned,
    2
  ),
  "%)"
)

message(
  "Next: R/12_compare_C2_UREH.R"
)