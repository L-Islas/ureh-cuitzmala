# ============================================================
# 10. GOWER DISTANCE + PAM CLUSTERING
# ============================================================
#
# Purpose
# -------
# Evaluate multivariate similarity for the two independent
# clustering configurations prepared in R/09_prepare_clustering.R.
#
# C1
# --
# Local bio-physical similarity:
#
#   USyV
#   Suelos
#   Grad_Pend
#   P_sum
#   Esc_sum
#
# C2
# --
# C1 + upstream contribution:
#
#   Agua_D_log = log1p(Agua_D)
#
# Procedure
# ---------
# 1. Calculate Gower distances.
# 2. Run PAM for k = 2,...,20.
# 3. Calculate mean and median silhouette.
# 4. Calculate proportion of negative silhouettes.
# 5. Record minimum, median and maximum cluster sizes.
# 6. Preserve all unit-level assignments.
#
# IMPORTANT
# ---------
# This script DOES NOT automatically select an optimal k.
#
# The maximum observed silhouette is reported only as a
# diagnostic. If the maximum occurs at the upper boundary
# (k = 20), it must not automatically be interpreted as an
# optimum.
#
# Inputs
# ------
# temp/clustering_C1.rds
# temp/clustering_C2.rds
#
# Intermediate outputs
# --------------------
# temp/gower_C1.rds
# temp/gower_C2.rds
#
# Final transparent outputs
# -------------------------
# results/tables/clustering_k_evaluation.csv
# results/tables/clustering_k_maxima.csv
# results/tables/clustering_assignments_k2_20.csv
#
# ============================================================


# ------------------------------------------------------------
# 0. Configuration and packages
# ------------------------------------------------------------

source("R/00_config.R")

library(cluster)
library(dplyr)

results_tables_dir <- here::here(
  "results",
  "tables"
)

temp_dir <- here::here(
  "temp"
)

dir.create(
  results_tables_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  temp_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Input paths
# ------------------------------------------------------------

C1_file <- here::here(
  "temp",
  "clustering_C1.rds"
)

C2_file <- here::here(
  "temp",
  "clustering_C2.rds"
)

required_files <- c(
  C1_file,
  C2_file
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
    "Required clustering inputs not found:\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    "\nRun R/09_prepare_clustering.R first."
  )
}


# ------------------------------------------------------------
# 2. Read prepared clustering data
# ------------------------------------------------------------

C1 <- readRDS(
  C1_file
)

C2 <- readRDS(
  C2_file
)


# ------------------------------------------------------------
# 3. General validation
# ------------------------------------------------------------

if (
  nrow(C1) != n_expected_units ||
  nrow(C2) != n_expected_units
) {
  stop(
    "C1 and C2 must each contain ",
    n_expected_units,
    " units."
  )
}


if (
  anyNA(C1$id) ||
  anyNA(C2$id) ||
  anyDuplicated(C1$id) ||
  anyDuplicated(C2$id)
) {
  stop(
    "Invalid ids in C1 or C2."
  )
}


if (!identical(
  as.integer(C1$id),
  as.integer(C2$id)
)) {
  stop(
    "C1 and C2 ids are not identical or are not in the same order."
  )
}


if (!identical(
  sort(
    as.integer(
      C1$id
    )
  ),
  seq_len(
    n_expected_units
  )
)) {
  stop(
    "Unexpected ids in clustering inputs."
  )
}


# ------------------------------------------------------------
# 4. Validate exact variable definitions
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
  names(C1),
  expected_C1_fields
)) {
  stop(
    "Unexpected C1 fields."
  )
}


if (!identical(
  names(C2),
  expected_C2_fields
)) {
  stop(
    "Unexpected C2 fields."
  )
}


# ------------------------------------------------------------
# 5. Validate variable classes
# ------------------------------------------------------------

if (
  !is.factor(C1$USyV) ||
  is.ordered(C1$USyV)
) {
  stop(
    "USyV must be a nominal factor."
  )
}


if (
  !is.factor(C1$Suelos) ||
  is.ordered(C1$Suelos)
) {
  stop(
    "Suelos must be a nominal factor."
  )
}


if (!is.ordered(
  C1$Grad_Pend
)) {
  stop(
    "Grad_Pend must be an ordered factor."
  )
}


if (
  !is.numeric(C1$P_sum) ||
  !is.numeric(C1$Esc_sum)
) {
  stop(
    "P_sum and Esc_sum must be numeric."
  )
}


if (!is.numeric(
  C2$Agua_D_log
)) {
  stop(
    "Agua_D_log must be numeric."
  )
}


# ------------------------------------------------------------
# 6. Validate missing-value structure
# ------------------------------------------------------------

C1_NA <- sapply(
  C1[
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

C2_NA <- sapply(
  C2[
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


allowed_NA_variable <- "P_sum"

if (
  any(
    C1_NA[
      names(C1_NA) !=
      allowed_NA_variable
    ] > 0
  )
) {
  stop(
    "C1 contains missing values outside P_sum."
  )
}


if (
  any(
    C2_NA[
      names(C2_NA) !=
      allowed_NA_variable
    ] > 0
  )
) {
  stop(
    "C2 contains missing values outside P_sum."
  )
}


if (
  C1_NA["P_sum"] !=
  C2_NA["P_sum"]
) {
  stop(
    "C1 and C2 have different P_sum missing-value patterns."
  )
}


if (!identical(
  which(
    is.na(
      C1$P_sum
    )
  ),
  which(
    is.na(
      C2$P_sum
    )
  )
)) {
  stop(
    "P_sum missing units differ between C1 and C2."
  )
}


message(
  "Clustering inputs validated."
)

message(
  "P_sum missing values retained: ",
  C1_NA["P_sum"],
  "."
)


# ------------------------------------------------------------
# 7. Prepare matrices/data frames for Gower
# ------------------------------------------------------------
#
# id does not participate in the distance.
#
# Every variable receives the default equal Gower weight.
#
# Missing P_sum values are retained. cluster::daisy() computes
# each pairwise distance using the variables available for that
# pair.
#
# ------------------------------------------------------------

X_C1 <- C1 |>
  select(
    -id
  )

X_C2 <- C2 |>
  select(
    -id
  )


# Preserve ids as row labels for diagnostic traceability.

rownames(X_C1) <- as.character(
  C1$id
)

rownames(X_C2) <- as.character(
  C2$id
)


# ------------------------------------------------------------
# 8. Calculate Gower distance
# ------------------------------------------------------------

message("")
message(
  "============================================"
)
message(
  "CALCULATING GOWER DISTANCE: C1"
)
message(
  "============================================"
)

gower_C1 <- cluster::daisy(
  X_C1,
  metric = "gower"
)


message("")
message(
  "============================================"
)
message(
  "CALCULATING GOWER DISTANCE: C2"
)
message(
  "============================================"
)

gower_C2 <- cluster::daisy(
  X_C2,
  metric = "gower"
)


# ------------------------------------------------------------
# 9. Validate Gower distances
# ------------------------------------------------------------

validate_gower <- function(
    distance_object,
    object_name
) {
  
  distance_values <- as.vector(
    distance_object
  )
  
  
  if (
    anyNA(
      distance_values
    ) ||
    any(
      !is.finite(
        distance_values
      )
    )
  ) {
    stop(
      object_name,
      " contains NA or non-finite distances."
    )
  }
  
  
  if (
    any(
      distance_values <
      -1e-12
    ) ||
    any(
      distance_values >
      1 + 1e-12
    )
  ) {
    stop(
      object_name,
      " contains distances outside [0,1]."
    )
  }
  
  
  if (
    length(
      attr(
        distance_object,
        "Labels"
      )
    ) !=
    n_expected_units
  ) {
    stop(
      object_name,
      " does not contain the expected labels."
    )
  }
  
  
  invisible(TRUE)
}


validate_gower(
  gower_C1,
  "Gower C1"
)

validate_gower(
  gower_C2,
  "Gower C2"
)


message(
  "Gower distances calculated and validated."
)


# ------------------------------------------------------------
# 10. Save distance objects
# ------------------------------------------------------------

gower_C1_file <- here::here(
  "temp",
  "gower_C1.rds"
)

gower_C2_file <- here::here(
  "temp",
  "gower_C2.rds"
)


saveRDS(
  gower_C1,
  gower_C1_file
)

saveRDS(
  gower_C2,
  gower_C2_file
)


# ------------------------------------------------------------
# 11. PAM evaluation function
# ------------------------------------------------------------

evaluate_PAM <- function(
    distance_object,
    ids,
    scenario,
    k_values = 2:20
) {
  
  summary_list <- vector(
    "list",
    length(
      k_values
    )
  )
  
  assignment_list <- vector(
    "list",
    length(
      k_values
    )
  )
  
  names(
    summary_list
  ) <- as.character(
    k_values
  )
  
  names(
    assignment_list
  ) <- as.character(
    k_values
  )
  
  
  for (
    k in k_values
  ) {
    
    message(
      scenario,
      " | PAM k = ",
      k
    )
    
    
    # --------------------------------------------------------
    # 11.1 PAM
    # --------------------------------------------------------
    
    fit <- cluster::pam(
      x =
        distance_object,
      k =
        k,
      diss =
        TRUE,
      keep.diss =
        FALSE,
      keep.data =
        FALSE
    )
    
    
    cluster_assignment <- as.integer(
      fit$clustering
    )
    
    
    if (
      length(
        cluster_assignment
      ) !=
      length(ids)
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": unexpected number of PAM assignments."
      )
    }
    
    
    if (
      anyNA(
        cluster_assignment
      )
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": PAM generated missing assignments."
      )
    }
    
    
    if (
      length(
        unique(
          cluster_assignment
        )
      ) !=
      k
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": PAM did not generate exactly k clusters."
      )
    }
    
    
    # --------------------------------------------------------
    # 11.2 Silhouette
    # --------------------------------------------------------
    
    silhouette_object <- cluster::silhouette(
      cluster_assignment,
      distance_object
    )
    
    
    silhouette_df <- as.data.frame(
      silhouette_object
    )
    
    
    if (
      nrow(
        silhouette_df
      ) !=
      length(ids)
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": unexpected number of silhouette values."
      )
    }
    
    
    if (
      anyNA(
        silhouette_df$sil_width
      ) ||
      any(
        !is.finite(
          silhouette_df$sil_width
        )
      )
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": invalid silhouette values."
      )
    }
    
    
    if (
      any(
        silhouette_df$sil_width <
        -1 - 1e-12
      ) ||
      any(
        silhouette_df$sil_width >
        1 + 1e-12
      )
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": silhouette values outside [-1,1]."
      )
    }
    
    
    # --------------------------------------------------------
    # 11.3 Cluster sizes
    # --------------------------------------------------------
    
    cluster_sizes <- table(
      cluster_assignment
    )
    
    
    if (
      sum(
        cluster_sizes
      ) !=
      length(ids)
    ) {
      stop(
        scenario,
        " k=",
        k,
        ": cluster sizes do not sum to the expected number of units."
      )
    }
    
    
    # --------------------------------------------------------
    # 11.4 Summary for k
    # --------------------------------------------------------
    
    summary_list[[as.character(k)]] <-
      data.frame(
        
        escenario =
          scenario,
        
        k =
          k,
        
        silhouette_media =
          mean(
            silhouette_df$sil_width
          ),
        
        silhouette_mediana =
          median(
            silhouette_df$sil_width
          ),
        
        proporcion_silhouette_negativa =
          mean(
            silhouette_df$sil_width < 0
          ),
        
        cluster_min_n =
          min(
            cluster_sizes
          ),
        
        cluster_mediana_n =
          median(
            as.numeric(
              cluster_sizes
            )
          ),
        
        cluster_max_n =
          max(
            cluster_sizes
          ),
        
        stringsAsFactors = FALSE
      )
    
    
    # --------------------------------------------------------
    # 11.5 Unit-level assignment for k
    # --------------------------------------------------------
    
    assignment_list[[as.character(k)]] <-
      data.frame(
        
        id =
          ids,
        
        escenario =
          scenario,
        
        k =
          k,
        
        cluster =
          cluster_assignment,
        
        silhouette =
          silhouette_df$sil_width,
        
        stringsAsFactors = FALSE
      )
  }
  
  
  list(
    
    summary =
      bind_rows(
        summary_list
      ),
    
    assignments =
      bind_rows(
        assignment_list
      )
  )
}


# ------------------------------------------------------------
# 12. Evaluate C1: k = 2,...,20
# ------------------------------------------------------------

message("")
message(
  "============================================"
)
message(
  "EVALUATING PAM: C1"
)
message(
  "============================================"
)

PAM_C1 <- evaluate_PAM(
  distance_object =
    gower_C1,
  
  ids =
    C1$id,
  
  scenario =
    "C1",
  
  k_values =
    2:20
)


# ------------------------------------------------------------
# 13. Evaluate C2: k = 2,...,20
# ------------------------------------------------------------

message("")
message(
  "============================================"
)
message(
  "EVALUATING PAM: C2"
)
message(
  "============================================"
)

PAM_C2 <- evaluate_PAM(
  distance_object =
    gower_C2,
  
  ids =
    C2$id,
  
  scenario =
    "C2",
  
  k_values =
    2:20
)


# ------------------------------------------------------------
# 14. Combine k-evaluation results
# ------------------------------------------------------------

clustering_k_evaluation <-
  bind_rows(
    PAM_C1$summary,
    PAM_C2$summary
  ) |>
  arrange(
    escenario,
    k
  ) |>
  group_by(
    escenario
  ) |>
  mutate(
    
    incremento_silhouette =
      silhouette_media -
      lag(
        silhouette_media
      )
    
  ) |>
  ungroup()


# ------------------------------------------------------------
# 15. Combine all unit-level assignments
# ------------------------------------------------------------

clustering_assignments <-
  bind_rows(
    PAM_C1$assignments,
    PAM_C2$assignments
  ) |>
  arrange(
    escenario,
    k,
    id
  )


# ------------------------------------------------------------
# 16. Validate combined results
# ------------------------------------------------------------

expected_k_values <- 2:20

expected_evaluation_rows <-
  length(
    expected_k_values
  ) *
  2L

expected_assignment_rows <-
  n_expected_units *
  length(
    expected_k_values
  ) *
  2L


if (
  nrow(
    clustering_k_evaluation
  ) !=
  expected_evaluation_rows
) {
  stop(
    "Unexpected number of rows in clustering k evaluation."
  )
}


if (
  nrow(
    clustering_assignments
  ) !=
  expected_assignment_rows
) {
  stop(
    "Unexpected number of clustering assignment rows."
  )
}


if (
  anyDuplicated(
    clustering_assignments[
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


assignment_counts <-
  clustering_assignments |>
  count(
    escenario,
    k,
    name =
      "n"
  )


if (
  any(
    assignment_counts$n !=
    n_expected_units
  )
) {
  stop(
    "One or more scenario-k combinations do not contain ",
    n_expected_units,
    " assignments."
  )
}


# ------------------------------------------------------------
# 17. Maximum observed silhouette
# ------------------------------------------------------------
#
# This is a diagnostic only.
#
# A maximum at k = 20 lies on the evaluated upper boundary and
# therefore does not automatically define an optimal k.
#
# ------------------------------------------------------------

clustering_k_maxima <-
  clustering_k_evaluation |>
  group_by(
    escenario
  ) |>
  slice_max(
    order_by =
      silhouette_media,
    n =
      1,
    with_ties =
      FALSE
  ) |>
  transmute(
    
    escenario,
    
    k_max_silhouette =
      k,
    
    silhouette_maxima =
      silhouette_media,
    
    silhouette_mediana =
      silhouette_mediana,
    
    proporcion_silhouette_negativa =
      proporcion_silhouette_negativa,
    
    cluster_min_n =
      cluster_min_n,
    
    cluster_mediana_n =
      cluster_mediana_n,
    
    cluster_max_n =
      cluster_max_n,
    
    maximo_en_limite_superior =
      k == max(
        expected_k_values
      )
    
  ) |>
  ungroup()


# ------------------------------------------------------------
# 18. Extract k = 12 diagnostic preview
# ------------------------------------------------------------
#
# k = 12 is NOT selected by this script.
#
# It is displayed here because it is the finalized solution used
# in the subsequent C1/C2 comparison.
#
# ------------------------------------------------------------

k12_preview <-
  clustering_k_evaluation |>
  filter(
    k == 12
  ) |>
  select(
    escenario,
    k,
    silhouette_media,
    silhouette_mediana,
    proporcion_silhouette_negativa,
    cluster_min_n,
    cluster_mediana_n,
    cluster_max_n
  ) |>
  arrange(
    escenario
  )


# ------------------------------------------------------------
# 19. Output paths
# ------------------------------------------------------------

evaluation_file <- here::here(
  "results",
  "tables",
  "clustering_k_evaluation.csv"
)

maxima_file <- here::here(
  "results",
  "tables",
  "clustering_k_maxima.csv"
)

assignments_file <- here::here(
  "results",
  "tables",
  "clustering_assignments_k2_20.csv"
)


# ------------------------------------------------------------
# 20. Write outputs
# ------------------------------------------------------------

write.csv(
  clustering_k_evaluation,
  evaluation_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  clustering_k_maxima,
  maxima_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  clustering_assignments,
  assignments_file,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------
# 21. Read-back validation
# ------------------------------------------------------------

evaluation_check <- read.csv(
  evaluation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

maxima_check <- read.csv(
  maxima_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

assignments_check <- read.csv(
  assignments_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


if (
  nrow(
    evaluation_check
  ) !=
  expected_evaluation_rows
) {
  stop(
    "Clustering evaluation failed read-back validation."
  )
}


if (
  nrow(
    maxima_check
  ) != 2L
) {
  stop(
    "Clustering maxima table failed read-back validation."
  )
}


if (
  nrow(
    assignments_check
  ) !=
  expected_assignment_rows
) {
  stop(
    "Clustering assignments failed read-back validation."
  )
}


if (
  anyDuplicated(
    assignments_check[
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
    "Duplicated clustering assignments after read-back."
  )
}


# ------------------------------------------------------------
# 22. Console summary
# ------------------------------------------------------------

cat(
  "\n============================================\n"
)

cat(
  "GOWER + PAM EVALUATION COMPLETE\n"
)

cat(
  "============================================\n"
)

cat(
  "Units per configuration:",
  n_expected_units,
  "\n"
)

cat(
  "Evaluated k range: 2-20\n"
)

cat(
  "Total scenario-k solutions:",
  nrow(
    clustering_k_evaluation
  ),
  "\n\n"
)


cat(
  "k = 12 diagnostic:\n"
)

print(
  as.data.frame(
    k12_preview
  ),
  row.names = FALSE
)


cat(
  "\nMaximum observed mean silhouette:\n"
)

print(
  as.data.frame(
    clustering_k_maxima
  ),
  row.names = FALSE
)


if (
  any(
    clustering_k_maxima$
    maximo_en_limite_superior
  )
) {
  
  boundary_scenarios <- paste(
    clustering_k_maxima$escenario[
      clustering_k_maxima$
        maximo_en_limite_superior
    ],
    collapse = ", "
  )
  
  cat(
    "\nNOTE: maximum silhouette occurs at k = 20 for: ",
    boundary_scenarios,
    ".\n",
    sep = ""
  )
  
  cat(
    "This boundary maximum is not automatically interpreted as an optimal k.\n"
  )
}


cat(
  "\nNo automatic k selection was performed.\n"
)

cat(
  "============================================\n\n"
)


# ------------------------------------------------------------
# 23. Completion
# ------------------------------------------------------------

message(
  "10_gower_PAM.R completed successfully."
)

message(
  "Generated: temp/gower_C1.rds"
)

message(
  "Generated: temp/gower_C2.rds"
)

message(
  "Generated: results/tables/clustering_k_evaluation.csv"
)

message(
  "Generated: results/tables/clustering_k_maxima.csv"
)

message(
  "Generated: results/tables/clustering_assignments_k2_20.csv"
)

message(
  "Next: R/11_compare_C1_C2.R"
)