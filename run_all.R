# ============================================================
# RUN COMPLETE UREH REPRODUCIBILITY WORKFLOW
# ============================================================
#
# Purpose
# -------
# Reproduce the complete analytical workflow from the
# analysis-ready input data through the final spatial products,
# followed by regression/reproducibility checks.
#
# Publication figures are intentionally outside this automated
# workflow because their final composition may involve R, QGIS
# and manual graphical editing.
#
# ============================================================


run_ureh_workflow <- function() {
  
  
  # ----------------------------------------------------------
  # 0. Confirm repository root
  # ----------------------------------------------------------
  
  if (
    !file.exists("ureh-cuitzmala.Rproj") ||
    !file.exists("R/00_config.R")
  ) {
    stop(
      paste0(
        "run_all.R must be executed from the root of the ",
        "ureh-cuitzmala repository."
      )
    )
  }
  
  
  workflow_start <- Sys.time()
  
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    "UREH REPRODUCIBILITY WORKFLOW\n"
  )
  
  cat(
    "============================================\n"
  )
  
  cat(
    "Started: ",
    format(
      workflow_start,
      "%Y-%m-%d %H:%M:%S"
    ),
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # 1. Workflow scripts
  # ----------------------------------------------------------
  
  workflow_scripts <- c(
    
    "R/01_prepare_data.R",
    
    "R/02_soft_runoff_evidence.R",
    
    "R/03_build_bayesian_network.R",
    
    "R/04_infer_M0.R",
    
    "R/05_summarize_M0.R",
    
    "R/06_sensitivity_analysis.R",
    
    "R/07_summarize_sensitivity.R",
    
    "R/08_TM_diagnostic.R",
    
    "R/09_prepare_clustering.R",
    
    "R/10_gower_PAM.R",
    
    "R/11_compare_C1_C2.R",
    
    "R/12_compare_C2_UREH.R",
    
    "R/13_build_spatial_results.R",
    
    "checks/check_reproducibility.R"
  )
  
  
  # ----------------------------------------------------------
  # 2. Confirm all scripts exist
  # ----------------------------------------------------------
  
  missing_scripts <- workflow_scripts[
    !file.exists(
      workflow_scripts
    )
  ]
  
  
  if (
    length(
      missing_scripts
    ) > 0L
  ) {
    
    stop(
      "Missing workflow scripts:\n",
      paste(
        missing_scripts,
        collapse = "\n"
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 3. Execution log
  # ----------------------------------------------------------
  
  script_log <- data.frame(
    
    step =
      seq_along(
        workflow_scripts
      ),
    
    script =
      workflow_scripts,
    
    status =
      rep(
        NA_character_,
        length(
          workflow_scripts
        )
      ),
    
    elapsed_seconds =
      rep(
        NA_real_,
        length(
          workflow_scripts
        )
      ),
    
    stringsAsFactors =
      FALSE
  )
  
  
  # ----------------------------------------------------------
  # 4. Execute workflow
  # ----------------------------------------------------------
  #
  # Each script gets its own environment.
  #
  # This prevents variables created inside analytical scripts
  # (for example i, j, x, result, etc.) from overwriting the
  # workflow controller.
  #
  # Analytical communication between scripts occurs through
  # explicitly generated files, not through hidden objects left
  # in the R global environment.
  #
  # ----------------------------------------------------------
  
  for (
    step_index in seq_along(
      workflow_scripts
    )
  ) {
    
    script_path <- workflow_scripts[
      step_index
    ]
    
    
    cat(
      "\n\n============================================================\n"
    )
    
    cat(
      "[",
      step_index,
      "/",
      length(
        workflow_scripts
      ),
      "] ",
      script_path,
      "\n",
      sep = ""
    )
    
    cat(
      "============================================================\n\n"
    )
    
    
    script_start <- Sys.time()
    
    
    # --------------------------------------------------------
    # Fresh isolated environment for this script
    # --------------------------------------------------------
    
    script_environment <- new.env(
      parent =
        globalenv()
    )
    
    
    execution_error <- NULL
    
    
    tryCatch(
      
      {
        
        source(
          script_path,
          local =
            script_environment,
          echo =
            FALSE,
          chdir =
            FALSE
        )
        
      },
      
      error = function(
    e
      ) {
        
        execution_error <<- e
      }
    )
    
    
    script_end <- Sys.time()
    
    
    elapsed_seconds <- as.numeric(
      difftime(
        script_end,
        script_start,
        units =
          "secs"
      )
    )
    
    
    script_log$elapsed_seconds[
      step_index
    ] <- elapsed_seconds
    
    
    # --------------------------------------------------------
    # Handle failure
    # --------------------------------------------------------
    
    if (
      !is.null(
        execution_error
      )
    ) {
      
      script_log$status[
        step_index
      ] <- "FAIL"
      
      
      cat(
        "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      )
      
      cat(
        "WORKFLOW FAILED\n"
      )
      
      cat(
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
      )
      
      cat(
        "Step: ",
        step_index,
        " / ",
        length(
          workflow_scripts
        ),
        "\n",
        sep = ""
      )
      
      cat(
        "Script: ",
        script_path,
        "\n",
        sep = ""
      )
      
      cat(
        "Error: ",
        conditionMessage(
          execution_error
        ),
        "\n",
        sep = ""
      )
      
      cat(
        "Elapsed: ",
        round(
          elapsed_seconds,
          2
        ),
        " s\n",
        sep = ""
      )
      
      
      cat(
        "\nWorkflow log up to failure:\n\n"
      )
      
      
      print(
        script_log[
          seq_len(
            step_index
          ),
        ],
        row.names =
          FALSE
      )
      
      
      stop(
        paste0(
          "UREH reproducibility workflow failed at ",
          script_path,
          "."
        ),
        call. =
          FALSE
      )
    }
    
    
    # --------------------------------------------------------
    # Successful step
    # --------------------------------------------------------
    
    script_log$status[
      step_index
    ] <- "PASS"
    
    
    cat(
      "\nCompleted: ",
      script_path,
      "\n",
      sep = ""
    )
    
    
    cat(
      "Elapsed: ",
      round(
        elapsed_seconds,
        2
      ),
      " s\n",
      sep = ""
    )
    
    
    # Explicitly discard the script environment before moving
    # to the next analytical step.
    
    rm(
      script_environment
    )
  }
  
  
  # ----------------------------------------------------------
  # 5. Verify reproducibility report
  # ----------------------------------------------------------
  
  reproducibility_report_file <- file.path(
    "checks",
    "reproducibility_report.csv"
  )
  
  
  if (
    !file.exists(
      reproducibility_report_file
    )
  ) {
    
    stop(
      paste0(
        "Workflow finished but ",
        "checks/reproducibility_report.csv was not generated."
      )
    )
  }
  
  
  reproducibility_report <- read.csv(
    reproducibility_report_file,
    stringsAsFactors =
      FALSE,
    check.names =
      FALSE
  )
  
  
  if (!(
    "status" %in%
    names(
      reproducibility_report
    )
  )) {
    
    stop(
      paste0(
        "Reproducibility report does not contain ",
        "the expected 'status' field."
      )
    )
  }
  
  
  n_controls <- nrow(
    reproducibility_report
  )
  
  
  n_controls_pass <- sum(
    reproducibility_report$status ==
      "PASS"
  )
  
  
  n_controls_fail <- sum(
    reproducibility_report$status !=
      "PASS"
  )
  
  
  if (
    n_controls != 44L ||
    n_controls_fail != 0L
  ) {
    
    stop(
      paste0(
        "Final reproducibility report is inconsistent: ",
        n_controls_pass,
        " / ",
        n_controls,
        " controls passed."
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 6. Final timing
  # ----------------------------------------------------------
  
  workflow_end <- Sys.time()
  
  
  total_elapsed_seconds <- as.numeric(
    difftime(
      workflow_end,
      workflow_start,
      units =
        "secs"
    )
  )
  
  
  # ----------------------------------------------------------
  # 7. Final report
  # ----------------------------------------------------------
  
  cat(
    "\n\n============================================\n"
  )
  
  cat(
    "UREH WORKFLOW COMPLETED SUCCESSFULLY\n"
  )
  
  cat(
    "============================================\n\n"
  )
  
  
  print(
    script_log,
    row.names =
      FALSE
  )
  
  
  cat(
    "\nAnalytical scripts passed: ",
    sum(
      script_log$status ==
        "PASS"
    ),
    " / ",
    nrow(
      script_log
    ),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Reproducibility controls passed: ",
    n_controls_pass,
    " / ",
    n_controls,
    "\n",
    sep = ""
  )
  
  
  cat(
    "\nStarted: ",
    format(
      workflow_start,
      "%Y-%m-%d %H:%M:%S"
    ),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Finished: ",
    format(
      workflow_end,
      "%Y-%m-%d %H:%M:%S"
    ),
    "\n",
    sep = ""
  )
  
  
  cat(
    "Total elapsed time: ",
    round(
      total_elapsed_seconds,
      2
    ),
    " s\n",
    sep = ""
  )
  
  
  cat(
    "\nMain generated products:\n"
  )
  
  cat(
    " - results/tables/\n"
  )
  
  cat(
    " - results/spatial/ureh_results.gpkg\n"
  )
  
  cat(
    " - checks/reproducibility_report.csv\n"
  )
  
  
  cat(
    paste0(
      "\nPublication figures are intentionally outside ",
      "this automated workflow.\n"
    )
  )
  
  
  cat(
    paste0(
      "The reproducibility controls verify computational ",
      "consistency of the specified outputs; they do not ",
      "constitute external model validation.\n"
    )
  )
  
  
  cat(
    "============================================\n\n"
  )
  
  
  invisible(
    script_log
  )
}


# ============================================================
# Execute workflow
# ============================================================

run_ureh_workflow()