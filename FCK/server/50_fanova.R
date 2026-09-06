# ==========================================================================
# server/50_fanova.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 4425-6227  (group UIs, fANOVA (between + repeated measures))
#   WaPaa1_3.R lines 6524-6976  (post-hoc pairwise outputs)
# ==========================================================================
  # Group variable UI
  output$group_variable_ui <- renderUI({
    # Only depend on selected_group_vars and group_variables, NOT group_labels
    # This prevents re-rendering when the user switches group variables
    req(values$selected_group_vars)

    if(!is.null(values$group_variables)) {
      # Build UI elements
      ui_elements <- tagList(
        h5("Group Variable Status")
      )

      # If multiple group variables are available, show selector
      if(length(values$selected_group_vars) > 1) {
        # Preserve current selection if it exists and is valid
        current_selection <- isolate(input$fanova_group_var)
        if(is.null(current_selection) || !(current_selection %in% values$selected_group_vars)) {
          current_selection <- values$selected_group_vars[1]
        }

        ui_elements <- tagList(
          ui_elements,
          selectInput("fanova_group_var", "Select Group Variable for fANOVA:",
                      choices = values$selected_group_vars,
                      selected = current_selection),
          hr()
        )
      }

      ui_elements
    } else if(!is.null(values$group_labels)) {
      # Single group variable case (backwards compatibility)
      tagList(
        h5("Group Variable Status"),
        p(paste("Groups loaded:", length(unique(values$group_labels)), "groups")),
        p(paste("Total subjects:", length(values$group_labels)))
      )
    } else {
      tagList(
        p("No group variable detected."),
        p("To use Functional ANOVA, ensure you select a Group Variable during data import."),
        hr(),
        h5("Create Groups Manually"),
        numericInput("manual_n_groups", "Number of groups:", value = 2, min = 2, max = 10),
        actionButton("create_groups", "Create Groups", class = "btn-warning")
      )
    }
  })

  # Separate UI output for group info display (updates when fanova_group_var changes)
  output$fanova_group_info <- renderUI({
    # Get current group variable
    current_var <- input$fanova_group_var
    if(is.null(current_var) && !is.null(values$selected_group_vars)) {
      current_var <- values$selected_group_vars[1]
    }

    current_groups <- if(!is.null(current_var) && !is.null(values$group_variables) &&
                         current_var %in% colnames(values$group_variables)) {
      values$group_variables[[current_var]]
    } else {
      values$group_labels
    }

    if(is.null(current_groups)) return(NULL)

    n_groups <- length(unique(current_groups))

    ui_elements <- tagList(
      p(paste("Current group variable:", current_var)),
      p(paste("Number of groups:", n_groups)),
      p(paste("Total subjects:", length(current_groups)))
    )

    # If more than 2 groups, show group selector
    if(n_groups > 2) {
      ui_elements <- tagList(
        ui_elements,
        hr(),
        h5("Select Groups to Include in Analysis"),
        pickerInput(
          inputId = "fanova_groups_to_include",
          label = "Groups to include:",
          choices = levels(as.factor(current_groups)),
          selected = levels(as.factor(current_groups)),
          options = list(
            `actions-box` = TRUE,
            `selected-text-format` = "count > 3"
          ),
          multiple = TRUE
        ),
        helpText("Select at least 2 groups to compare. Deselect groups you want to exclude.")
      )
    }

    ui_elements
  })
  
  # Create groups manually
  observeEvent(input$create_groups, {
    req(values$data)
    n_subjects <- nrow(values$data)
    n_groups <- input$manual_n_groups

    group_labels <- rep(paste0("Group", 1:n_groups), length.out = n_subjects)
    values$group_labels <- factor(group_labels)

    showNotification("Groups created successfully!", type = "message", duration = 3)
  })

  # Helper function to get current fANOVA group labels based on selection
  get_fanova_group_labels <- reactive({
    # If user selected a specific group variable for fANOVA, use that
    if(!is.null(input$fanova_group_var) && !is.null(values$group_variables) &&
       input$fanova_group_var %in% colnames(values$group_variables)) {
      return(values$group_variables[[input$fanova_group_var]])
    }
    # Otherwise fall back to primary group_labels
    return(values$group_labels)
  })
  
  # UI for subject ID selection (for RM-ANOVA)
  output$subject_id_ui <- renderUI({
    req(values$data)
    
    # Get column names from original uploaded data if available
    col_names <- if(!is.null(values$uploaded_data)) {
      colnames(values$uploaded_data)
    } else {
      NULL
    }
    
    if(!is.null(col_names)) {
      # Filter out time-related columns (those matching the data columns)
      time_cols <- colnames(values$data)
      non_time_cols <- setdiff(col_names, time_cols)
      
      if(length(non_time_cols) > 0) {
        selectInput("rm_subject_id_var", 
                    "Subject ID variable:",
                    choices = c("", non_time_cols),
                    selected = "")
      } else {
        tagList(
          p("No ID variables found.", style = "color: orange;"),
          p("For repeated measures ANOVA, your data should include a subject ID column."),
          helpText("You can create subject IDs manually by numbering rows.")
        )
      }
    } else {
      tagList(
        p("Upload data to select subject ID variable.", style = "color: orange;")
      )
    }
  })
  
  # UI for repeated measures factor selection
  output$rm_factor_ui <- renderUI({
    req(values$data)

    # Get column names from original uploaded data if available
    col_names <- if(!is.null(values$uploaded_data)) {
      colnames(values$uploaded_data)
    } else {
      NULL
    }

    if(!is.null(col_names)) {
      # Filter out time-related columns
      time_cols <- colnames(values$data)
      non_time_cols <- setdiff(col_names, time_cols)

      if(length(non_time_cols) > 0) {
        selectInput("rm_factor_var",
                    "Repeated measures factor (visit/condition/time):",
                    choices = c("", non_time_cols),
                    selected = "")
      } else {
        tagList(
          p("No factor variables found.", style = "color: orange;"),
          p("For repeated measures ANOVA, you need a column indicating the repeated condition/visit/time."),
          helpText("Example: A 'visit' column with values like 'baseline', 'week1', 'week2', etc.")
        )
      }
    } else {
      tagList(
        p("Upload data to select repeated measures factor.", style = "color: orange;")
      )
    }
  })

  # UI for selecting which levels of the RM factor to include
  output$rm_factor_levels_ui <- renderUI({
    req(input$rm_factor_var)
    req(values$uploaded_data)

    if(input$rm_factor_var == "" || !(input$rm_factor_var %in% colnames(values$uploaded_data))) {
      return(NULL)
    }

    # Get the levels of the selected factor
    rm_factor_data <- fck_rm_column(values, input$rm_factor_var)
    factor_levels <- unique(as.character(rm_factor_data))
    factor_levels <- factor_levels[!is.na(factor_levels)]

    if(length(factor_levels) > 2) {
      tagList(
        hr(),
        h5("Select Conditions/Visits to Include"),
        pickerInput(
          inputId = "rm_levels_to_include",
          label = "Levels to include in analysis:",
          choices = factor_levels,
          selected = factor_levels,
          options = list(
            `actions-box` = TRUE,
            `selected-text-format` = "count > 3"
          ),
          multiple = TRUE
        ),
        helpText("Select at least 2 levels to compare. Deselect levels you want to exclude.")
      )
    } else {
      # Only 2 or fewer levels - no need for selection
      tagList(
        hr(),
        p(paste("Levels detected:", paste(factor_levels, collapse = ", "))),
        p(paste("Number of levels:", length(factor_levels)))
      )
    }
  })
  
  # Pointwise repeated-measures functional ANOVA with within-subject,
  # curve-wise permutation. (Named for what it does: it does NOT call the
  # rmfanova package -- see P1.b below.)
  # AUDIT (P3.5): the exported script writes this function out with deparse(),
  # so anything Shiny-only in its body becomes an undefined symbol in a plain R
  # session. It referenced input$rm_global_test, showNotification(),
  # withProgress() and incProgress(), which means the exported script for a
  # repeated-measures design could not run AT ALL -- it failed at
  # "object 'input' not found" the moment the permutation loop started. The
  # Shiny pieces are now arguments with inert defaults: the app passes the real
  # ones, the exported script gets a function that needs nothing but base R,
  # fda and (optionally) rmfanova.
  perform_rm_fanova <- function(fd_obj, subject_id, rm_factor, n_permutations = 200,
                                alpha = 0.05, run_global_test = FALSE,
                                progress = NULL, notify = NULL) {
    if (is.null(progress)) progress <- function(frac, detail = NULL) invisible(NULL)
    if (is.null(notify))   notify   <- function(msg, type = "message") {
      message(msg); invisible(NULL)
    }
    
    # AUDIT (P1.b): this used to refuse to run without the rmfanova package,
    # which the code then never successfully called. The procedure below has no
    # such dependency.
    
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    # Evaluate functional data at time points
    curves <- eval.fd(time_points, fd_obj)  # n_time x n_curves
    
    # Prepare the design
    subject_id <- as.factor(subject_id)
    rm_factor <- as.factor(rm_factor)
    
    # Get unique levels
    visits <- levels(rm_factor)
    n_visits <- length(visits)
    
    cat("Running repeated measures functional ANOVA...\n")
    cat("Number of subjects:", length(unique(subject_id)), "\n")
    cat("Number of visits/conditions:", n_visits, "\n")
    cat("Visits:", paste(visits, collapse = ", "), "\n")
    
    # AUDIT (P1.b): this used to try rmfanova::rmfanova(curves, id=, visit=),
    # then the same call positionally, then exists() probes for rm.fanova and
    # rmfANOVA, before falling through to the implementation below. None of
    # those signatures is the package's: rmfanova() takes a LIST of
    # condition-specific n x p matrices. So the package branch could never
    # succeed, and the app ran its own procedure while the UI said otherwise.
    #
    # Guessing at an API in a tryCatch cascade is not a fallback, it is a way of
    # not knowing which estimator produced your numbers. The speculative calls
    # are removed and the procedure is named for what it is. Adapting the data
    # to the real rmfanova API is a worthwhile separate change; it is not this
    # one, and pretending it already happened was the defect.
    # AUDIT (P3.3): the repeated-measures sums of squares were computed at every
    # time point, used for F, and then discarded -- and the effect size reported
    # downstream was rebuilt from a DIFFERENT, between-subjects decomposition
    # (SSB / SST over all curves, ignoring the subject margin). That is
    # classical eta-squared for an independent-groups design, reported under a
    # repeated-measures F. These vectors keep the RM sums of squares so partial
    # eta-squared can be formed from the same decomposition the test used:
    # SS_condition / (SS_condition + SS_error). Declared out here so the effect
    # size block below can tell "not computed" from "computed and zero".
    SS_visit_t    <- rep(NA_real_, n_time)
    SS_residual_t <- rep(NA_real_, n_time)
    n_complete_t  <- rep(NA_real_, n_time)

    tryCatch({
      rm_results <- NULL
      if(is.null(rm_results)) {
        cat("Performing manual repeated measures functional ANOVA...\n")
        
        # Manual RM-ANOVA implementation
        # This is a simplified version that accounts for repeated measures structure

        # Calculate mean curves for each visit
        visit_means <- matrix(NA, n_time, n_visits)
        visit_sds <- matrix(NA, n_time, n_visits)
        visit_sizes <- numeric(n_visits)

        for(i in 1:n_visits) {
          visit_idx <- which(rm_factor == visits[i])
          visit_sizes[i] <- length(unique(subject_id[visit_idx]))
          visit_curves <- curves[, visit_idx, drop = FALSE]
          visit_means[, i] <- rowMeans(visit_curves, na.rm = TRUE)

          if(length(visit_idx) > 1) {
            visit_sds[, i] <- apply(visit_curves, 1, sd, na.rm = TRUE)
          } else {
            visit_sds[, i] <- 0
          }
        }

        # === DIAGNOSTIC OUTPUT ===
        cat("\n=== VISIT MEANS DIAGNOSTIC ===\n")
        cat("Number of visits/conditions:", n_visits, "\n")
        cat("Visit names:", paste(visits, collapse=", "), "\n")
        cat("Visit means matrix dimensions:", dim(visit_means), "\n")
        cat("Visit sizes (n subjects per visit):", paste(visit_sizes, collapse=", "), "\n")

        for(i in 1:n_visits) {
          cat("\nVisit", i, "(", visits[i], "):\n")
          cat("  First 5 time points:", paste(round(visit_means[1:5, i], 3), collapse=", "), "\n")
          cat("  Mean across all time points:", round(mean(visit_means[, i]), 3), "\n")
          cat("  SD across all time points:", round(sd(visit_means[, i]), 3), "\n")
        }

        if(n_visits == 2) {
          cat("\nComparison between visits:\n")
          cat("  Are means identical?", all.equal(visit_means[,1], visit_means[,2]), "\n")
          cat("  Max absolute difference:", round(max(abs(visit_means[,1] - visit_means[,2])), 6), "\n")
        }
        cat("==============================\n\n")

        # Calculate F-statistics accounting for within-subject correlation
        # Use subject-specific deviations
        unique_subjects <- unique(subject_id)
        n_subjects <- length(unique_subjects)
        
        # First pass: Calculate observed F-statistics and build data matrices
        F_stat <- numeric(n_time)
        Y_matrices <- vector("list", n_time)  # Store for permutation testing
        Y_rows     <- vector("list", n_time)  # subject indices behind each Y_complete

        cat("=== Y_MATRIX CONSTRUCTION DIAGNOSTIC ===\n")

        for(t in 1:n_time) {
          # Create data matrix: subjects × visits
          Y_matrix <- matrix(NA, n_subjects, n_visits)
          
          for(i in 1:n_subjects) {
            subj <- unique_subjects[i]
            for(j in 1:n_visits) {
              idx <- which(subject_id == subj & rm_factor == visits[j])
              if(length(idx) > 0) {
                Y_matrix[i, j] <- curves[t, idx[1]]
              }
            }
          }
          
          # Remove subjects with missing data at this time point
          complete_rows <- complete.cases(Y_matrix)
          Y_complete <- Y_matrix[complete_rows, , drop = FALSE]

          # Diagnostic output for first time point only
          if(t == 1) {
            cat("\nTime point 1 - Y_matrix construction:\n")
            cat("  Y_matrix dimensions:", dim(Y_matrix), "(subjects × visits)\n")
            cat("  Y_complete dimensions:", dim(Y_complete), "\n")
            cat("  First 5 subjects, visit 1:", paste(round(Y_complete[1:min(5, nrow(Y_complete)), 1], 3), collapse=", "), "\n")
            cat("  First 5 subjects, visit 2:", paste(round(Y_complete[1:min(5, nrow(Y_complete)), 2], 3), collapse=", "), "\n")
            cat("  Column means:", paste(round(colMeans(Y_complete), 3), collapse=", "), "\n")
            cat("  Are columns identical?", all.equal(Y_complete[,1], Y_complete[,2]), "\n")
          }

          # Store for permutation testing. The SUBJECT INDICES are kept too:
          # the complete-case set differs between time points, and a curve-wise
          # permutation has to apply one subject's relabelling wherever that
          # subject appears (P1.a).
          Y_matrices[[t]] <- Y_complete
          Y_rows[[t]] <- which(complete_rows)

          if(nrow(Y_complete) >= 2 && n_visits >= 2) {
            # Perform repeated measures ANOVA at this time point
            # Remove subject mean (within-subject centering)
            # AUDIT (P0.4): this used to be
            #     Y_centered  <- Y_complete - rowMeans(Y_complete)
            #     SS_residual <- sum(Y_centered^2)
            # Removing the SUBJECT means leaves the VISIT effect in the
            # residual, so SS_residual was SS_error + SS_visit exactly. The
            # denominator was inflated and F came out at 22.8% of its correct
            # value on a 12x4 fixture (8.49 where stats::aov gives 37.26), so
            # the test was badly conservative and would miss real effects.
            # The residual of a one-factor repeated-measures design removes
            # BOTH margins and adds the grand mean back.
            grand_mean <- mean(Y_complete, na.rm = TRUE)
            visit_means_t <- colMeans(Y_complete, na.rm = TRUE)
            subject_means <- rowMeans(Y_complete, na.rm = TRUE)

            E <- sweep(sweep(Y_complete, 1, subject_means, "-"),
                       2, visit_means_t, "-") + grand_mean

            SS_visit <- sum(nrow(Y_complete) * (visit_means_t - grand_mean)^2)
            SS_residual <- sum(E^2, na.rm = TRUE)

            SS_visit_t[t]    <- SS_visit
            SS_residual_t[t] <- SS_residual
            n_complete_t[t]  <- nrow(Y_complete)

            df_visit <- n_visits - 1
            df_residual <- (nrow(Y_complete) - 1) * (n_visits - 1)
            
            if(df_residual > 0 && SS_residual > 0) {
              F_stat[t] <- (SS_visit / df_visit) / (SS_residual / df_residual)
            } else {
              F_stat[t] <- NA
            }

            # Diagnostic output for first time point
            if(t == 1) {
              cat("\nTime point 1 - F-statistic calculation:\n")
              cat("  Grand mean:", round(grand_mean, 3), "\n")
              cat("  Visit means:", paste(round(visit_means_t, 3), collapse=", "), "\n")
              cat("  SS_visit:", round(SS_visit, 3), "\n")
              cat("  SS_residual:", round(SS_residual, 3), "\n")
              cat("  df_visit:", df_visit, ", df_residual:", df_residual, "\n")
              cat("  F-statistic:", round(F_stat[t], 3), "\n")
            }
          } else {
            F_stat[t] <- NA
          }
        }
        cat("=========================================\n\n")
        
        # Handle NAs in F-statistics
        # AUDIT (P5.7): this said F_stat[is.na(F_stat)] <- 0, and the
        # permutation branch and the p-values did the same thing in their own
        # way (F* set to 0, p set to 1). Undefined is not the same as "no
        # effect". A time point where the repeated-measures residual sum of
        # squares is exactly zero has NO F ratio -- reporting F = 0, p = 1
        # asserts that the conditions were tested there and found equal, which
        # is a claim the data cannot support. P4.7 fixed exactly this in the
        # between-subjects branch; the RM branch kept the old behaviour.
        # NA is now carried through to the p-value and out to the readout.
        n_undefined_F <- sum(!is.finite(F_stat))
        
        # PERMUTATION TEST for p-values (matching pairwise approach)
        cat("Computing permutation-based p-values (", n_permutations, " permutations)...\n")
        
        F_stat_perm <- matrix(NA, n_time, n_permutations)
        
        # P3.5: this loop used to sit inside withProgress({...}), which
        # evaluates its expression in the CALLER's frame, so the assignments
        # below reached F_stat_perm. local() would not -- it makes a new
        # environment and the writes would be discarded. The progress sink is
        # a plain callback now, so no wrapper is needed at all.
          for(perm in 1:n_permutations) {
            if(perm %% 20 == 0) progress(20 / n_permutations, "Permutation testing")

            # AUDIT (P1.a): the condition relabelling used to be drawn INSIDE
            # the time loop, so one subject could be read as 2-1-3 at t1 and
            # 3-2-1 at t2. The exchangeable unit in a functional design is the
            # whole trajectory, not the value at a time point. Each pointwise
            # p-value was still marginally valid (permuting conditions within a
            # subject is the right null at a single t), so no p-value the app
            # has reported is wrong because of this -- but the permuted curves
            # were temporally scrambled, which makes the joint null across time
            # meaningless and blocks any global max-F or L2 statistic.
            #
            # One relabelling per subject per replicate, applied at every time
            # point where that subject appears.
            perm_map <- matrix(NA_integer_, n_subjects, n_visits)
            for(i in seq_len(n_subjects)) perm_map[i, ] <- sample.int(n_visits)

            for(t in 1:n_time) {
              Y_complete <- Y_matrices[[t]]
              
              if(!is.null(Y_complete) && nrow(Y_complete) >= 2 && n_visits >= 2) {
                rows_t <- Y_rows[[t]]
                Y_perm <- Y_complete
                for(i in seq_len(nrow(Y_complete))) {
                  Y_perm[i, ] <- Y_complete[i, perm_map[rows_t[i], ]]
                }
                
                # Calculate permuted F-statistic
                subject_means_perm <- rowMeans(Y_perm, na.rm = TRUE)
                grand_mean_perm <- mean(Y_perm, na.rm = TRUE)
                visit_means_perm <- colMeans(Y_perm, na.rm = TRUE)
                
                E_perm <- sweep(sweep(Y_perm, 1, subject_means_perm, "-"),
                                2, visit_means_perm, "-") + grand_mean_perm
                SS_visit_perm <- sum(nrow(Y_perm) * (visit_means_perm - grand_mean_perm)^2)
                SS_residual_perm <- sum(E_perm^2, na.rm = TRUE)
                
                df_visit <- n_visits - 1
                df_residual <- (nrow(Y_perm) - 1) * (n_visits - 1)
                
                if(df_residual > 0 && SS_residual_perm > 0) {
                  F_stat_perm[t, perm] <- (SS_visit_perm / df_visit) / (SS_residual_perm / df_residual)
                } else {
                  # P5.7: a permutation draw with no residual variation supplies
                  # no null value. It is dropped from the reference set (which
                  # the p-value below counts), not entered as a zero -- a zero
                  # is the SMALLEST possible F, so it never exceeds the
                  # observed statistic and silently biases every p-value down.
                  F_stat_perm[t, perm] <- NA_real_
                }
              } else {
                F_stat_perm[t, perm] <- NA_real_
              }
            }
          }
        
        # Calculate permutation-based p-values
        p_values_pointwise <- numeric(n_time)
        for(t in 1:n_time) {
          # AUDIT (P1.1): mean(perm >= obs) can be exactly 0, reporting p = 0
          # for a Monte Carlo test of finitely many draws. The observed
          # statistic counts as one of its own null draws.
          # P5.7: an undefined OBSERVED statistic has no test either.
          .np <- sum(is.finite(F_stat_perm[t, ]))
          p_values_pointwise[t] <- if (!is.finite(F_stat[t]) || .np < 1) NA_real_ else
            (1 + sum(F_stat_perm[t, ] >= F_stat[t], na.rm = TRUE)) / (1 + .np)
        }
        
        # P5.7: NAs are NOT converted to p = 1. A time point with no defined
        # F, or with no usable null draws, has no test, and says so.
        
        cat("Permutation testing complete.\n")
        
      } else {
        # Assemble the results
        # The structure depends on the package version
        
        # Try to extract pointwise statistics
        if(!is.null(rm_results$pointwise)) {
          F_stat <- rm_results$pointwise$stat
          p_values_pointwise <- rm_results$pointwise$pval
        } else if(!is.null(rm_results$stat)) {
          F_stat <- rm_results$stat
          p_values_pointwise <- rm_results$pval
        } else {
          # Calculate manually as fallback
          visit_means <- matrix(NA, n_time, n_visits)
          for(i in 1:n_visits) {
            visit_idx <- which(rm_factor == visits[i])
            visit_means[, i] <- rowMeans(curves[, visit_idx, drop = FALSE])
          }
          
          # Simple F-statistic calculation
          overall_mean <- rowMeans(curves)
          SSB <- rowSums((visit_means - overall_mean)^2) * length(unique(subject_id))
          SST <- rowSums((curves - overall_mean)^2)
          
          F_stat <- SSB / (SST - SSB + 1e-10)
          p_values_pointwise <- pf(F_stat, n_visits - 1, ncol(curves) - n_visits, lower.tail = FALSE)
        }
        
        # Calculate visit means for plotting
        visit_means <- matrix(NA, n_time, n_visits)
        visit_sds <- matrix(NA, n_time, n_visits)
        visit_sizes <- numeric(n_visits)

        for(i in 1:n_visits) {
          visit_idx <- which(rm_factor == visits[i])
          visit_sizes[i] <- length(unique(subject_id[visit_idx]))
          visit_curves <- curves[, visit_idx, drop = FALSE]
          visit_means[, i] <- rowMeans(visit_curves, na.rm = TRUE)

          if(length(visit_idx) > 1) {
            visit_sds[, i] <- apply(visit_curves, 1, sd, na.rm = TRUE)
          } else {
            visit_sds[, i] <- 0
          }
        }
      }
      
      # Continue with common processing...
      # Adjust p-values for multiple comparisons
      p_values_adjusted <- p.adjust(p_values_pointwise, method = "fdr")
      # P5.7: NA is not significant and is not non-significant. It is excluded
      # and counted, the same rule the between-subjects branch uses (P4.7).
      sig_regions <- !is.na(p_values_adjusted) & p_values_adjusted < alpha
      n_undefined <- sum(is.na(p_values_pointwise))
      if (n_undefined > 0)
        notify(sprintf(
          paste("%d of %d time points have no defined repeated-measures F (no residual",
                "variation there). They are reported as NA, not as p = 1."),
          n_undefined, length(p_values_pointwise)), "warning")
      
      # ---- Effect size --------------------------------------------------
      # AUDIT (P3.3): this block used to compute
      #     SSB[t] = sum_i n_i (visit_mean - overall_mean)^2
      #     SST[t] = sum_j (curve_j - overall_mean)^2
      #     eta_squared = SSB / SST
      # over ALL curves, with no subject margin anywhere in it. That is
      # classical eta-squared for a BETWEEN-subjects one-way design. In a
      # repeated-measures design the subject differences sit inside that SST,
      # so the denominator carries between-subject variance the F test has
      # already partialled out, and the reported effect size is biased DOWN --
      # arbitrarily far down, since it shrinks as between-subject spread grows.
      # It was also computed on a different case set from the test: the F used
      # complete cases per time point, this used every curve.
      #
      # The correct companion to a repeated-measures F is PARTIAL eta-squared,
      # from the decomposition the test itself used:
      #     partial eta^2 = SS_condition / (SS_condition + SS_error)
      # which is exactly df_visit * F / (df_visit * F + df_residual). Those sums
      # of squares are recorded per time point in the loop above.
      eta_squared <- rep(NA_real_, n_time)
      ok_ss <- is.finite(SS_visit_t) & is.finite(SS_residual_t) &
               (SS_visit_t + SS_residual_t) > 0
      eta_squared[ok_ss] <- SS_visit_t[ok_ss] /
                            (SS_visit_t[ok_ss] + SS_residual_t[ok_ss])
      # A time point where both sums of squares are zero has no variation at
      # all to explain; 0 is the honest value there, NA where the RM
      # decomposition could not be formed (fewer than two complete subjects).
      eta_squared[is.finite(SS_visit_t) & is.finite(SS_residual_t) &
                  (SS_visit_t + SS_residual_t) == 0] <- 0
      eta_squared_type <- "partial"

      # Kept for reporting: the classical (between-subjects) eta-squared this
      # tab used to show, so a user comparing against an earlier run can see
      # both and see why they differ.
      overall_mean <- rowMeans(curves)
      SSB_classical <- numeric(n_time)
      SST_classical <- numeric(n_time)
      for(t in 1:n_time) {
        for(i in 1:n_visits) {
          SSB_classical[t] <- SSB_classical[t] +
            visit_sizes[i] * (visit_means[t, i] - overall_mean[t])^2
        }
        SST_classical[t] <- sum((curves[t, ] - overall_mean[t])^2, na.rm = TRUE)
      }
      eta_squared_classical <- ifelse(SST_classical > 0,
                                      SSB_classical / SST_classical, 0)
      
      # Global test statistic (L2 norm)
      # P1.3: integrate over time rather than summing over grid points
      # P3.3 (continued): this used to read the same between-subjects SSB the
      # old effect size was built from -- a global summary of the RM tab built
      # out of a decomposition the RM test does not use. It is now the
      # integrated CONDITION sum of squares from the repeated-measures
      # decomposition, per complete subject. THIS NUMBER CHANGES relative to
      # earlier versions of the app. It is descriptive only: p_value_L2 is NA
      # here because this procedure computes no global p-value for a
      # repeated-measures design (rmfanova does -- see run_global_test).
      L2_stat <- fck_l2_norm(
        ifelse(is.finite(SS_visit_t) & is.finite(n_complete_t) & n_complete_t > 0,
               SS_visit_t / n_complete_t, 0),
        time_points)
      L2_stat_classical <- fck_l2_norm(SSB_classical / length(subject_id), time_points)
      p_value_L2 <- NA  # this app's own procedure computes no global p-value

      # AUDIT (P2.6): the pointwise procedure above answers "where do the
      # conditions differ", never "do they differ at all". rmfanova (Kurylo &
      # Smaga 2023) supplies the global test, on a complete balanced design.
      # It is optional, off unless asked for, and reports which subjects it had
      # to drop. See server/09b_helpers_rmfanova.R for why only three of its
      # fifteen outputs are shown by default.
      rm_global <- NULL
      if (isTRUE(run_global_test)) {
        rm_global <- tryCatch(
          fck_rmfanova_global(curves, subject_id, rm_factor,
                              n_perm = min(n_permutations, 2000),
                              n_boot = min(n_permutations, 2000)),
          error = function(e) list(ok = FALSE, reason = conditionMessage(e)))
        if (isTRUE(rm_global$ok)) {
          notify(sprintf("Global RM test on %d of %d subjects with a complete design.",
                         rm_global$n_complete, rm_global$n_total), "message")
        } else {
          notify(paste("Global RM test not run:", rm_global$reason), "warning")
        }
      }
      
      # Pointwise bootstrap percentile intervals (NOT simultaneous bands)
      # AUDIT (P1.2): 100 replicates gives a 2.5% quantile estimated from the
      # 2nd-3rd order statistic. These are POINTWISE percentile intervals, not
      # simultaneous functional bands -- see the label below.
      n_boot <- 2000
      visit_means_boot <- array(NA, dim = c(n_time, n_visits, n_boot))
      
      unique_subjects <- unique(subject_id)
      
      for(boot in 1:n_boot) {
        # Resample subjects (not observations)
        boot_subjects <- sample(unique_subjects, replace = TRUE)
        
        for(i in 1:n_visits) {
          boot_curves_list <- list()
          for(subj in boot_subjects) {
            # Get curves for this subject in this visit
            subj_visit_idx <- which(subject_id == subj & rm_factor == visits[i])
            if(length(subj_visit_idx) > 0) {
              boot_curves_list[[length(boot_curves_list) + 1]] <- curves[, subj_visit_idx[1], drop = FALSE]
            }
          }
          
          if(length(boot_curves_list) > 0) {
            boot_curves_matrix <- do.call(cbind, boot_curves_list)
            visit_means_boot[, i, boot] <- rowMeans(boot_curves_matrix)
          } else {
            visit_means_boot[, i, boot] <- visit_means[, i]
          }
        }
      }
      
      # Pointwise 95% percentile intervals
      visit_means_lower <- matrix(NA, n_time, n_visits)
      visit_means_upper <- matrix(NA, n_time, n_visits)
      
      for(i in 1:n_visits) {
        for(t in 1:n_time) {
          quantiles <- quantile(visit_means_boot[t, i, ], probs = c(0.025, 0.975), na.rm = TRUE)
          visit_means_lower[t, i] <- quantiles[1]
          visit_means_upper[t, i] <- quantiles[2]
        }
      }
      
      # DIAGNOSTIC: Verify F_stat before return
      cat("\n=== FINAL F-STAT DIAGNOSTIC ===\n")
      cat("F_stat vector length:", length(F_stat), "\n")
      cat("First 5 F-statistics:", paste(round(F_stat[1:min(5, length(F_stat))], 3), collapse=", "), "\n")
      cat("Max F-statistic:", round(max(F_stat, na.rm = TRUE), 3), "\n")
      cat("Mean F-statistic:", round(mean(F_stat, na.rm = TRUE), 3), "\n")
      cat("Number of significant time points (raw p < 0.05):", sum(p_values_pointwise < 0.05, na.rm = TRUE), "/", length(p_values_pointwise), "\n")
      cat("Number of significant time points (FDR adjusted):", sum(p_values_adjusted < 0.05, na.rm = TRUE), "/", length(p_values_adjusted), "\n")
      cat("Min p-value:", format.pval(min(p_values_pointwise, na.rm = TRUE), digits = 4), "\n")
      cat("================================\n\n")

      # Return results in same format as between-subjects ANOVA
      return(list(
        design = "within",
        time_points = time_points,
        group_means = visit_means,
        group_sds = visit_sds,
        group_means_lower = visit_means_lower,
        group_means_upper = visit_means_upper,
        group_labels = rm_factor,
        groups = visits,
        n_groups = n_visits,
        group_sizes = visit_sizes,
        F_stat = F_stat,
        p_values_pointwise = p_values_pointwise,
        p_values_adjusted = p_values_adjusted,
        sig_regions = sig_regions,
        L2_stat = L2_stat,
        p_value_L2 = p_value_L2,
        rm_global = rm_global,     # P2.6: optional rmfanova global test, or NULL
        eta_squared = eta_squared,
        df_between = n_visits - 1,
        # AUDIT (P0.4): was length(unique(subject_id)) * (n_visits - 1), which
        # counts n*(k-1) instead of the (n-1)*(k-1) the F test actually uses.
        df_within = (length(unique(subject_id)) - 1) * (n_visits - 1),
        eta_squared_type = eta_squared_type,
        n_undefined = n_undefined,
        L2_stat_classical = L2_stat_classical,
        eta_squared_classical = eta_squared_classical,
        SS_condition = SS_visit_t,
        SS_error = SS_residual_t,
        n_complete_per_time = n_complete_t,
        alpha = alpha,
        n_permutations = n_permutations,
        subject_id = subject_id,
        rm_factor = rm_factor
      ))
      
    }, error = function(e) {
      stop(paste("Error in RM-ANOVA:", e$message))
    })
  }
  
  # Functional ANOVA function - ENHANCED
  perform_functional_anova <- function(fd_obj, group_labels, n_permutations = 200, 
                                       test_type = "both", alpha = 0.05) {
    
    n_curves <- ncol(fd_obj$coefs)
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    curves <- eval.fd(time_points, fd_obj)
    
    group_labels <- as.factor(group_labels)
    groups <- levels(group_labels)
    n_groups <- length(groups)
    
    group_means <- matrix(NA, n_time, n_groups)
    group_sds <- matrix(NA, n_time, n_groups)
    group_sizes <- numeric(n_groups)
    
    for(i in 1:n_groups) {
      group_idx <- which(group_labels == groups[i])
      group_sizes[i] <- length(group_idx)
      group_curves <- curves[, group_idx, drop = FALSE]
      group_means[, i] <- rowMeans(group_curves)
      
      # Calculate standard deviation for each time point
      if(length(group_idx) > 1) {
        group_sds[, i] <- apply(group_curves, 1, sd)
      } else {
        group_sds[, i] <- 0
      }
    }
    
    overall_mean <- rowMeans(curves)
    
    # Calculate F-statistics
    SSB <- numeric(n_time)
    SSW <- numeric(n_time)
    
    for(t in 1:n_time) {
      for(i in 1:n_groups) {
        SSB[t] <- SSB[t] + group_sizes[i] * (group_means[t, i] - overall_mean[t])^2
      }
      
      for(i in 1:n_groups) {
        group_idx <- which(group_labels == groups[i])
        for(j in group_idx) {
          SSW[t] <- SSW[t] + (curves[t, j] - group_means[t, i])^2
        }
      }
    }
    
    df_between <- n_groups - 1
    df_within <- n_curves - n_groups
    
    # AUDIT (P4.7): (SSB/df) / (SSW/df) with no guard. Constant curves at a
    # time point give SSW = 0, so F is Inf when SSB > 0 and NaN when both are
    # zero -- and NaN then propagates into the permutation comparison, where
    # `NaN >= NaN` is NA and the p-value silently comes from fewer draws than it
    # claims. The repeated-measures branch already refuses this case
    # (df_residual > 0 && SS_residual > 0); the between-subjects branch now does
    # too, and NA is carried through to the p-value rather than smuggled in as
    # a number. Groups that genuinely differ with zero within-group spread are
    # a degenerate design, not evidence of infinite strength.
    ss_floor <- .Machine$double.eps * max(1, max(abs(curves), na.rm = TRUE))^2 * n_curves
    F_stat <- ifelse(SSW > ss_floor, (SSB / df_between) / (SSW / df_within), NA_real_)
    if (any(!is.finite(F_stat)))
      warning(sprintf(
        "F is undefined at %d of %d time points (no within-group variation there); those points are reported as NA.",
        sum(!is.finite(F_stat)), length(F_stat)))
    # P1.3: a functional L2 norm, not a grid-density-dependent vector norm
    L2_stat <- fck_l2_norm(SSB / n_curves, time_points)
    
    # Permutation test for p-values
    F_stat_perm <- matrix(NA, n_time, n_permutations)
    L2_stat_perm <- numeric(n_permutations)
    
    # AUDIT (P2.5): this was a four-deep interpreted loop --
    #     for(perm) for(t) for(group) for(curve in group)
    # which is O(n_perm x n_time x n_curves). At the P1.1 default of 5,000
    # permutations, 100 evaluation points and the Circaflex sample that is ~650
    # million R-level iterations, i.e. not runnable. The same two sums are
    # matrix operations: SSB is a weighted row-sum over group means, and SSW is
    # the row-sums of squared deviations from each curve's own permuted group
    # mean. Identical arithmetic, verified against the loop in
    # tests/testthat/test-p2-corrections.R.
    for(perm in 1:n_permutations) {
      perm_labels <- sample(group_labels)
      gidx <- match(as.character(perm_labels), groups)   # group of each curve

      perm_means <- vapply(seq_len(n_groups), function(i)
        rowMeans(curves[, gidx == i, drop = FALSE]), numeric(n_time))
      n_per_group <- tabulate(gidx, nbins = n_groups)

      SSB_perm <- as.vector((perm_means - overall_mean)^2 %*% n_per_group)
      SSW_perm <- rowSums((curves - perm_means[, gidx, drop = FALSE])^2)

      F_stat_perm[, perm] <- ifelse(SSW_perm > ss_floor,
                                    (SSB_perm / df_between) / (SSW_perm / df_within),
                                    NA_real_)
      L2_stat_perm[perm] <- fck_l2_norm(SSB_perm / n_curves, time_points)
    }
    
    # Calculate p-values
    p_values_pointwise <- numeric(n_time)
    for(t in 1:n_time) {
      # AUDIT (P1.1): (1 + #{T* >= T}) / (1 + B) -- a Monte Carlo p is never 0.
      # P4.7: a time point whose observed F is undefined has no test, and must
      # not be handed a p-value of 1 as though it had one.
      .np <- sum(is.finite(F_stat_perm[t, ]))
      p_values_pointwise[t] <- if (!is.finite(F_stat[t]) || .np < 1) NA_real_ else
        (1 + sum(F_stat_perm[t, ] >= F_stat[t], na.rm = TRUE)) / (1 + .np)
    }
    
    p_value_L2 <- (1 + sum(L2_stat_perm >= L2_stat, na.rm = TRUE)) /
                      (1 + sum(is.finite(L2_stat_perm)))
    
    p_values_adjusted <- p.adjust(p_values_pointwise, method = "fdr")
    # P4.7: a time point with no defined test is not significant, and is not
    # non-significant either. It is excluded, and counted so the readout can
    # say how many were dropped rather than folding them into "not significant".
    sig_regions <- !is.na(p_values_adjusted) & p_values_adjusted < alpha
    n_undefined <- sum(is.na(p_values_pointwise))
    
    SST <- SSB + SSW
    # Classical eta-squared: correct for this INDEPENDENT-groups design, where
    # SST carries no subject margin to partial out. The repeated-measures
    # branch reports partial eta-squared instead (P3.3); both carry an
    # eta_squared_type label so the readout never has to guess which it holds.
    eta_squared <- ifelse(SST > 0, SSB / SST, 0)
    
    # Pointwise bootstrap percentile intervals (NOT simultaneous bands)
    # AUDIT (P1.2): pointwise percentile intervals, 100 -> 2000 replicates.
    n_boot <- 2000
    group_means_boot <- array(NA, dim = c(n_time, n_groups, n_boot))
    
    for(boot in 1:n_boot) {
      for(i in 1:n_groups) {
        group_idx <- which(group_labels == groups[i])
        if(length(group_idx) > 0) {
          boot_idx <- sample(group_idx, replace = TRUE)
          group_means_boot[, i, boot] <- rowMeans(curves[, boot_idx, drop = FALSE])
        }
      }
    }
    
    # Pointwise 95% percentile intervals
    group_means_lower <- matrix(NA, n_time, n_groups)
    group_means_upper <- matrix(NA, n_time, n_groups)
    
    for(i in 1:n_groups) {
      for(t in 1:n_time) {
        quantiles <- quantile(group_means_boot[t, i, ], probs = c(0.025, 0.975), na.rm = TRUE)
        group_means_lower[t, i] <- quantiles[1]
        group_means_upper[t, i] <- quantiles[2]
      }
    }
    
    return(list(
      time_points = time_points,
      group_means = group_means,
      group_sds = group_sds,  # Added standard deviations
      group_means_lower = group_means_lower,
      group_means_upper = group_means_upper,
      group_labels = group_labels,
      groups = groups,
      n_groups = n_groups,
      group_sizes = group_sizes,
      F_stat = F_stat,
      p_values_pointwise = p_values_pointwise,
      p_values_adjusted = p_values_adjusted,
      sig_regions = sig_regions,
      n_undefined = n_undefined,
      L2_stat = L2_stat,
      p_value_L2 = p_value_L2,
      eta_squared = eta_squared,
      eta_squared_type = "classical",
      df_between = df_between,
      df_within = df_within,
      alpha = alpha,
      n_permutations = n_permutations
    ))
  }
  
  # Run functional ANOVA - ENHANCED for both between and within designs
  observeEvent(input$run_fanova, {
    if(is.null(values$data)) {
      showNotification("Please load data first!", type = "error", duration = 5)
      return()
    }
    
    # Validate based on design type
    if(input$fanova_design == "between") {
      # Between-subjects validation - REQUIRES group labels from preprocessing
      fanova_groups <- get_fanova_group_labels()
      if(is.null(fanova_groups)) {
        showNotification("For between-subjects ANOVA: Please define group labels in the Data Preprocessing tab first!",
                         type = "error", duration = 5)
        return()
      }

      if(length(unique(fanova_groups)) < 2) {
        showNotification("Need at least 2 groups for ANOVA!", type = "error", duration = 5)
        return()
      }

      cat("Running BETWEEN-SUBJECTS functional ANOVA...\n")
      cat("Using group variable:", if(!is.null(input$fanova_group_var)) input$fanova_group_var else "primary", "\n")
      
    } else if(input$fanova_design == "within") {
      # Within-subjects (repeated measures) validation
      if(is.null(input$rm_subject_id_var) || input$rm_subject_id_var == "") {
        showNotification("Please select a Subject ID variable for repeated measures ANOVA!", 
                         type = "error", duration = 5)
        return()
      }
      
      if(is.null(input$rm_factor_var) || input$rm_factor_var == "") {
        showNotification("Please select a repeated measures factor (visit/condition/time)!", 
                         type = "error", duration = 5)
        return()
      }
      
      # Extract subject ID and RM factor from uploaded data
      if(is.null(values$uploaded_data)) {
        showNotification("Original data with ID and factor variables not available!", 
                         type = "error", duration = 5)
        return()
      }
      
      # Get the variables
      subject_id_col <- input$rm_subject_id_var
      rm_factor_col <- input$rm_factor_var
      
      if(!(subject_id_col %in% colnames(values$uploaded_data)) || 
         !(rm_factor_col %in% colnames(values$uploaded_data))) {
        showNotification("Selected variables not found in data!", 
                         type = "error", duration = 5)
        return()
      }
      
      subject_id_data <- values$uploaded_data[[subject_id_col]]
      rm_factor_data <- values$uploaded_data[[rm_factor_col]]
      
      # Check that we have multiple levels of the RM factor
      if(length(unique(rm_factor_data)) < 2) {
        showNotification("Need at least 2 levels of the repeated measures factor!", 
                         type = "error", duration = 5)
        return()
      }
      
      cat("Running WITHIN-SUBJECTS (Repeated Measures) functional ANOVA...\n")
      cat("Subject ID variable:", subject_id_col, "\n")
      cat("RM factor variable:", rm_factor_col, "\n")
      cat("Number of unique subjects:", length(unique(subject_id_data)), "\n")
      cat("Number of conditions/visits:", length(unique(rm_factor_data)), "\n")
      cat("Total number of observations:", length(subject_id_data), "\n")
      
      # Validate this is truly repeated measures
      n_subjects <- length(unique(subject_id_data))
      n_obs <- length(subject_id_data)
      if(n_subjects == n_obs) {
        cat("\n*** WARNING ***\n")
        cat("Number of subjects equals number of observations!\n")
        cat("This suggests you may not have true repeated measures data.\n")
        cat("Each subject should appear multiple times (once per condition/visit).\n")
        cat("Consider using between-subjects design instead.\n")
        cat("***************\n\n")
      } else {
        cat("Expected observations per subject:", round(n_obs / n_subjects, 2), "\n")
      }
    }
    
    tryCatch({
      # Determine which data to use based on user selection
      fd_to_use <- NULL
      
      # Check if user wants to use warped curves and if they're available
      if(input$fanova_data_source == "warped" && !is.null(values$warping_results)) {
        cat("Using time-warped curves for FANOVA\n")
        
        # Try to get the registered fd object from warping results
        if(!is.null(values$warping_results$regfd)) {
          fd_to_use <- values$warping_results$regfd
          showNotification("Using time-warped curves for FANOVA", type = "message", duration = 3)
        } else if(!is.null(values$warping_results$registered_curves)) {
          # If no regfd but registered curves exist, create fd from them
          n_time <- nrow(values$warping_results$registered_curves)
          time_points <- values$warping_results$time_points
          if(is.null(time_points)) {
            time_points <- seq(0, 1, length.out = n_time)
          }
          basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time-2))
          # Use lambda=0: registered curves are already processed, just need fd representation
          fd_to_use <- smooth.basis(time_points, values$warping_results$registered_curves, 
                                    fdPar(basis, 2, 0))$fd
          showNotification("Using time-warped curves for FANOVA", type = "message", duration = 3)
        } else {
          showNotification("No warped curves available, using original data", type = "warning", duration = 5)
        }
      }
      
      # If no warped data was used or available, use original data
      if(is.null(fd_to_use)) {
        cat("Using original curves for FANOVA\n")
        
        # CRITICAL: Check if data has already been smoothed in Data Preprocessing.
        # MERGED APP: when it has not, fck_ensure_fd_obj() builds an
        # INTERPOLATING representation and says so, instead of quietly
        # smoothing onto min(20, n_time - 2) basis functions as WaPaa did.
        # Nothing here re-smooths data that the preprocessing step already
        # smoothed: values$fd_obj is reused as-is.
        if(is.null(values$fd_obj)) {
          if(!fck_ensure_fd_obj(values)) return()
        }
        fd_to_use <- values$fd_obj
        showNotification("Using original curves for FANOVA", type = "message", duration = 3)
      }
      
      # Perform FANOVA based on design type
      if(input$fanova_design == "between") {
        # Between-subjects ANOVA

        # Get current group labels from the selected fANOVA group variable
        current_group_labels <- get_fanova_group_labels()

        # Filter by selected groups if user has made a selection
        groups_to_include <- input$fanova_groups_to_include
        if(!is.null(groups_to_include) && length(groups_to_include) >= 2) {
          # Find indices of subjects in selected groups
          include_idx <- which(current_group_labels %in% groups_to_include)

          if(length(include_idx) < 2) {
            showNotification("Not enough subjects in selected groups!", type = "error", duration = 5)
            return()
          }

          # Subset the fd object - need to extract curves, subset, and recreate
          n_time <- 100
          time_points_eval <- seq(0, 1, length.out = n_time)
          all_curves <- eval.fd(time_points_eval, fd_to_use)
          subset_curves <- all_curves[, include_idx, drop = FALSE]

          # Recreate fd object from subset
          basis <- fd_to_use$basis
          fd_to_use <- smooth.basis(time_points_eval, subset_curves, basis)$fd

          # Subset group labels and drop unused levels
          current_group_labels <- droplevels(current_group_labels[include_idx])

          cat("Filtered to groups:", paste(groups_to_include, collapse = ", "), "\n")
          cat("Subjects included:", length(include_idx), "\n")

          # Store which groups were selected for reference
          values$fanova_selected_groups <- groups_to_include
        } else if(!is.null(groups_to_include) && length(groups_to_include) < 2) {
          showNotification("Please select at least 2 groups for comparison!", type = "error", duration = 5)
          return()
        }

        values$fanova_results <- perform_functional_anova(
          fd_obj = fd_to_use,
          group_labels = current_group_labels,
          n_permutations = input$n_permutations,
          test_type = input$fanova_test_type,
          alpha = input$alpha_level
        )

        values$fanova_results$design <- "between"
        
      } else if(input$fanova_design == "within") {
        # Within-subjects (repeated measures) ANOVA
        # MERGED APP: aligned to the curves (see FCK/server/52_posthoc_source.R).
        subject_id_data <- fck_rm_column(values, input$rm_subject_id_var)
        rm_factor_data  <- fck_rm_column(values, input$rm_factor_var)

        # Filter by selected levels if user has made a selection
        levels_to_include <- input$rm_levels_to_include
        if(!is.null(levels_to_include) && length(levels_to_include) >= 2) {
          # Find indices of observations with selected levels
          include_idx <- which(as.character(rm_factor_data) %in% levels_to_include)

          if(length(include_idx) < 2) {
            showNotification("Not enough observations in selected levels!", type = "error", duration = 5)
            return()
          }

          # Subset the fd object
          n_time <- 100
          time_points_eval <- seq(0, 1, length.out = n_time)
          all_curves <- eval.fd(time_points_eval, fd_to_use)
          subset_curves <- all_curves[, include_idx, drop = FALSE]

          # Recreate fd object from subset
          basis <- fd_to_use$basis
          fd_to_use <- smooth.basis(time_points_eval, subset_curves, basis)$fd

          # Subset the subject ID and RM factor data
          subject_id_data <- subject_id_data[include_idx]
          rm_factor_data <- droplevels(as.factor(rm_factor_data[include_idx]))

          cat("Filtered to RM levels:", paste(levels_to_include, collapse = ", "), "\n")
          cat("Observations included:", length(include_idx), "\n")
        } else if(!is.null(levels_to_include) && length(levels_to_include) < 2) {
          showNotification("Please select at least 2 levels for comparison!", type = "error", duration = 5)
          return()
        }

        # CRITICAL VALIDATION: Check data structure matches fd_obj
        n_curves_in_fd <- ncol(fd_to_use$coefs)

        # Additional validation: lengths must match
        if(length(subject_id_data) != n_curves_in_fd || length(rm_factor_data) != n_curves_in_fd) {
          showNotification(
            sprintf("ERROR: Variable lengths don't match! Subject ID: %d, RM factor: %d, Curves: %d",
                    length(subject_id_data), length(rm_factor_data), n_curves_in_fd),
            type = "error", duration = 10)
          return()
        }

        cat("Within-subjects validation passed:\n")
        cat("  Number of curves:", n_curves_in_fd, "\n")
        cat("  Subject ID length:", length(subject_id_data), "\n")
        cat("  RM factor length:", length(rm_factor_data), "\n")
        cat("  Unique levels:", length(unique(rm_factor_data)), "\n")

        # P3.5: the estimator no longer reads input$ or calls showNotification
        # itself -- the Shiny pieces are passed in here, so the same function
        # object runs unchanged in the exported script with inert defaults.
        values$fanova_results <- withProgress(
          message = "Repeated-measures functional ANOVA...", value = 0,
          perform_rm_fanova(
            fd_obj = fd_to_use,
            subject_id = subject_id_data,
            rm_factor = rm_factor_data,
            n_permutations = input$n_permutations,
            alpha = input$alpha_level,
            run_global_test = isTRUE(input$rm_global_test),
            progress = function(frac, detail = NULL) incProgress(frac, detail = detail),
            notify = function(msg, type = "message")
              showNotification(msg, type = type,
                               duration = if(identical(type, "message")) 8 else 12)
          ))

        values$fanova_results$design <- "within"
        # P3.8: the exported script referenced subject_id and rm_factor without
        # ever defining them, so the within-design section could not run. Keep
        # the design vectors the fit actually used so the export can write them.
        values$fanova_results$subject_id <- as.character(subject_id_data)
        values$fanova_results$rm_factor  <- as.character(rm_factor_data)
      }
      
      # Store which data source was used
      values$fanova_results$data_source <- input$fanova_data_source

      # MERGED APP: the post-hoc tab reads these back so it compares the same
      # variable, the same levels and the same curves the omnibus did.
      values$fanova_results$fd_used <- fd_to_use
      values$fanova_results$group_var <- if(input$fanova_design == "within") {
        input$rm_factor_var
      } else if(!is.null(input$fanova_group_var) && nzchar(input$fanova_group_var)) {
        input$fanova_group_var
      } else if(!is.null(values$selected_group_vars)) {
        values$selected_group_vars[1]
      } else NULL
      
      showNotification("Functional ANOVA completed!", type = "message", duration = 3)
      
    }, error = function(e) {
      cat("FANOVA error:", e$message, "\n")
      showNotification(paste("FANOVA error:", e$message), type = "error", duration = 10)
    })
  })
  
  # FANOVA outputs (enhanced for both between and within designs)
  output$fanova_global_results <- renderPrint({
    req(values$fanova_results)
    
    res <- values$fanova_results
    
    cat("========================================\n")
    cat("    FUNCTIONAL ANOVA RESULTS\n")
    cat("========================================\n\n")
    
    # Show design type
    design_type <- if(!is.null(res$design)) {
      if(res$design == "between") "Between-Subjects" else "Within-Subjects (Repeated Measures)"
    } else {
      "Between-Subjects"
    }
    cat("Design:", design_type, "\n\n")
    
    # Show group/condition information
    if(!is.null(res$design) && res$design == "within") {
      cat("Number of conditions:", res$n_groups, "\n")
      cat("Conditions:", paste(res$groups, collapse = ", "), "\n")
      cat("Subjects per condition:", paste(res$group_sizes, collapse = ", "), "\n")
    } else {
      cat("Number of groups:", res$n_groups, "\n")
      cat("Groups:", paste(res$groups, collapse = ", "), "\n")
      cat("Group sizes:", paste(res$group_sizes, collapse = ", "), "\n")
    }
    
    cat("\nL2 statistic:", round(res$L2_stat, 4), "\n")
    
    if(!is.na(res$p_value_L2)) {
      cat("L2 p-value:", round(res$p_value_L2, 4), "\n")
    } else if (identical(res$design, "within")) {
      cat("L2 p-value: not computed by the pointwise procedure.\n")
    }

    # P2.6: the optional global repeated-measures test
    if (!is.null(res$rm_global)) {
      cat("\n---------------- GLOBAL TEST (rmfanova) ----------------\n")
      g <- res$rm_global
      if (!isTRUE(g$ok)) {
        cat("Not run:", g$reason, "\n")
      } else {
        cat(sprintf("Complete balanced design: %d of %d subjects (%d condition%s).\n",
                    g$n_complete, g$n_total, length(g$visits),
                    if (length(g$visits) == 1) "" else "s"))
        if (length(g$dropped))
          cat("  Dropped for an incomplete design:",
              paste(head(g$dropped, 12), collapse = ", "),
              if (length(g$dropped) > 12) sprintf("... (%d total)", length(g$dropped)) else "", "\n")
        cat(sprintf("  %d permutations, %d bootstrap replicates.\n", g$n_perm, g$n_boot))
        cat("\nTest statistics:\n")
        print(round(g$stats, 4))
        cat("\nGlobal p-values (permutation scheme P1):\n")
        for (nm in names(g$p_default))
          cat(sprintf("  %-8s %.4f\n", nm, g$p_default[[nm]]))
        cat("\nThese three are shown because a 400-replicate simulation put them\n")
        cat("closest to nominal with full power. Of the other twelve outputs,\n")
        cat(sprintf("  %s rejected true nulls at 17.5%% and 14.2%% (nominal 5%%)\n",
                    paste(FCK_RMFANOVA_INFLATED, collapse = " and ")))
        cat(sprintf("  %s had no power at all in that setting.\n",
                    paste(FCK_RMFANOVA_NOPOWER, collapse = " and ")))
        cat("All fifteen are in the results object as rm_global$p_values.\n")
      }
    }

    cat("\n========================================\n")
  })
  
  output$fanova_summary_table <- renderDT({
    req(values$fanova_results)
    
    res <- values$fanova_results
    
    # Calculate additional statistics for each group/condition
    curves <- eval.fd(res$time_points, values$fd_obj)
    
    # Determine label for first column
    group_label <- if(!is.null(res$design) && res$design == "within") {
      "Condition"
    } else {
      "Group"
    }
    
    summary_df <- data.frame(
      GroupOrCondition = res$groups,
      N = res$group_sizes,
      Mean_Area = round(apply(res$group_means, 2, function(x) mean(x)), 3),
      SD_Area = round(apply(res$group_means, 2, function(x) sd(x)), 3),
      Max_Value = round(apply(res$group_means, 2, max), 3),
      Min_Value = round(apply(res$group_means, 2, min), 3),
      Range = round(apply(res$group_means, 2, function(x) diff(range(x))), 3)
    )
    
    # Rename first column
    colnames(summary_df)[1] <- group_label
    
    datatable(summary_df, 
              options = list(pageLength = 10, dom = 't', scrollX = TRUE), 
              rownames = FALSE) %>%
      formatStyle("N",
                  backgroundColor = styleInterval(c(5, 10), 
                                                  c("#ffcccc", "#ffffcc", "#ccffcc")))
  })
  
  # FANOVA plots - ENHANCED with SD bands and toggles
  output$fanova_mean_plot <- renderPlotly({
    req(values$fanova_results)

    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    # Get toggle values (with defaults if not yet initialized)
    show_sd_bands <- if(!is.null(input$fanova_show_sd_bands)) input$fanova_show_sd_bands else TRUE
    show_sig_regions <- if(!is.null(input$fanova_show_sig_regions)) input$fanova_show_sig_regions else TRUE

    # Create color palette that scales with number of groups
    base_cols <- c("red","blue","green","orange","purple","brown","cyan","magenta","darkgray","gold")
    colors <- colorRampPalette(base_cols)(res$n_groups)

    n_time <- length(res$time_points)

    # Calculate SD for each group
    curves <- eval.fd(res$time_points, values$fd_obj)
    group_sds <- matrix(NA, n_time, res$n_groups)

    for(i in 1:res$n_groups) {
      group_idx <- which(res$group_labels == res$groups[i])
      if(length(group_idx) > 1) {
        group_curves <- curves[, group_idx, drop = FALSE]
        group_sds[,i] <- apply(group_curves, 1, sd)
      } else {
        group_sds[,i] <- 0
      }
    }

    p <- plot_ly(type = 'scatter', mode = 'lines')

    # Add SD bands and mean curves for each group
    # NOTE: No legendgroup linking - each trace toggles independently
    for(i in 1:res$n_groups) {
      # Build safe rgba fill color with alpha
      col_rgb <- col2rgb(colors[i])
      fillcol <- sprintf("rgba(%d,%d,%d,%.2f)", col_rgb[1], col_rgb[2], col_rgb[3], 0.2)

      # Add ±1 SD band (only if toggle is on)
      if(show_sd_bands) {
        p <- p %>% add_trace(
          x = c(res$time_points, rev(res$time_points)),
          y = c(res$group_means[,i] + group_sds[,i],
                rev(res$group_means[,i] - group_sds[,i])),
          fill = 'toself',
          fillcolor = fillcol,
          line = list(color = 'transparent'),
          showlegend = TRUE,
          name = paste0(res$groups[i], " ±SD"),
          hoverinfo = 'skip'
          # No legendgroup - toggles independently from mean curve
        )
      }

      # Add mean curve
      p <- p %>% add_trace(
        x = res$time_points,
        y = res$group_means[,i],
        line = list(color = colors[i], width = 3),
        name = res$groups[i],
        # No legendgroup - toggles independently from SD band
        hovertemplate = paste(res$groups[i], "<br>Time: %{customdata}<br>Mean: %{y:.3f}<extra></extra>"),
        customdata = hover_times
      )
    }

    # Add significant regions as background (only if toggle is on)
    # All significant regions share a legendgroup so they toggle together
    if(show_sig_regions && any(res$sig_regions)) {
      sig_starts <- which(diff(c(0, res$sig_regions)) == 1)
      sig_ends <- which(diff(c(res$sig_regions, 0)) == -1)

      for(idx in 1:length(sig_starts)) {
        y_range <- range(c(res$group_means + group_sds, res$group_means - group_sds))

        p <- p %>% add_trace(
          x = c(res$time_points[sig_starts[idx]], res$time_points[sig_ends[idx]],
                res$time_points[sig_ends[idx]], res$time_points[sig_starts[idx]]),
          y = c(min(y_range), min(y_range), max(y_range), max(y_range)),
          fill = 'toself',
          fillcolor = 'rgba(200, 200, 200, 0.2)',
          line = list(color = 'transparent'),
          showlegend = (idx == 1),
          name = 'Significant',
          legendgroup = 'significant_regions',  # All sig regions toggle together
          hoverinfo = 'skip'
        )
      }
    }
    
    # Dynamic title based on settings
    plot_title <- if(show_sd_bands) "Group Mean Functions ±1 SD" else "Group Mean Functions"

    p <- p %>% layout(
      title = plot_title,
      yaxis = list(title = "Value"),
      hovermode = 'x',
      legend = list(x = 0.02, y = 0.98, tracegroupgap = 5)
    )
    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))

    # Enable editable mode so legend can be dragged
    p <- p %>% config(editable = TRUE)
    p
  })
  
  output$fanova_fstat_plot <- renderPlotly({
    req(values$fanova_results)
    
    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    # Calculate critical value
    crit_val <- qf(1 - res$alpha, res$df_between, res$df_within)
    
    p <- plot_ly(type = 'scatter', mode = 'lines') %>%
      add_trace(x = res$time_points, y = res$F_stat,
                line = list(color = 'blue', width = 2),
                name = 'F-statistic',
                hovertemplate = "Time: %{customdata}<br>F-stat: %{y:.2f}<extra></extra>",
                customdata = hover_times) %>%
      add_trace(x = c(0, 1), y = c(crit_val, crit_val),
                line = list(color = 'red', width = 2, dash = 'dash'),
                name = paste('Critical value (α =', res$alpha, ')'),
                hovertemplate = paste("Critical F =", round(crit_val, 2), "<extra></extra>")) %>%
      layout(title = paste("Pointwise F-statistics (", res$n_groups, "groups)"),
             yaxis = list(title = "F-statistic"))

    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))
    p
  })
  
  output$fanova_pvalue_plot <- renderPlotly({
    req(values$fanova_results)
    
    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    p <- plot_ly(type = 'scatter', mode = 'lines') %>%
      add_trace(x = res$time_points, y = res$p_values_adjusted,
                line = list(color = 'darkgreen', width = 2),
                name = 'Adjusted p-values (FDR)',
                hovertemplate = "Time: %{customdata}<br>P-value: %{y:.4f}<extra></extra>",
                customdata = hover_times) %>%
      add_trace(x = res$time_points, y = res$p_values_pointwise,
                line = list(color = 'lightgreen', width = 1, dash = 'dot'),
                name = 'Raw p-values',
                hovertemplate = "Time: %{customdata}<br>P-value: %{y:.4f}<extra></extra>",
                customdata = hover_times) %>%
      add_trace(x = c(0, 1), y = c(res$alpha, res$alpha),
                line = list(color = 'red', width = 2, dash = 'dash'),
                name = paste('α =', res$alpha),
                hovertemplate = paste("Alpha =", res$alpha, "<extra></extra>")) %>%
      add_trace(x = c(0, 1), y = c(0.01, 0.01),
                line = list(color = 'orange', width = 1, dash = 'dot'),
                name = 'p = 0.01',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      add_trace(x = c(0, 1), y = c(0.001, 0.001),
                line = list(color = 'darkred', width = 1, dash = 'dot'),
                name = 'p = 0.001',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      layout(title = paste("Pointwise p-values (", res$n_groups, "groups)"),
             yaxis = list(title = "p-value", type = 'log',
                          range = c(log10(0.0001), log10(1))))
    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))
    p
  })
  
  output$fanova_effect_size_plot <- renderPlotly({
    req(values$fanova_results)
    
    res <- values$fanova_results
    hover_times <- hover_time_labels(res$time_points)

    # P3.3: a repeated-measures design reports PARTIAL eta-squared, an
    # independent-groups design classical eta-squared. They are different
    # quantities on the same axis, so the label has to say which one this is.
    eta_partial <- identical(res$eta_squared_type, "partial")
    eta_lab     <- if(eta_partial) "partial \u03b7\u00b2" else "\u03b7\u00b2"

    # Add background for effect size interpretation
    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    # Add reference lines for effect size interpretation
    p <- p %>% 
      add_trace(x = c(0, 1), y = c(0.01, 0.01),
                line = list(color = 'lightgray', width = 1, dash = 'dot'),
                name = 'Small effect',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      add_trace(x = c(0, 1), y = c(0.06, 0.06),
                line = list(color = 'gray', width = 1, dash = 'dot'),
                name = 'Medium effect',
                showlegend = FALSE,
                hoverinfo = 'skip') %>%
      add_trace(x = c(0, 1), y = c(0.14, 0.14),
                line = list(color = 'darkgray', width = 1, dash = 'dot'),
                name = 'Large effect',
                showlegend = FALSE,
                hoverinfo = 'skip')
    
    # Add eta-squared curve
    p <- p %>% add_trace(x = res$time_points, y = res$eta_squared,
                         fill = 'tozeroy',
                         fillcolor = 'rgba(100, 100, 255, 0.3)',
                         line = list(color = 'purple', width = 2),
                         name = paste(eta_lab, '(Effect size)'),
                         hovertemplate = paste0("Time: %{customdata}<br>", eta_lab,
                                                ": %{y:.3f}<extra></extra>"),
                         customdata = hover_times)

    p <- p %>% layout(title = list(text = paste0(
                        "Effect Size (", eta_lab, ") across Time \u2014 ",
                        res$n_groups, if(eta_partial) " conditions" else " groups",
                        if(eta_partial)
                          "<br><sub>SS<sub>condition</sub> / (SS<sub>condition</sub> + SS<sub>error</sub>), from the same repeated-measures decomposition as F</sub>"
                        else
                          "<br><sub>SS<sub>between</sub> / SS<sub>total</sub></sub>")),
                      yaxis = list(title = paste(eta_lab, "(proportion of variance explained)"),
                                   range = c(0, max(0.3, max(res$eta_squared, na.rm = TRUE) * 1.1))),
                      annotations = list(
                        list(x = 0.95, y = 0.01, text = "Small", showarrow = FALSE, 
                             font = list(size = 10, color = "gray")),
                        list(x = 0.95, y = 0.06, text = "Medium", showarrow = FALSE,
                             font = list(size = 10, color = "gray")),
                        list(x = 0.95, y = 0.14, text = "Large", showarrow = FALSE,
                             font = list(size = 10, color = "gray"))
                      ))

    p <- format_plotly_time_axis(p, res$time_points, tick_step_hours = as.numeric(input$tick_freq_fanova))
    p
  })
  
  output$fanova_effect_summary <- renderPrint({
    req(values$fanova_results)
    
    res <- values$fanova_results
    eta_partial <- identical(res$eta_squared_type, "partial")
    eta_lab     <- if(eta_partial) "partial \u03b7\u00b2" else "\u03b7\u00b2"
    mean_eta    <- mean(res$eta_squared, na.rm = TRUE)

    cat("Effect Size Summary (", eta_lab, "):\n", sep = "")
    cat("Mean ", eta_lab, ": ", round(mean_eta, 4), "\n", sep = "")
    cat("Max  ", eta_lab, ": ", round(max(res$eta_squared, na.rm = TRUE), 4),
        " at ", hover_time_labels(res$time_points)[which.max(res$eta_squared)], "\n", sep = "")
    if(any(!is.finite(res$eta_squared)))
      cat("Not computable at ", sum(!is.finite(res$eta_squared)),
          " of ", length(res$eta_squared),
          " time points (fewer than two subjects with a complete design there).\n", sep = "")

    if(eta_partial) {
      cat("\nDefinition: SS_condition / (SS_condition + SS_error), from the same\n")
      cat("repeated-measures decomposition the F test used -- equivalently\n")
      cat("df_cond * F / (df_cond * F + df_resid). Between-subject differences are\n")
      cat("partialled out of the denominator, as they are out of the F test.\n")
      if(!is.null(res$eta_squared_classical)) {
        cat("\nFor comparison, classical eta-squared (SS_between / SS_total over all\n")
        cat("curves, subject margin NOT removed) has mean ",
            round(mean(res$eta_squared_classical, na.rm = TRUE), 4), ".\n", sep = "")
        cat("Versions of this app before the P3.3 correction reported THAT number\n")
        cat("under a repeated-measures F. It is biased downward here and is shown\n")
        cat("only so earlier output can be reconciled; do not report it.\n")
      }
    } else {
      cat("\nDefinition: SS_between / SS_total, the classical eta-squared for an\n")
      cat("independent-groups design.\n")
    }
    cat("\nConventional benchmarks (Cohen): 0.01 small, 0.06 medium, 0.14 large.\n")
  })
  
  # Check if FANOVA completed
  output$fanova_completed <- reactive({
    !is.null(values$fanova_results)
  })
  outputOptions(output, "fanova_completed", suspendWhenHidden = FALSE)
  
  # Pairwise comparison functions - ENHANCED
  # Repeated Measures Pairwise Comparisons (for within-subjects designs)
  perform_pairwise_comparisons_rm <- function(fd_obj, subject_id, rm_factor, n_permutations = 200,
                                              correction_method = "bonferroni", alpha = 0.05) {
    
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    curves <- eval.fd(time_points, fd_obj)
    
    subject_id <- as.factor(subject_id)
    rm_factor <- as.factor(rm_factor)
    
    conditions <- levels(rm_factor)
    n_conditions <- length(conditions)
    
    n_pairs <- choose(n_conditions, 2)
    pairs <- combn(conditions, 2, simplify = FALSE)
    pair_names <- sapply(pairs, function(p) paste(p[1], "vs", p[2]))
    
    pairwise_results <- list()
    
    cat("Performing paired comparisons (within-subjects)...\n")
    cat("Number of condition pairs:", n_pairs, "\n")
    
    withProgress(message = 'Performing pairwise comparisons', value = 0, {
      for(pair_idx in 1:n_pairs) {
        incProgress(1/n_pairs, detail = paste("Comparing", pair_names[pair_idx]))
        
        pair <- pairs[[pair_idx]]
        
        # Get indices for each condition
        idx1 <- which(rm_factor == pair[1])
        idx2 <- which(rm_factor == pair[2])
        
        # Match subjects across conditions
        subjects_in_1 <- subject_id[idx1]
        subjects_in_2 <- subject_id[idx2]
        
        # Find subjects present in both conditions (for paired comparison)
        common_subjects <- intersect(subjects_in_1, subjects_in_2)
        n_pairs_subj <- length(common_subjects)
        
        cat("  Pair:", pair_names[pair_idx], "- Matched subjects:", n_pairs_subj, "\n")
        
        if(n_pairs_subj < 2) {
          warning("Not enough paired subjects for comparison: ", pair_names[pair_idx])
          next
        }
        
        # Extract matched curves
        matched_curves1 <- matrix(NA, n_time, n_pairs_subj)
        matched_curves2 <- matrix(NA, n_time, n_pairs_subj)
        
        for(i in 1:n_pairs_subj) {
          subj <- common_subjects[i]
          
          # Find this subject's curve in condition 1
          subj_idx1 <- idx1[which(subjects_in_1 == subj)[1]]
          matched_curves1[, i] <- curves[, subj_idx1]
          
          # Find this subject's curve in condition 2
          subj_idx2 <- idx2[which(subjects_in_2 == subj)[1]]
          matched_curves2[, i] <- curves[, subj_idx2]
        }
        
        # Calculate paired differences
        paired_diffs <- matched_curves1 - matched_curves2
        mean_diff <- rowMeans(paired_diffs)
        
        # Calculate paired t-statistics
        se_diff <- apply(paired_diffs, 1, sd) / sqrt(n_pairs_subj)
        t_stat <- mean_diff / se_diff
        
        # Handle NaN/Inf (when se_diff is 0)
        t_stat[!is.finite(t_stat)] <- 0
        
        # L2 norm statistic
        L2_stat <- fck_l2_norm(mean_diff, time_points)
        
        # Permutation test for paired data
        # Randomly flip signs of differences
        t_stat_perm <- matrix(NA, n_time, n_permutations)
        L2_stat_perm <- numeric(n_permutations)
        
        for(perm in 1:n_permutations) {
          # Random sign flips for each subject
          sign_flips <- sample(c(-1, 1), n_pairs_subj, replace = TRUE)
          
          perm_diffs <- paired_diffs * rep(sign_flips, each = n_time)
          perm_mean_diff <- rowMeans(perm_diffs)
          perm_se_diff <- apply(perm_diffs, 1, sd) / sqrt(n_pairs_subj)
          
          t_stat_perm[, perm] <- perm_mean_diff / perm_se_diff
          t_stat_perm[!is.finite(t_stat_perm[, perm]), perm] <- 0
          
          L2_stat_perm[perm] <- fck_l2_norm(perm_mean_diff, time_points)
        }
        
        # Calculate p-values
        p_values_pointwise <- numeric(n_time)
        for(t in 1:n_time) {
          p_values_pointwise[t] <- mean(abs(t_stat_perm[t, ]) >= abs(t_stat[t]), na.rm = TRUE)
        }
        
        p_value_L2 <- (1 + sum(L2_stat_perm >= L2_stat, na.rm = TRUE)) /
                      (1 + sum(is.finite(L2_stat_perm)))
        
        # Bootstrap confidence intervals for paired differences
        # AUDIT (P1.2): 100 replicates gives a 2.5% quantile estimated from the
      # 2nd-3rd order statistic. These are POINTWISE percentile intervals, not
      # simultaneous functional bands -- see the label below.
      n_boot <- 2000
        diff_boot <- matrix(NA, n_time, n_boot)
        
        for(boot in 1:n_boot) {
          boot_idx <- sample(1:n_pairs_subj, replace = TRUE)
          diff_boot[, boot] <- rowMeans(paired_diffs[, boot_idx, drop = FALSE])
        }
        
        ci_lower <- apply(diff_boot, 1, quantile, probs = 0.025, na.rm = TRUE)
        ci_upper <- apply(diff_boot, 1, quantile, probs = 0.975, na.rm = TRUE)
        
        # Cohen's d for paired samples (using SD of differences)
        sd_diff <- apply(paired_diffs, 1, sd)
        cohens_d <- mean_diff / sd_diff
        cohens_d[!is.finite(cohens_d)] <- 0
        
        # Calculate means for each condition (for plotting)
        mean1 <- rowMeans(matched_curves1)
        mean2 <- rowMeans(matched_curves2)
        
        pairwise_results[[pair_names[pair_idx]]] <- list(
          group1 = pair[1],
          group2 = pair[2],
          n1 = n_pairs_subj,  # Number of matched pairs
          n2 = n_pairs_subj,
          n_matched = n_pairs_subj,
          design = "within",
          mean1 = mean1,
          mean2 = mean2,
          mean_diff = mean_diff,
          t_stat = t_stat,
          p_values_pointwise = p_values_pointwise,
          L2_stat = L2_stat,
          p_value_L2 = p_value_L2,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          cohens_d = cohens_d,
          se_diff = se_diff
        )
      }
    })
    
    # Apply multiple comparison correction
    all_p_values_L2 <- sapply(pairwise_results, function(x) x$p_value_L2)
    adjusted_p_values_L2 <- p.adjust(all_p_values_L2, method = correction_method)
    
    # Apply correction to pointwise p-values
    for(i in 1:length(pairwise_results)) {
      pairwise_results[[i]]$p_values_adjusted <- p.adjust(pairwise_results[[i]]$p_values_pointwise, 
                                                          method = correction_method)
      pairwise_results[[i]]$p_value_L2_adjusted <- adjusted_p_values_L2[i]
      pairwise_results[[i]]$sig_regions <- pairwise_results[[i]]$p_values_adjusted < alpha
      pairwise_results[[i]]$sig_global <- pairwise_results[[i]]$p_value_L2_adjusted < alpha
    }
    
    return(list(
      results = pairwise_results,
      time_points = time_points,
      correction_method = correction_method,
      alpha = alpha,
      n_permutations = n_permutations,
      groups = conditions,
      n_groups = n_conditions,
      pair_names = pair_names,
      design = "within"
    ))
  }
  
  # Between-Subjects Pairwise Comparisons (original function)
  perform_pairwise_comparisons <- function(fd_obj, group_labels, n_permutations = 200,
                                           correction_method = "bonferroni", alpha = 0.05) {
    
    n_curves <- ncol(fd_obj$coefs)
    n_time <- 100
    time_points <- seq(0, 1, length.out = n_time)
    
    curves <- eval.fd(time_points, fd_obj)
    
    group_labels <- as.factor(group_labels)
    groups <- levels(group_labels)
    n_groups <- length(groups)
    
    n_pairs <- choose(n_groups, 2)
    pairs <- combn(groups, 2, simplify = FALSE)
    pair_names <- sapply(pairs, function(p) paste(p[1], "vs", p[2]))
    
    pairwise_results <- list()
    
    withProgress(message = 'Performing pairwise comparisons', value = 0, {
      for(pair_idx in 1:n_pairs) {
        incProgress(1/n_pairs, detail = paste("Comparing", pair_names[pair_idx]))
        
        pair <- pairs[[pair_idx]]
        idx1 <- which(group_labels == pair[1])
        idx2 <- which(group_labels == pair[2])
        
        curves1 <- curves[, idx1, drop = FALSE]
        curves2 <- curves[, idx2, drop = FALSE]
        
        n1 <- length(idx1)
        n2 <- length(idx2)
        
        mean1 <- rowMeans(curves1)
        mean2 <- rowMeans(curves2)
        mean_diff <- mean1 - mean2
        
        # Calculate proper t-statistics
        pooled_var <- ((n1 - 1) * apply(curves1, 1, var) + 
                         (n2 - 1) * apply(curves2, 1, var)) / (n1 + n2 - 2)
        se_diff <- sqrt(pooled_var * (1/n1 + 1/n2))
        t_stat <- mean_diff / se_diff
        
        # L2 norm statistic
        L2_stat <- fck_l2_norm(mean_diff, time_points)
        
        # Permutation test
        t_stat_perm <- matrix(NA, n_time, n_permutations)
        L2_stat_perm <- numeric(n_permutations)
        
        combined_curves <- cbind(curves1, curves2)
        combined_labels <- c(rep(1, n1), rep(2, n2))
        
        for(perm in 1:n_permutations) {
          perm_labels <- sample(combined_labels)
          
          perm_curves1 <- combined_curves[, perm_labels == 1, drop = FALSE]
          perm_curves2 <- combined_curves[, perm_labels == 2, drop = FALSE]
          
          perm_mean1 <- rowMeans(perm_curves1)
          perm_mean2 <- rowMeans(perm_curves2)
          perm_diff <- perm_mean1 - perm_mean2
          
          perm_pooled_var <- ((n1 - 1) * apply(perm_curves1, 1, var) + 
                                (n2 - 1) * apply(perm_curves2, 1, var)) / (n1 + n2 - 2)
          perm_se_diff <- sqrt(perm_pooled_var * (1/n1 + 1/n2))
          
          t_stat_perm[, perm] <- perm_diff / perm_se_diff
          L2_stat_perm[perm] <- fck_l2_norm(perm_diff, time_points)
        }
        
        # Calculate p-values
        p_values_pointwise <- numeric(n_time)
        for(t in 1:n_time) {
          p_values_pointwise[t] <- mean(abs(t_stat_perm[t, ]) >= abs(t_stat[t]), na.rm = TRUE)
        }
        
        p_value_L2 <- (1 + sum(L2_stat_perm >= L2_stat, na.rm = TRUE)) /
                      (1 + sum(is.finite(L2_stat_perm)))
        
        # Bootstrap confidence intervals
        # AUDIT (P1.2): 100 replicates gives a 2.5% quantile estimated from the
      # 2nd-3rd order statistic. These are POINTWISE percentile intervals, not
      # simultaneous functional bands -- see the label below.
      n_boot <- 2000
        diff_boot <- matrix(NA, n_time, n_boot)
        
        for(boot in 1:n_boot) {
          boot_idx1 <- sample(idx1, replace = TRUE)
          boot_idx2 <- sample(idx2, replace = TRUE)
          
          boot_mean1 <- rowMeans(curves[, boot_idx1, drop = FALSE])
          boot_mean2 <- rowMeans(curves[, boot_idx2, drop = FALSE])
          diff_boot[, boot] <- boot_mean1 - boot_mean2
        }
        
        ci_lower <- apply(diff_boot, 1, quantile, probs = 0.025)
        ci_upper <- apply(diff_boot, 1, quantile, probs = 0.975)
        
        # Cohen's d effect size
        cohens_d <- mean_diff / sqrt(pooled_var)
        
        pairwise_results[[pair_names[pair_idx]]] <- list(
          group1 = pair[1],
          group2 = pair[2],
          n1 = n1,
          n2 = n2,
          design = "between",
          mean_diff = mean_diff,
          t_stat = t_stat,
          p_values_pointwise = p_values_pointwise,
          L2_stat = L2_stat,
          p_value_L2 = p_value_L2,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          cohens_d = cohens_d,
          se_diff = se_diff
        )
      }
    })
    
    # Apply multiple comparison correction
    all_p_values_L2 <- sapply(pairwise_results, function(x) x$p_value_L2)
    adjusted_p_values_L2 <- p.adjust(all_p_values_L2, method = correction_method)
    
    # Apply correction to pointwise p-values
    for(i in 1:n_pairs) {
      pairwise_results[[i]]$p_values_adjusted <- p.adjust(pairwise_results[[i]]$p_values_pointwise, 
                                                          method = correction_method)
      pairwise_results[[i]]$p_value_L2_adjusted <- adjusted_p_values_L2[i]
      pairwise_results[[i]]$sig_regions <- pairwise_results[[i]]$p_values_adjusted < alpha
      pairwise_results[[i]]$sig_global <- pairwise_results[[i]]$p_value_L2_adjusted < alpha
    }
    
    return(list(
      results = pairwise_results,
      time_points = time_points,
      correction_method = correction_method,
      alpha = alpha,
      n_permutations = n_permutations,
      groups = groups,
      n_groups = n_groups,
      pair_names = pair_names,
      design = "between"
    ))
  }
  
  # Run pairwise comparisons - ENHANCED for both between and within designs
  observeEvent(input$run_pairwise, {
    if(is.null(values$fanova_results)) {
      showNotification("Please run Functional ANOVA first!", type = "error", duration = 5)
      return()
    }
    
    cat("Running pairwise comparisons...\n")
    
    # Check if this is a within-subjects design
    design_type <- if(!is.null(values$fanova_results$design)) {
      values$fanova_results$design
    } else {
      "between"
    }
    
    cat("Design type:", design_type, "\n")
    
    tryCatch({
      # Use the same data source as FANOVA
      fd_to_use <- NULL
      
      # Check what data source was used in FANOVA
      if(!is.null(values$fanova_results$data_source) && 
         values$fanova_results$data_source == "warped" && 
         !is.null(values$warping_results)) {
        
        cat("Using time-warped curves for pairwise comparisons\n")
        
        if(!is.null(values$warping_results$regfd)) {
          fd_to_use <- values$warping_results$regfd
        } else if(!is.null(values$warping_results$registered_curves)) {
          n_time <- nrow(values$warping_results$registered_curves)
          time_points <- values$warping_results$time_points
          if(is.null(time_points)) {
            time_points <- seq(0, 1, length.out = n_time)
          }
          basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = min(20, n_time-2))
          # Use lambda=0: registered curves already processed, just need fd representation
          fd_to_use <- smooth.basis(time_points, values$warping_results$registered_curves, 
                                    fdPar(basis, 2, 0))$fd
        }
      }
      
      # Fallback to original data if warped not available
      if(is.null(fd_to_use)) {
        cat("Using original curves for pairwise comparisons\n")
        fd_to_use <- values$fd_obj
      }
      
      # Call appropriate pairwise function based on design
      # IMPORTANT: Use the same n_permutations that was used in the omnibus FANOVA
      n_perm_to_use <- if(!is.null(values$fanova_results$n_permutations)) {
        values$fanova_results$n_permutations
      } else {
        input$pairwise_permutations  # Fallback
      }
      
      cat("Using", n_perm_to_use, "permutations (same as omnibus test)\n")
      
      # MERGED APP: what to compare is resolved in one place
      # (FCK/server/52_posthoc_source.R). Before this, the between-subjects
      # branch passed values$group_labels -- the FIRST scalar variable selected
      # at import -- while the omnibus test had run on input$fanova_group_var,
      # optionally restricted to a subset of its levels on a matching fd
      # object. With more than one scalar variable selected, the omnibus and
      # its "post-hoc" tests were silently testing different variables.
      spec <- fck_posthoc_spec(input, values)
      if (!isTRUE(spec$ok)) {
        showNotification(spec$message, type = "error", duration = 10)
        return()
      }
      fd_for_pairs <- if (!is.null(spec$fd)) spec$fd else fd_to_use
      cat("Post-hoc source:", spec$description, "\n")

      if(spec$design == "within") {
        cat("Performing PAIRED comparisons for within-subjects design\n")

        values$pairwise_results <- perform_pairwise_comparisons_rm(
          fd_obj = fd_for_pairs,
          subject_id = spec$subject_id,
          rm_factor = spec$rm_factor,
          n_permutations = n_perm_to_use,
          correction_method = input$pairwise_correction,
          alpha = input$pairwise_alpha
        )
      } else {
        cat("Performing INDEPENDENT comparisons for between-subjects design\n")

        values$pairwise_results <- perform_pairwise_comparisons(
          fd_obj = fd_for_pairs,
          group_labels = spec$group_labels,
          n_permutations = n_perm_to_use,
          correction_method = input$pairwise_correction,
          alpha = input$pairwise_alpha
        )
      }

      # Travels with the result, so the summary and the export can say what was
      # compared rather than leaving the reader to assume it was the omnibus.
      values$pairwise_results$posthoc_source <- spec$description
      values$pairwise_results$matches_omnibus <- isTRUE(spec$matches_omnibus)
      
      showNotification("Pairwise comparisons completed!", type = "message", duration = 3)
      
    }, error = function(e) {
      cat("Pairwise error:", e$message, "\n")
      showNotification(paste("Pairwise comparison error:", e$message), type = "error", duration = 10)
    })
  })

  # ALL REMAINING PAIRWISE AND EXPORT OUTPUTS - included as-is from the original
  # (These are all already correctly written - just ensuring they're included)
  
  output$pairwise_selector <- renderUI({
    req(values$pairwise_results)
    
    selectInput("selected_pair", "Select pairwise comparison:",
                choices = values$pairwise_results$pair_names,
                selected = values$pairwise_results$pair_names[1])
  })
  
  output$pairwise_summary <- renderPrint({
    req(values$pairwise_results)
    
    res <- values$pairwise_results
    
    cat("========================================\n")
    cat("  PAIRWISE COMPARISONS SUMMARY\n")
    cat("========================================\n\n")
    
    # Show design type
    design_type <- if(!is.null(res$design)) {
      if(res$design == "within") "Within-Subjects (PAIRED)" else "Between-Subjects (INDEPENDENT)"
    } else {
      "Between-Subjects (INDEPENDENT)"
    }
    cat("Design:", design_type, "\n\n")
    
    cat("Number of groups/conditions:", res$n_groups, "\n")
    cat("Number of comparisons:", length(res$pair_names), "\n")
    cat("Correction method:", res$correction_method, "\n")
    cat("Significance level:", res$alpha, "\n")
    cat("Permutations:", res$n_permutations, "\n\n")
    
    n_sig_global <- sum(sapply(res$results, function(x) x$sig_global))
    
    cat("Results Summary:\n")
    cat("----------------\n")
    cat("Globally significant comparisons:", n_sig_global, "/", 
        length(res$pair_names), 
        sprintf("(%.1f%%)", 100 * n_sig_global / length(res$pair_names)), "\n\n")
    
    if(n_sig_global > 0) {
      cat("Significant pairs (global test):\n")
      for(pair_name in res$pair_names) {
        if(res$results[[pair_name]]$sig_global) {
          p_val <- res$results[[pair_name]]$p_value_L2_adjusted
          stars <- if(p_val < 0.001) "***" else if(p_val < 0.01) "**" else if(p_val < 0.05) "*" else ""
          cat(sprintf("   - %-20s p = %s %s\n", 
                      pair_name, 
                      format.pval(p_val, digits = 4),
                      stars))
        }
      }
    } else {
      cat("No significant pairwise differences found.\n")
    }
    
    cat("\n========================================\n")
  })
  
  output$pairwise_global_table <- renderDT({
    req(values$pairwise_results)
    
    results_df <- data.frame(
      Comparison = values$pairwise_results$pair_names,
      N1 = sapply(values$pairwise_results$results, function(x) x$n1),
      N2 = sapply(values$pairwise_results$results, function(x) x$n2),
      L2_Stat = round(sapply(values$pairwise_results$results, function(x) x$L2_stat), 3),
      P_Raw = sapply(values$pairwise_results$results, function(x) 
        format.pval(x$p_value_L2, digits = 4)),
      P_Adj = sapply(values$pairwise_results$results, function(x) 
        format.pval(x$p_value_L2_adjusted, digits = 4)),
      Sig = ifelse(sapply(values$pairwise_results$results, function(x) x$sig_global), 
                   "Yes", "No"),
      Mean_Cohen_d = round(sapply(values$pairwise_results$results, function(x) 
        mean(abs(x$cohens_d))), 3)
    )
    
    datatable(results_df, 
              options = list(pageLength = 15, scrollX = TRUE), 
              rownames = FALSE) %>%
      formatStyle("Sig",
                  backgroundColor = styleEqual("Yes", "#d4edda")) %>%
      formatStyle("Mean_Cohen_d",
                  backgroundColor = styleInterval(c(0.2, 0.5, 0.8), 
                                                  c("white", "#ffffcc", "#ffcccc", "#ff9999")))
  })
  
  output$pairwise_difference_plot <- renderPlotly({
    req(values$pairwise_results, input$selected_pair)

    pair_result <- values$pairwise_results$results[[input$selected_pair]]
    time_points <- values$pairwise_results$time_points
    hover_times <- hover_time_labels(time_points)

    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    if(any(pair_result$sig_regions)) {
      sig_starts <- which(diff(c(0, pair_result$sig_regions)) == 1)
      sig_ends <- which(diff(c(pair_result$sig_regions, 0)) == -1)
      
      for(i in 1:length(sig_starts)) {
        y_range <- range(c(pair_result$ci_lower, pair_result$ci_upper))
        y_expand <- diff(y_range) * 0.1
        
        p <- p %>% add_trace(
          x = c(time_points[sig_starts[i]], time_points[sig_ends[i]], 
                time_points[sig_ends[i]], time_points[sig_starts[i]]),
          y = c(min(y_range) - y_expand, min(y_range) - y_expand,
                max(y_range) + y_expand, max(y_range) + y_expand),
          fill = 'toself',
          fillcolor = 'rgba(255, 200, 200, 0.3)',
          line = list(color = 'transparent'),
          showlegend = (i == 1),
          name = 'Significant',
          hoverinfo = 'skip'
        )
      }
    }
    
    if(input$pairwise_confidence_bands) {
      p <- p %>% add_trace(
        x = c(time_points, rev(time_points)),
        y = c(pair_result$ci_lower, rev(pair_result$ci_upper)),
        fill = 'toself',
        fillcolor = 'rgba(100, 100, 255, 0.2)',
        line = list(color = 'transparent'),
        name = '95% pointwise CI',   # P1.2: pointwise, not simultaneous
        hoverinfo = 'skip'
      )
      
      p <- p %>% add_trace(
        x = time_points,
        y = pair_result$ci_lower,
        line = list(color = 'lightblue', width = 1, dash = 'dot'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      ) %>% add_trace(
        x = time_points,
        y = pair_result$ci_upper,
        line = list(color = 'lightblue', width = 1, dash = 'dot'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      )
    }
    
    p <- p %>% add_trace(
      x = time_points,
      y = pair_result$mean_diff,
      line = list(color = 'blue', width = 3),
      name = 'Mean Difference',
      hovertemplate = paste("Time: %{customdata}",
                            "<br>Difference: %{y:.3f}",
                            "<br>Cohen's d: ", round(pair_result$cohens_d, 2),
                            "<extra></extra>"),
      customdata = hover_times
    )
    
    p <- p %>% add_trace(
      x = c(0, 1),
      y = c(0, 0),
      line = list(color = 'black', width = 1, dash = 'dash'),
      name = 'Zero',
      hoverinfo = 'none'
    )
    
    mean_effect <- mean(abs(pair_result$cohens_d))
    effect_text <- if(mean_effect < 0.2) "Negligible" else if(mean_effect < 0.5) "Small" else if(mean_effect < 0.8) "Medium" else "Large"
    
    p <- p %>% layout(
      title = paste("Difference:", pair_result$group1, "-", pair_result$group2,
                    "<br><sub>Mean |Cohen's d| =", round(mean_effect, 2), 
                    "(", effect_text, "effect)</sub>"),
      yaxis = list(title = "Mean Difference"),
      hovermode = 'x',
      legend = list(x = 0.02, y = 0.98)
    )

    p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_pairwise))
    p
  })
  
  output$pairwise_pvalue_plot <- renderPlotly({
    req(values$pairwise_results, input$selected_pair)

    pair_result <- values$pairwise_results$results[[input$selected_pair]]
    time_points <- values$pairwise_results$time_points
    hover_times <- hover_time_labels(time_points)
    
    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    if(any(pair_result$sig_regions)) {
      sig_starts <- which(diff(c(0, pair_result$sig_regions)) == 1)
      sig_ends <- which(diff(c(pair_result$sig_regions, 0)) == -1)
      
      for(i in 1:length(sig_starts)) {
        p <- p %>% add_trace(
          x = c(time_points[sig_starts[i]], time_points[sig_ends[i]], 
                time_points[sig_ends[i]], time_points[sig_starts[i]]),
          y = c(0.0001, 0.0001, 1, 1),
          fill = 'toself',
          fillcolor = 'rgba(200, 255, 200, 0.3)',
          line = list(color = 'transparent'),
          showlegend = (i == 1),
          name = 'Significant region',
          hoverinfo = 'skip'
        )
      }
    }
    
    p <- p %>% add_trace(
      x = time_points,
      y = pair_result$p_values_adjusted,
      line = list(color = 'darkgreen', width = 2),
      name = paste('Adjusted p-values (', values$pairwise_results$correction_method, ')', sep = ''),
      hovertemplate = "Time: %{customdata}<br>Adjusted p: %{y:.4f}<extra></extra>",
      customdata = hover_times
    ) %>%
      add_trace(
        x = time_points,
        y = pair_result$p_values_pointwise,
        line = list(color = 'lightgreen', width = 1, dash = 'dot'),
        name = 'Raw p-values',
        hovertemplate = "Time: %{customdata}<br>Raw p: %{y:.4f}<extra></extra>",
        customdata = hover_times
      ) %>%
      add_trace(
        x = c(0, 1),
        y = c(input$pairwise_alpha, input$pairwise_alpha),
        line = list(color = 'red', width = 2, dash = 'dash'),
        name = paste('α =', input$pairwise_alpha),
        hoverinfo = 'none'
      )
    
    p <- p %>% add_trace(
      x = c(0, 1),
      y = c(0.01, 0.01),
      line = list(color = 'orange', width = 1, dash = 'dot'),
      name = 'p = 0.01',
      hoverinfo = 'none'
    ) %>%
      add_trace(
        x = c(0, 1),
        y = c(0.001, 0.001),
        line = list(color = 'darkred', width = 1, dash = 'dot'),
        name = 'p = 0.001',
        hoverinfo = 'none'
      )
    
    p <- p %>% layout(
      title = paste("P-values:", input$selected_pair),
      yaxis = list(
        title = "p-value",
        type = 'log',
        range = c(log10(0.0001), log10(1)),
        tickvals = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
        ticktext = c("0.001", "0.01", "0.05", "0.1", "0.5", "1")
      ),
      hovermode = 'x',
      legend = list(x = 0.02, y = 0.02)
    )

    p <- format_plotly_time_axis(p, time_points, tick_step_hours = as.numeric(input$tick_freq_pairwise))
    p
  })
  
  output$pairwise_heatmap <- renderPlotly({
    req(values$pairwise_results)
    
    n_groups <- values$pairwise_results$n_groups
    groups <- values$pairwise_results$groups
    
    p_matrix <- matrix(1, nrow = n_groups, ncol = n_groups)
    rownames(p_matrix) <- groups
    colnames(p_matrix) <- groups
    
    for(pair_name in values$pairwise_results$pair_names) {
      pair_result <- values$pairwise_results$results[[pair_name]]
      i <- which(groups == pair_result$group1)
      j <- which(groups == pair_result$group2)
      p_val <- pair_result$p_value_L2_adjusted
      p_matrix[i, j] <- p_val
      p_matrix[j, i] <- p_val
    }
    
    z_matrix <- -log10(p_matrix)
    z_matrix[is.infinite(z_matrix)] <- 5
    
    hover_text <- matrix("", nrow = n_groups, ncol = n_groups)
    for(i in 1:n_groups) {
      for(j in 1:n_groups) {
        if(i == j) {
          hover_text[i, j] <- paste(groups[i], "(same group)")
        } else {
          p_val_text <- if(p_matrix[i, j] < 0.001) {
            "< 0.001"
          } else {
            format(round(p_matrix[i, j], 4), nsmall = 4)
          }
          sig_text <- if(p_matrix[i, j] < values$pairwise_results$alpha) "***" else "n.s."
          hover_text[i, j] <- paste(groups[i], "vs", groups[j],
                                    "<br>p =", p_val_text,
                                    "<br>", sig_text)
        }
      }
    }
    
    annotations <- list()
    for(i in 1:n_groups) {
      for(j in 1:n_groups) {
        if(i != j) {
          star_text <- if(p_matrix[i, j] < 0.001) "***" else 
            if(p_matrix[i, j] < 0.01) "**" else 
              if(p_matrix[i, j] < 0.05) "*" else ""
          
          if(star_text != "") {
            annotations <- append(annotations, list(list(
              x = groups[j],
              y = groups[i],
              text = star_text,
              showarrow = FALSE,
              font = list(color = 'white', size = 14)
            )))
          }
        }
      }
    }
    
    plot_ly(
      z = z_matrix,
      x = groups,
      y = groups,
      type = 'heatmap',
      colorscale = list(
        c(0, 'white'),
        c(0.3, 'lightblue'),
        c(0.6, 'blue'),
        c(1, 'darkred')
      ),
      hovertemplate = "%{text}<extra></extra>",
      text = hover_text,
      colorbar = list(
        title = "Significance",
        tickmode = "array",
        tickvals = c(0, -log10(0.05), -log10(0.01), -log10(0.001), 5),
        ticktext = c("1", "0.05", "0.01", "0.001", "<0.00001")
      )
    ) %>%
      layout(
        title = paste("Pairwise Comparison P-values (", values$pairwise_results$correction_method, "correction)"),
        xaxis = list(title = "", tickangle = 45),
        yaxis = list(title = "", autorange = 'reversed'),
        annotations = annotations
      )
  })
  
  output$pairwise_significance_timeline <- renderPlotly({
    req(values$pairwise_results)
    
    time_points <- values$pairwise_results$time_points
    n_pairs <- length(values$pairwise_results$pair_names)
    hover_times <- hover_time_labels(time_points)
    
    p <- plot_ly(type = 'scatter', mode = 'lines')
    
    base_cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00',
                   '#ffff33','#a65628','#f781bf','#999999','#66c2a5',
                   '#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f')
    colors <- colorRampPalette(base_cols)(n_pairs)
    
    for(i in 1:n_pairs) {
      p <- p %>% add_trace(
        x = time_points,
        y = rep(i, length(time_points)),
        line = list(color = 'lightgray', width = 0.5),
        showlegend = FALSE,
        hoverinfo = 'skip'
      )
    }
    
    for(i in 1:n_pairs) {
      pair_name <- values$pairwise_results$pair_names[i]
      pair_result <- values$pairwise_results$results[[pair_name]]
      
      if(any(pair_result$sig_regions)) {
        sig_starts <- which(diff(c(0, pair_result$sig_regions)) == 1)
        sig_ends <- which(diff(c(pair_result$sig_regions, 0)) == -1)
        
        for(j in 1:length(sig_starts)) {
          p <- p %>% add_trace(
            x = time_points[sig_starts[j]:sig_ends[j]],
            y = rep(i, sig_ends[j] - sig_starts[j] + 1),
            line = list(color = colors[i], width = 8),
            showlegend = (j == 1),
            name = pair_name,
            legendgroup = pair_name,
            hovertemplate = paste(pair_name,
                                  "<br>Time: %{customdata}",
                                  "<br>Significant<extra></extra>"),
            customdata = hover_times[sig_starts[j]:sig_ends[j]]
          )
        }
      } else {
        p <- p %>% add_trace(
          x = c(NA),
          y = c(NA),
          name = paste(pair_name, "(n.s.)"),
          line = list(color = 'gray'),
          showlegend = TRUE
        )
      }
    }
    
    for(t in seq(0, 1, by = 0.25)) {
      p <- p %>% add_trace(
        x = c(t, t),
        y = c(0.5, n_pairs + 0.5),
        line = list(color = 'lightgray', width = 0.5, dash = 'dot'),
        showlegend = FALSE,
        hoverinfo = 'skip'
      )
    }
    
    p <- p %>% layout(
      title = "Timeline of Significant Regions",
      xaxis = list(title = "Time", range = c(0, 1)),
      yaxis = list(
        title = "Comparison",
        tickmode = "array",
        tickvals = 1:n_pairs,
        ticktext = values$pairwise_results$pair_names,
        range = c(0.5, n_pairs + 0.5),
        autorange = 'reversed'
      ),
      hovermode = 'x',
      plot_bgcolor = 'rgba(240, 240, 240, 0.5)',
      legend = list(x = 1.02, y = 1, font = list(size = 10))
    )
    
    # Apply time label formatting
    p <- format_plotly_time_axis(p, tick_step_hours = as.numeric(input$tick_freq_pairwise))
    p
  })
  
  output$pairwise_regions_table <- renderDT({
    req(values$pairwise_results)
    
    # Placeholder table - would need implementation
    datatable(data.frame(Message = "Regions table not yet implemented"),
              options = list(dom = 't'),
              rownames = FALSE)
  })
