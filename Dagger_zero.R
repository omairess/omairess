# Enhanced Bayesian Network DAG Analysis Shiny App with Multi-Source and Multi-Target Constraint Support
# Modified to allow multiple selections in both "FROM:" and "TO:" dropdowns for blacklist and whitelist

# Load required packages
packages <- c("shiny", "bnlearn", "qgraph", "RColorBrewer", "sortable",
              "DT", "shinycssloaders", "shinyWidgets", "plotly", "igraph",
              "visNetwork", "shinyjs", "jsonlite", "colourpicker",
              "readxl", "haven")
new_pkgs <- setdiff(packages, rownames(installed.packages()))
if(length(new_pkgs)) install.packages(new_pkgs)

library(shiny)
library(bnlearn)
library(qgraph)
library(RColorBrewer)
library(sortable)
library(DT)
library(shinycssloaders)
library(shinyWidgets)
library(visNetwork)
library(shinyjs)
library(jsonlite)
library(colourpicker)
library(readxl)
library(haven)

# ============================================================================
# CONSTANTS AND CONFIGURATION
# ============================================================================

ALGORITHMS <- list(
  "Hill-Climbing" = "hc",
  "Tabu Search" = "tabu",
  "Max-Min Hill Climbing" = "mmhc",
  "PC Stable" = "pc.stable",
  "Grow-Shrink" = "gs",
  "Incremental Association" = "iamb",
  "RSMAX2" = "rsmax2"
)

ALGORITHM_INFO <- list(
  "hc" = "Hill-climbing: greedy search with random restarts",
  "tabu" = "Tabu search: explores larger neighborhoods, avoids cycles",
  "mmhc" = "Max-Min Hill Climbing: hybrid constraint + score-based approach",
  "pc.stable" = "PC Stable: modern constraint-based algorithm using conditional independence tests. Order-independent and more stable than original PC. Gold standard for causal discovery.",
  "gs" = "Grow-Shrink: constraint-based Markov blanket discovery",
  "iamb" = "Incremental Association: constraint-based causal discovery",
  "rsmax2" = "RSMAX2: constraint-based with statistical tests"
)

SCORE_TYPES <- list(
  # === DISCRETE DATA SCORES ===
  "BIC (Discrete)" = "bic",
  "AIC (Discrete)" = "aic",
  "Extended BIC (Discrete)" = "ebic",
  "Log-Likelihood (Discrete)" = "loglik",
  "Predictive Log-Likelihood (Discrete)" = "pred-loglik",
  "BDe - Bayesian Dirichlet Equivalent" = "bde",
  "BDs - Bayesian Dirichlet Sparse" = "bds",
  "BDj - Bayesian Dirichlet Jeffrey's" = "bdj",
  "K2 - Cooper & Herskovits" = "k2",
  "MBDE - Mixed Experimental Data" = "mbde",
  "BDLA - Locally Averaged Bayesian" = "bdla",
  "FNML - Factorized Normalized ML" = "fnml",
  "QNML - Quotient Normalized ML" = "qnml",
  "NAL - Node-Average Likelihood" = "nal",
  "PNAL - Penalized Node-Average" = "pnal",

  # === GAUSSIAN DATA SCORES ===
  "BIC (Gaussian)" = "bic-g",
  "AIC (Gaussian)" = "aic-g",
  "Extended BIC (Gaussian)" = "ebic-g",
  "Log-Likelihood (Gaussian)" = "loglik-g",
  "Predictive Log-Likelihood (Gaussian)" = "pred-loglik-g",
  "BGE - Bayesian Gaussian Equivalent" = "bge",
  "NAL (Gaussian)" = "nal-g",
  "PNAL (Gaussian)" = "pnal-g",

  # === CONDITIONAL GAUSSIAN (MIXED) SCORES ===
  "BIC (Conditional Gaussian)" = "bic-cg",
  "AIC (Conditional Gaussian)" = "aic-cg",
  "Extended BIC (Conditional Gaussian)" = "ebic-cg",
  "Log-Likelihood (Conditional Gaussian)" = "loglik-cg",
  "Predictive Log-Likelihood (Conditional Gaussian)" = "pred-loglik-cg",
  "NAL (Conditional Gaussian)" = "nal-cg",
  "PNAL (Conditional Gaussian)" = "pnal-cg"
)

SCORE_INFO <- list(
  # === DISCRETE DATA SCORES ===
  "bic" = "BIC (Discrete): Bayesian Information Criterion for discrete data. Penalizes complexity heavily, good for larger datasets. Prefers simpler models.",
  "aic" = "AIC (Discrete): Akaike Information Criterion for discrete data. More forgiving than BIC, better for smaller datasets. Allows slightly more complex models.",
  "ebic" = "Extended BIC (Discrete): Adds extra penalty to BIC for dense networks. Best when you want very sparse networks with few edges.",
  "loglik" = "Log-Likelihood (Discrete): Raw goodness-of-fit measure. No complexity penalty. Equivalent to entropy measure used in Weka.",
  "pred-loglik" = "Predictive Log-Likelihood (Discrete): Evaluated on separate test set. Good for assessing generalization. Requires test data.",
  "bde" = "BDe - Bayesian Dirichlet Equivalent: Bayesian score with uniform prior. Score equivalent (arc direction doesn't affect score). Good with domain knowledge.",
  "bds" = "BDs - Bayesian Dirichlet Sparse: Sparsity-inducing Bayesian score. Encourages fewer edges. Good for discovering sparse networks.",
  "bdj" = "BDj - Bayesian Dirichlet Jeffrey's: Uses Jeffrey's non-informative prior. Less sensitive to prior assumptions than BDe.",
  "k2" = "K2 - Cooper & Herskovits: NOT score equivalent - arc direction matters. Good when causal direction is important to preserve.",
  "mbde" = "MBDE - Mixed Experimental Data: For datasets mixing experimental and observational data. Handles interventional studies.",
  "bdla" = "BDLA - Locally Averaged Bayesian: Averages over local structures. Good for datasets with local dependencies.",
  "fnml" = "FNML - Factorized Normalized ML: Advanced normalized likelihood approach. Good theoretical properties.",
  "qnml" = "QNML - Quotient Normalized ML: Newer normalized ML variant. Often better than FNML in practice.",
  "nal" = "NAL - Node-Average Likelihood: Handles incomplete data well. Averages likelihood over nodes. Good with missing values.",
  "pnal" = "PNAL - Penalized Node-Average: NAL with complexity penalty. Better than NAL for model selection with incomplete data.",

  # === GAUSSIAN DATA SCORES ===
  "bic-g" = "BIC (Gaussian): BIC for continuous/Gaussian data. Standard choice for continuous networks. Penalizes complexity appropriately.",
  "aic-g" = "AIC (Gaussian): AIC for continuous data. More permissive than BIC-G. Good for smaller continuous datasets.",
  "ebic-g" = "Extended BIC (Gaussian): Extra penalty for dense continuous networks. Use when you want sparse continuous networks.",
  "loglik-g" = "Log-Likelihood (Gaussian): Raw fit for continuous data. No complexity penalty. Shows pure goodness-of-fit.",
  "pred-loglik-g" = "Predictive Log-Likelihood (Gaussian): Evaluated on test set. Good for assessing continuous model generalization.",
  "bge" = "BGE - Bayesian Gaussian Equivalent: Bayesian score for continuous data. Score equivalent. Good with prior knowledge about continuous relationships.",
  "nal-g" = "NAL (Gaussian): Node-average likelihood for continuous data with missing values. Handles incomplete continuous data.",
  "pnal-g" = "PNAL (Gaussian): Penalized NAL for continuous data. Better model selection with incomplete continuous data.",

  # === CONDITIONAL GAUSSIAN (MIXED) SCORES ===
  "bic-cg" = "BIC (Conditional Gaussian): For mixed discrete + continuous data. Standard choice for hybrid networks.",
  "aic-cg" = "AIC (Conditional Gaussian): More permissive mixed-data score. Good for smaller mixed datasets.",
  "ebic-cg" = "Extended BIC (Conditional Gaussian): Sparse mixed networks. Use when you want few edges in mixed-type networks.",
  "loglik-cg" = "Log-Likelihood (Conditional Gaussian): Raw fit for mixed data. Shows pure goodness-of-fit without penalty.",
  "pred-loglik-cg" = "Predictive Log-Likelihood (Conditional Gaussian): Test set evaluation for mixed data. Assesses mixed model generalization.",
  "nal-cg" = "NAL (Conditional Gaussian): Handles missing values in mixed discrete/continuous data very well.",
  "pnal-cg" = "PNAL (Conditional Gaussian): Penalized version for model selection with incomplete mixed data."
)

LAYOUTS <- list("Spring" = "spring", "Circle" = "circle", "Cascade" = "cascade")

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

# Auto-detect data types based on strict heuristics
# Only classify variables as discrete if they have exactly 2 unique values

auto_detect_data_types <- function(df, selected_vars) {
  results <- list(continuous = character(0), discrete = character(0))

  for (var in selected_vars) {
    col_data <- df[[var]]
    col_data <- col_data[!is.na(col_data)]

    if (length(col_data) == 0) next

    unique_vals <- length(unique(col_data))

    # Discrete if factor/character with 2 unique values, or numeric with exactly 2 unique values
    if ((is.factor(col_data) || is.character(col_data)) && unique_vals == 2) {
      results$discrete <- c(results$discrete, var)
    } else if (is.numeric(col_data) && unique_vals == 2) {
      results$discrete <- c(results$discrete, var)
    } else {
      results$continuous <- c(results$continuous, var)
    }
  }

  return(results)
}

# Read a data file based on its extension
read_data_file <- function(filepath, filename) {
  ext <- tolower(tools::file_ext(filename))
  df <- switch(ext,
    "csv"  = read.csv(filepath, header = TRUE, stringsAsFactors = FALSE),
    "txt"  = read.table(filepath, header = TRUE, stringsAsFactors = FALSE, sep = "\t", fill = TRUE),
    "tsv"  = read.delim(filepath, header = TRUE, stringsAsFactors = FALSE),
    "xlsx" = as.data.frame(readxl::read_excel(filepath)),
    "xls"  = as.data.frame(readxl::read_excel(filepath)),
    "sav"  = as.data.frame(haven::read_sav(filepath)),
    stop(paste0("Unsupported file format: .", ext,
                ". Supported formats: .csv, .txt, .tsv, .xlsx, .xls, .sav"))
  )
  # Convert haven labelled columns to plain R types
  if (any(sapply(df, inherits, "haven_labelled"))) {
    df <- as.data.frame(haven::as_factor(df))
  }
  df
}

# Reshape long-by-time data to wide format.
# Each subject gets one row; measurement variables renamed var_<time>.
long_to_wide <- function(df, id_col, time_col) {
  vars      <- setdiff(names(df), c(id_col, time_col))
  if (length(vars) == 0) stop("No measurement columns found after removing id and time columns.")
  time_vals <- sort(unique(df[[time_col]]))
  if (length(time_vals) < 2) stop("Time column must have at least 2 distinct values.")

  subsets <- lapply(time_vals, function(t) {
    sub <- df[df[[time_col]] == t, c(id_col, vars), drop = FALSE]
    names(sub)[names(sub) != id_col] <- paste0(vars, "_", t)
    sub
  })

  result <- Reduce(function(a, b) merge(a, b, by = id_col, all = TRUE), subsets)
  result[[id_col]] <- NULL
  result
}

# Validate uploaded data
validate_data <- function(df) {
  errors <- character(0)

  if (nrow(df) < 10) {
    errors <- c(errors, "Dataset too small. Need at least 10 observations.")
  }

  if (ncol(df) < 2) {
    errors <- c(errors, "Dataset needs at least 2 variables.")
  }

  # Check for at least 2 analyzable variables (numeric, factor, or character)
  analyzable_vars <- sapply(df, function(x)
    is.numeric(x) || is.factor(x) || is.character(x) || all(is.na(x) | suppressWarnings(!is.na(as.numeric(x)))))
  if (sum(analyzable_vars) < 2) {
    errors <- c(errors, "Need at least 2 analyzable variables for DAG analysis.")
  }

  return(errors)
}

# Get analyzable variable names (numeric, factor, character)
get_analyzable_vars <- function(df) {
  names(df)[sapply(df, function(x)
    is.numeric(x) || is.factor(x) || is.character(x) || all(is.na(x) | suppressWarnings(!is.na(as.numeric(x)))))]
}

# Prepare mixed data for analysis
prepare_mixed_data <- function(df, continuous_vars, discrete_vars, omit_na = TRUE) {
  # Store original data before any processing
  original_df <- df

  all_vars <- c(continuous_vars, discrete_vars)
  if (length(all_vars) == 0) {
    stop("No variables selected")
  }

  # First subset to only the variables we need
  data <- df[, all_vars, drop = FALSE]

  # THEN remove NA if requested (only from the selected variables)
  if (omit_na) {
    data <- na.omit(data)
  }

  # Handle continuous variables
  if (length(continuous_vars) > 0) {
    cont_data <- as.matrix(data[, continuous_vars, drop = FALSE])
    storage.mode(cont_data) <- "numeric"

    # Check for zero variance (only for non-missing values)
    if (omit_na) {
      var_check <- apply(cont_data, 2, var, na.rm = TRUE)
      zero_vars <- names(var_check)[var_check == 0 | is.na(var_check)]

      if (length(zero_vars) > 0) {
        stop(paste("Continuous variables with zero variance:", paste(zero_vars, collapse = ", ")))
      }
    }

    data[, continuous_vars] <- as.data.frame(cont_data)
  }

  # Handle discrete variables
  if (length(discrete_vars) > 0) {
    for (var in discrete_vars) {
      if (!is.factor(data[[var]])) {
        data[[var]] <- as.factor(data[[var]])
      }
      # Check levels only for non-missing
      if (length(levels(data[[var]])) < 2) {
        stop(paste("Discrete variable", var, "must have at least 2 levels"))
      }
    }
  }

  # Only check row count if we removed NAs
  if (omit_na && nrow(data) < 5) {
    stop("Not enough observations after removing missing data.")
  }

  # Store original data as attribute for later reference
  attr(data, "original") <- original_df[, all_vars, drop = FALSE]

  return(data)
}

# FIX: Calculate common layout for all groups in split analysis
calculate_common_layout <- function(all_vars, split_results, layout_type = "spring") {
  n <- length(all_vars)

  if (layout_type == "circle") {
    # Circle layout - simple and consistent
    angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
    coords <- matrix(c(cos(angles), sin(angles)), ncol = 2)
    rownames(coords) <- all_vars
    return(coords)
  } else if (layout_type == "cascade") {
    # Combine all edges from all groups for cascade calculation
    all_edges_df <- do.call(rbind, lapply(split_results, function(res) {
      if (nrow(res$boot_selected) > 0) {
        data.frame(from = res$boot_selected$from,
                  to = res$boot_selected$to,
                  weight = res$boot_selected$strength,
                  stringsAsFactors = FALSE)
      } else {
        NULL
      }
    }))

    if (!is.null(all_edges_df) && nrow(all_edges_df) > 0 && requireNamespace("igraph", quietly = TRUE)) {
      # Average weights for duplicate edges
      edge_summary <- aggregate(weight ~ from + to, data = all_edges_df, FUN = mean)

      # Create adjacency matrix
      adj <- matrix(0, n, n, dimnames = list(all_vars, all_vars))
      for(i in seq_len(nrow(edge_summary))) {
        from_idx <- which(all_vars == edge_summary$from[i])
        to_idx <- which(all_vars == edge_summary$to[i])
        if(length(from_idx) > 0 && length(to_idx) > 0) {
          adj[from_idx, to_idx] <- edge_summary$weight[i]
        }
      }

      # Calculate Katz centrality
      g <- igraph::graph_from_adjacency_matrix(adj, mode = "directed", weighted = TRUE)

      # Handle disconnected nodes
      tryCatch({
        katz_vals <- igraph::alpha_centrality(g, alpha = 0.1)
      }, error = function(e) {
        # If Katz centrality fails, use degree centrality
        katz_vals <- igraph::degree(g, mode = "in")
      })

      if(length(unique(katz_vals)) > 1) {
        # Create hierarchical layout
        y_coords <- 1 - ((katz_vals - min(katz_vals)) / (max(katz_vals) - min(katz_vals)))
      } else {
        # All nodes have same centrality, distribute evenly
        y_coords <- rep(0.5, n)
      }

      coords <- matrix(nrow = n, ncol = 2)

      # Group by levels and distribute horizontally
      # Use 8 discrete levels for finer granularity
      y_levels <- round(y_coords * 8) / 8
      for (level in unique(y_levels)) {
        nodes_in_level <- which(abs(y_levels - level) < 0.07)
        if (length(nodes_in_level) > 1) {
          x_positions <- seq(-1, 1, length.out = length(nodes_in_level))
          coords[nodes_in_level, 1] <- x_positions
        } else {
          coords[nodes_in_level, 1] <- 0
        }
        # Snap all nodes in the same level to exactly the same y coordinate
        # to prevent within-level vertical overlap
        coords[nodes_in_level, 2] <- level
      }

      rownames(coords) <- all_vars
      return(coords)
    } else {
      # Fallback to circle if cascade fails
      angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
      coords <- matrix(c(cos(angles), sin(angles)), ncol = 2)
      rownames(coords) <- all_vars
      return(coords)
    }
  }

  # Default: Spring layout - use qgraph to calculate once
  # Create combined adjacency matrix from all groups
  adj <- matrix(0, n, n, dimnames = list(all_vars, all_vars))

  for (result in split_results) {
    if (nrow(result$boot_selected) > 0) {
      for(i in seq_len(nrow(result$boot_selected))) {
        from <- result$boot_selected$from[i]
        to <- result$boot_selected$to[i]
        w <- result$boot_selected$strength[i]
        from_idx <- which(all_vars == from)
        to_idx <- which(all_vars == to)
        if(length(from_idx) > 0 && length(to_idx) > 0) {
          adj[from_idx, to_idx] <- max(adj[from_idx, to_idx], w)  # Use max strength across groups
        }
      }
    }
  }

  # Calculate layout using qgraph
  set.seed(123)

  # Check if there are any edges
  if(sum(adj) > 0) {
    layout_result <- qgraph::qgraph(adj, DoNotPlot = TRUE, layout = "spring")$layout
  } else {
    # No edges, use circle layout as fallback
    angles <- seq(0, 2*pi, length.out = n + 1)[1:n]
    layout_result <- matrix(c(cos(angles), sin(angles)), ncol = 2)
  }

  rownames(layout_result) <- all_vars

  return(layout_result)
}

# Create DAG difference network plot comparing two groups
create_dag_diff_plot <- function(result1, result2, group1_name, group2_name, params, common_layout = NULL) {
  tryCatch({
    vars1 <- c(result1$stats$continuous_vars, result1$stats$discrete_vars)
    vars2 <- c(result2$stats$continuous_vars, result2$stats$discrete_vars)
    vars  <- union(vars1, vars2)
    n     <- length(vars)

    bs1 <- result1$boot_selected
    bs2 <- result2$boot_selected

    key <- function(df) if (nrow(df) > 0) paste(df$from, df$to, sep = "->") else character(0)
    k1 <- key(bs1); k2 <- key(bs2)

    shared_k <- intersect(k1, k2)
    only1_k  <- setdiff(k1, k2)
    only2_k  <- setdiff(k2, k1)

    make_edges <- function(df, keys, color, group_label) {
      if (length(keys) == 0) return(NULL)
      rows <- df[key(df) %in% keys, , drop = FALSE]
      data.frame(
        from  = rows$from,
        to    = rows$to,
        color = color,
        width = pmax(1, as.numeric(rows$strength) * (params$edge_width %||% 0.5) * 5),
        title = paste0(group_label, "<br>Strength: ", round(as.numeric(rows$strength), 3)),
        label = if (!is.null(params$show_labels) && params$show_labels)
                  as.character(round(as.numeric(rows$strength), 3)) else "",
        stringsAsFactors = FALSE
      )
    }

    show_shared <- isTRUE(params$diff_show_shared %||% TRUE)
    show_only1  <- isTRUE(params$diff_show_only1  %||% TRUE)
    show_only2  <- isTRUE(params$diff_show_only2  %||% TRUE)

    edges <- rbind(
      if (show_shared) make_edges(bs1, shared_k, "#5C8DB8", "Both groups")  else NULL,
      if (show_only1)  make_edges(bs1, only1_k,  "#E53935", paste0("Only: ", group1_name)) else NULL,
      if (show_only2)  make_edges(bs2, only2_k,  "#43A047", paste0("Only: ", group2_name)) else NULL
    )
    if (is.null(edges)) edges <- data.frame(from=character(0), to=character(0),
                                            color=character(0), width=numeric(0),
                                            title=character(0), label=character(0),
                                            stringsAsFactors=FALSE)

    # Curve bidirectional pairs so they don't overlap
    if (nrow(edges) > 0) {
      fwd     <- paste(edges$from, edges$to, sep = "->")
      rev_key <- paste(edges$to,   edges$from, sep = "->")
      bidir   <- fwd %in% rev_key
      if (any(bidir)) {
        edges$smooth.enabled   <- bidir
        edges$smooth.type      <- ifelse(bidir, "curvedCW", "continuous")
        edges$smooth.roundness <- ifelse(bidir, 0.25, 0)
      }
    }

    continuous_vars <- union(result1$stats$continuous_vars, result2$stats$continuous_vars)
    discrete_vars   <- setdiff(vars, continuous_vars)

    nodes <- data.frame(id = vars, label = vars, stringsAsFactors = FALSE)
    node_colors <- ifelse(vars %in% continuous_vars,
                          params$node_continuous_color %||% "#e8f4f8",
                          params$node_discrete_color   %||% "#fff3e0")
    nodes$color      <- node_colors
    nodes$size       <- (params$node_size %||% 10) * 3
    nodes$font.size  <- (params$label_size %||% 1.0) * 16
    nodes$font.color <- "black"
    nodes$font.strokeWidth <- if (!is.null(params$bold_labels) && params$bold_labels) 2 else 0
    nodes$font.strokeColor <- "white"
    nodes$borderWidth <- params$node_border_width %||% 1.5
    nodes$color.border <- params$node_border_color %||% "#000000"

    layout_type <- params$layout %||% "spring"

    if (!is.null(common_layout)) {
      node_positions <- common_layout[nodes$id, , drop = FALSE]
      nodes$x <- node_positions[, 1] * 300
      nodes$y <- node_positions[, 2] * 300
    }

    title_str <- paste0("DAG Difference: ", group1_name, " vs ", group2_name,
                        "<br><span style='color:#E53935'>&#9632;</span> Only in ", group1_name,
                        " &nbsp; <span style='color:#43A047'>&#9632;</span> Only in ", group2_name,
                        " &nbsp; <span style='color:#5C8DB8'>&#9632;</span> Shared")

    net <- visNetwork(nodes, edges, main = list(text = title_str, style = "font-size:13px;")) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
      visEdges(
        arrows = list(to = list(enabled = TRUE, scaleFactor = 0.8)),
        smooth = list(enabled = TRUE, type = "curvedCW",
                      roundness = if (layout_type == "cascade") 0.35 else 0.2)
      ) %>%
      visNodes(
        borderWidth = params$node_border_width %||% 1.5,
        color = list(highlight = list(border = "#2196f3", background = "#bbdefb"))
      )

    if (!is.null(common_layout)) {
      if (layout_type == "spring") {
        net <- net %>% visPhysics(
          enabled = TRUE,
          solver = "repulsion",
          repulsion = list(nodeDistance = 0),
          stabilization = list(enabled = TRUE, iterations = 0)
        )
      } else {
        net <- net %>% visPhysics(enabled = FALSE)
      }
    } else {
      if (layout_type == "circle") {
        net <- net %>% visIgraphLayout(layout = "layout_in_circle")
      } else if (layout_type == "cascade") {
        tryCatch(
          net <- net %>% visIgraphLayout(layout = "layout_with_sugiyama"),
          error = function(e) net <<- net %>% visLayout(randomSeed = 123)
        )
      } else {
        net <- net %>% visPhysics(
          solver = "forceAtlas2Based",
          forceAtlas2Based = list(gravitationalConstant = -50),
          stabilization = list(iterations = 100)
        ) %>% visLayout(randomSeed = 123)
      }
    }
    net
  }, error = function(e) {
    nodes <- data.frame(id = "error", label = paste("Error:", conditionMessage(e)), stringsAsFactors = FALSE)
    edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
    visNetwork(nodes, edges)
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Create interactive network plot with drag-and-drop positioning
create_interactive_network_plot <- function(boot_selected, continuous_vars, discrete_vars, data, params,
                                           title_suffix = "", common_layout = NULL, group_name = NULL,
                                           original_data = NULL) {

  tryCatch({
    vars <- c(continuous_vars, discrete_vars)
    n <- length(vars)

    # Filter edges by minimum size if specified
    if (!is.null(params$min_edge_size) && params$min_edge_size > 0) {
      boot_selected <- boot_selected[boot_selected$strength >= params$min_edge_size, ]
    }

    # Prepare nodes first
    nodes <- data.frame(
      id = vars,
      label = vars,
      stringsAsFactors = FALSE
    )

    # Calculate node means using colMeans with na.rm = TRUE
    node_means <- numeric(n)
    names(node_means) <- vars

    if (length(continuous_vars) > 0) {
      cont_means <- colMeans(data[continuous_vars], na.rm = TRUE)
      node_means[continuous_vars] <- cont_means
    }

    # Set discrete variables to neutral value
    if (length(discrete_vars) > 0) {
      if (length(continuous_vars) > 0) {
        node_means[discrete_vars] <- mean(node_means[continuous_vars], na.rm = TRUE)
      } else {
        node_means[discrete_vars] <- 0.5
      }
    }

    # Node styling with data type awareness
    if (params$scale_nodes && length(continuous_vars) > 0) {
      if(length(unique(node_means)) > 1) {
        normalized_means <- (node_means - min(node_means)) / (max(node_means) - min(node_means))
      } else {
        normalized_means <- rep(0.5, length(node_means))
      }
      node_sizes <- normalized_means * (params$max_node_size - params$min_node_size) + params$min_node_size

      # More subtle alpha variation for colors with data type differentiation
      alpha_vals <- normalized_means * 0.4 + 0.6
      node_colors <- character(n)
      # Convert hex to RGB with alpha for continuous variables
      cont_col_rgb <- col2rgb(params$node_continuous_color) / 255
      disc_col_rgb <- col2rgb(params$node_discrete_color) / 255
      node_colors[vars %in% continuous_vars] <- rgb(cont_col_rgb[1], cont_col_rgb[2], cont_col_rgb[3],
                                                     alpha = alpha_vals[vars %in% continuous_vars])
      node_colors[vars %in% discrete_vars] <- rgb(disc_col_rgb[1], disc_col_rgb[2], disc_col_rgb[3],
                                                   alpha = alpha_vals[vars %in% discrete_vars])

      nodes$size <- pmax(10, node_sizes * 3)  # Ensure minimum size and scale for visNetwork
      nodes$color <- node_colors

      # Add title attribute with mean values
      nodes$title <- paste0(nodes$id, "<br>Mean: ", round(node_means, 3))
      if (!is.null(group_name)) {
        nodes$title <- paste0(nodes$title, "<br>Group: ", group_name)
      }

    } else {
      nodes$size <- params$node_size * 3
      # Different colors for different data types using user-selected colors
      node_colors <- character(n)
      node_colors[vars %in% continuous_vars] <- params$node_continuous_color
      node_colors[vars %in% discrete_vars] <- params$node_discrete_color
      nodes$color <- node_colors

      # Add title for hover info
      nodes$title <- nodes$id
      if (!is.null(group_name)) {
        nodes$title <- paste0(nodes$title, "<br>Group: ", group_name)
      }
    }

    # FIX: Apply common layout if provided
    if (!is.null(common_layout)) {
      # Ensure nodes are in the same order as the common layout
      node_positions <- common_layout[nodes$id, , drop = FALSE]
      nodes$x <- node_positions[, 1] * 300  # Scale for visNetwork
      nodes$y <- node_positions[, 2] * 300
      # Note: physics settings are handled separately based on layout type
      # Circle/Cascade: physics disabled to prevent jiggling
      # Spring: minimal physics to allow natural movement
    }

    # Add font styling
    nodes$font.size <- params$label_size * 16
    nodes$font.color <- "black"
    nodes$font.strokeWidth <- if(params$bold_labels) 2 else 0
    nodes$font.strokeColor <- "white"

    # Handle empty edges case
    if (nrow(boot_selected) == 0) {
      edges <- data.frame(
        from = character(0),
        to = character(0),
        width = numeric(0),
        color = character(0),
        label = character(0),
        stringsAsFactors = FALSE
      )

      net <- visNetwork(nodes, edges, main = title_suffix) %>%
             visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
             visInteraction(dragNodes = TRUE,  # Always allow dragging
                           dragView = TRUE, zoomView = TRUE)

      if (is.null(common_layout)) {
        # Apply layout based on params when no common layout
        if (params$layout == "circle") {
          net <- net %>% visIgraphLayout(layout = "layout_in_circle")
        } else if (params$layout == "cascade") {
          net <- net %>% visIgraphLayout(layout = "layout_with_sugiyama")
        } else {
          net <- net %>% visLayout(randomSeed = 123)
        }
      } else {
        # Use appropriate physics based on layout type
        if (params$layout == "spring") {
          # Spring: minimal physics
          net <- net %>% visPhysics(
            enabled = TRUE,
            solver = "repulsion",
            repulsion = list(nodeDistance = 0),
            stabilization = list(enabled = TRUE, iterations = 0)
          )
        } else {
          # Circle/Cascade: check if user wants physics
          if (!is.null(params$enable_physics_for_fixed) && params$enable_physics_for_fixed) {
            net <- net %>% visPhysics(
              enabled = TRUE,
              solver = "repulsion",
              repulsion = list(nodeDistance = 0),
              stabilization = list(enabled = TRUE, iterations = 0)
            )
          } else {
            # No physics to prevent jiggling
            net <- net %>% visPhysics(enabled = FALSE)
          }
        }
      }

      return(net)
    }

    # Prepare edges - simplified approach
    edges <- data.frame(
      from = boot_selected$from,
      to = boot_selected$to,
      stringsAsFactors = FALSE
    )

    # Select which IC score column to use (absolute or relative)
    ic_use_relative <- !is.null(params$ic_score_type) && params$ic_score_type == "relative"
    ic_col <- if (ic_use_relative && "score_delta_relative" %in% names(boot_selected)) {
      "score_delta_relative"
    } else {
      "score_delta"
    }

    # Calculate edge widths based on edge_display_type
    edge_width_values <- if (!is.null(params$edge_display_type)) {
      switch(params$edge_display_type,
             "strength" = as.numeric(boot_selected$strength),
             "direction" = as.numeric(boot_selected$direction),
             "combined" = as.numeric(boot_selected$strength) * as.numeric(boot_selected$direction),
             "ic" = {
               if (ic_col %in% names(boot_selected) && sum(!is.na(boot_selected[[ic_col]])) > 0) {
                 abs_deltas <- abs(as.numeric(boot_selected[[ic_col]]))
                 abs_deltas_clean <- abs_deltas[!is.na(abs_deltas)]
                 if (length(abs_deltas_clean) > 0 && max(abs_deltas_clean) > min(abs_deltas_clean)) {
                   normalized <- (abs_deltas - min(abs_deltas_clean)) / (max(abs_deltas_clean) - min(abs_deltas_clean))
                   normalized[is.na(normalized)] <- 0
                   normalized
                 } else if (length(abs_deltas_clean) > 0) {
                   result <- rep(0.5, length(abs_deltas))
                   result[is.na(abs_deltas)] <- 0
                   result
                 } else {
                   as.numeric(boot_selected$strength)
                 }
               } else {
                 try(showNotification(
                   "BIC/AIC score data (score_delta) not available for this analysis. Edge width falling back to Bootstrap Strength.",
                   type = "warning", duration = 8, id = "ic_fallback_warn"
                 ), silent = TRUE)
                 as.numeric(boot_selected$strength)
               }
             },
             as.numeric(boot_selected$strength)  # default fallback
      )
    } else {
      as.numeric(boot_selected$strength)  # default if not specified
    }

    # Normalize to [0,1] across displayed edges so the full width range is used
    ev_min <- min(edge_width_values, na.rm = TRUE)
    ev_max <- max(edge_width_values, na.rm = TRUE)
    if (ev_max > ev_min) {
      edge_width_norm <- (edge_width_values - ev_min) / (ev_max - ev_min)
    } else {
      edge_width_norm <- rep(0.5, length(edge_width_values))
    }
    min_w <- 1
    max_w <- max(2, params$edge_width * 10)
    edges$width <- min_w + edge_width_norm * (max_w - min_w)


    # Calculate edge colors
    strength_vals <- as.numeric(boot_selected$strength)
    if(length(unique(strength_vals)) > 1) {
      pal <- brewer.pal(min(9, max(3, length(unique(strength_vals)))), params$palette)
      strength_breaks <- seq(min(strength_vals), max(strength_vals), length.out = length(pal))
      edge_colors <- pal[findInterval(strength_vals, strength_breaks, all.inside = TRUE)]
    } else {
      edge_colors <- rep(brewer.pal(3, params$palette)[2], length(strength_vals))
    }

    # Apply transparency to colors
    if (!is.null(params$edge_transparency) && params$edge_transparency < 1) {
      edge_colors_rgb <- col2rgb(edge_colors)
      edge_colors <- rgb(edge_colors_rgb[1,], edge_colors_rgb[2,], edge_colors_rgb[3,],
                        alpha = params$edge_transparency * 255, maxColorValue = 255)
    }

    edges$color <- edge_colors

    # Add edge labels if requested
    if (params$show_labels) {
      edges$label <- if (!is.null(params$edge_display_type)) {
        switch(params$edge_display_type,
          "strength" = round(as.numeric(boot_selected$strength), 3),
          "direction" = round(as.numeric(boot_selected$direction), 3),
          "combined" = round(as.numeric(boot_selected$strength) * as.numeric(boot_selected$direction), 3),
          "ic" = {
            raw <- if (ic_col %in% names(boot_selected)) as.numeric(boot_selected[[ic_col]]) else as.numeric(boot_selected$strength)
            ifelse(is.na(raw), NA_character_, as.character(round(raw, 3)))
          },
          round(as.numeric(boot_selected$strength), 3)  # default
        )
      } else {
        round(as.numeric(boot_selected$strength), 3)
      }
    } else {
      edges$label <- ""
    }

    # Add edge font settings for bold labels
    if(params$show_labels) {
      edges$font.size <- params$label_size * 14
      edges$font.color <- if (!is.null(params$edge_display_type)) switch(params$edge_display_type,
        "strength" = "darkblue", "direction" = "darkred", "combined" = "darkorchid", "ic" = "darkgreen", "darkblue"
      ) else "darkblue"
      edges$font.strokeWidth <- if(params$bold_edge_labels) 2 else 0
      edges$font.strokeColor <- "white"
      edges$font.bold <- params$bold_edge_labels
    }

    # Create the interactive network
    net <- visNetwork(nodes, edges, main = title_suffix) %>%
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE
      ) %>%
      visInteraction(
        dragNodes = TRUE,  # Always allow dragging, even with common layout
        dragView = TRUE,
        zoomView = TRUE,
        zoomSpeed = 0.2,
        selectConnectedEdges = FALSE
      ) %>%
      visEdges(
        arrows = list(to = list(enabled = TRUE, scaleFactor = max(0.5, params$arrow_size / 5))),
        smooth = list(enabled = TRUE, type = "curvedCW",
                      roundness = if (!is.null(params$layout) && params$layout == "cascade") 0.35 else 0.2)
      ) %>%
      visNodes(
        borderWidth = params$node_border_width,
        borderWidthSelected = params$node_border_width + 1,
        color = list(
          border = params$node_border_color,
          highlight = list(border = "#2196f3", background = "#bbdefb")
        )
      )

    # FIX: Don't apply random seed if using common layout
    if (is.null(common_layout)) {
      # Apply layout type when no common layout
      if (params$layout == "circle") {
        net <- net %>% visIgraphLayout(layout = "layout_in_circle")
      } else if (params$layout == "cascade") {
        if (requireNamespace("igraph", quietly = TRUE) && nrow(edges) > 0) {
          tryCatch({
            net <- net %>% visIgraphLayout(layout = "layout_with_sugiyama")
          }, error = function(e) {
            net <- net %>% visLayout(randomSeed = 123)
          })
        } else {
          net <- net %>% visLayout(randomSeed = 123)
        }
      } else {
        # Spring layout (default)
        net <- net %>% visPhysics(
          solver = "forceAtlas2Based",
          forceAtlas2Based = list(gravitationalConstant = -50),
          stabilization = list(iterations = 100)
        ) %>% visLayout(randomSeed = 123)
      }
    } else {
      # When using common layout, check layout type for physics settings
      if (params$layout == "spring") {
        # Spring layout: use minimal physics to allow dragging but maintain positions
        net <- net %>% visPhysics(
          enabled = TRUE,
          solver = "repulsion",
          repulsion = list(nodeDistance = 0),  # No repulsion to maintain positions
          stabilization = list(
            enabled = TRUE,
            iterations = 0  # Don't move nodes from their initial positions
          )
        )
      } else {
        # Circle or Cascade: check if user wants physics enabled
        if (!is.null(params$enable_physics_for_fixed) && params$enable_physics_for_fixed) {
          # User wants physics even for fixed layouts
          net <- net %>% visPhysics(
            enabled = TRUE,
            solver = "repulsion",
            repulsion = list(nodeDistance = 0),
            stabilization = list(
              enabled = TRUE,
              iterations = 0
            )
          )
        } else {
          # Default: disable physics to prevent jiggling
          net <- net %>% visPhysics(enabled = FALSE)
        }
      }
    }

    return(net)

  }, error = function(e) {
    # Fallback to simple network if there's an error
    vars <- c(continuous_vars, discrete_vars)
    nodes <- data.frame(
      id = vars,
      label = vars,
      size = 20,
      color = params$node_continuous_color,
      stringsAsFactors = FALSE
    )

    edges <- data.frame(
      from = character(0),
      to = character(0),
      stringsAsFactors = FALSE
    )

    return(visNetwork(nodes, edges, main = title_suffix) %>%
           visOptions(highlightNearest = TRUE) %>%
           visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
           visLayout(randomSeed = 123))
  })
}

# Compute per-edge curvature for cascade layout to avoid edges passing through nodes.
# Returns a matrix (same dims as adj) with curvature values: 0 = straight, +/-0.4 = curved.
compute_cascade_curves <- function(adj, coords) {
  n <- nrow(adj)
  curve_mat <- matrix(0, n, n)

  # Estimate node radius: coords are in [-1,1]; use 10% of range as threshold
  coord_range <- max(diff(range(coords[, 1])), diff(range(coords[, 2])), 0.5)
  node_radius <- coord_range * 0.10

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (adj[i, j] <= 0) next

      x1 <- coords[i, 1]; y1 <- coords[i, 2]
      x2 <- coords[j, 1]; y2 <- coords[j, 2]
      dx <- x2 - x1;      dy <- y2 - y1
      len_sq <- dx^2 + dy^2
      if (len_sq < 1e-10) next

      y_lo <- min(y1, y2); y_hi <- max(y1, y2)

      for (k in seq_len(n)) {
        if (k == i || k == j) next
        xk <- coords[k, 1]; yk <- coords[k, 2]

        # Only check nodes that lie in the y-range of this edge
        if (yk < y_lo - node_radius || yk > y_hi + node_radius) next

        # Distance from node k to the infinite line through i->j
        t <- max(0, min(1, ((xk - x1) * dx + (yk - y1) * dy) / len_sq))
        dist <- sqrt((xk - x1 - t * dx)^2 + (yk - y1 - t * dy)^2)

        if (dist < node_radius) {
          # Edge passes through node k - determine which side k is on
          # Cross product (i->j) x (i->k): positive => k is to the left of the edge
          cross <- dx * (yk - y1) - dy * (xk - x1)
          # Curve away from k: if k is to the left, bend right (positive), else left (negative)
          curve_mat[i, j] <- if (cross > 0) 0.4 else -0.4
          break
        }
      }
    }
  }

  curve_mat
}


# Create network plot
create_network_plot <- function(boot_selected, continuous_vars, discrete_vars, data, params,
                               title_suffix = "", common_layout = NULL, group_name = NULL) {
  vars <- c(continuous_vars, discrete_vars)
  n <- length(vars)
  adj <- matrix(0, n, n, dimnames = list(vars, vars))
  strengths <- numeric()
  directions <- numeric()

  # Filter edges by minimum size if specified
  if (!is.null(params$min_edge_size) && params$min_edge_size > 0) {
    boot_selected <- boot_selected[boot_selected$strength >= params$min_edge_size, ]
  }

  # Build adjacency matrix
  if (nrow(boot_selected) > 0) {
    for(i in seq_len(nrow(boot_selected))) {
      adj[boot_selected$from[i], boot_selected$to[i]] <- boot_selected$strength[i]
    }
  }

  # Rebuild edge vectors in qgraph's traversal order.
  # qgraph reads the adj matrix via which(adj!=0, arr.ind=TRUE), which returns
  # indices in R's native COLUMN-MAJOR order (ci outer, ri inner).
  score_deltas <- numeric()
  if (nrow(boot_selected) > 0) {
    ordered_rows <- integer(0)
    for (ci in seq_len(n)) {
      for (ri in seq_len(n)) {
        if (adj[ri, ci] > 0) {
          hit <- which(boot_selected$from == vars[ri] & boot_selected$to == vars[ci])
          if (length(hit) > 0) ordered_rows <- c(ordered_rows, hit[1])
        }
      }
    }
    boot_selected <- boot_selected[ordered_rows, , drop = FALSE]
    strengths    <- as.numeric(boot_selected$strength)
    directions   <- as.numeric(boot_selected$direction)
    use_relative <- !is.null(params$ic_score_type) && params$ic_score_type == "relative"
    if (use_relative && "score_delta_relative" %in% names(boot_selected)) {
      score_deltas <- boot_selected$score_delta_relative
    } else if ("score_delta" %in% names(boot_selected)) {
      score_deltas <- boot_selected$score_delta
    }
  }

  # Calculate edge width values based on edge_display_type
  edge_width_values <- if (!is.null(params$edge_display_type)) {
    switch(params$edge_display_type,
           "strength" = strengths,
           "direction" = directions,
           "combined" = strengths * directions,
           "ic" = {
             if (length(score_deltas) > 0 && sum(!is.na(score_deltas)) > 0) {
               # Use absolute values since score deltas can be negative
               # Normalize to 0-1 range for better visualization
               abs_deltas <- abs(score_deltas)
               # Remove NAs for min/max calculation
               abs_deltas_clean <- abs_deltas[!is.na(abs_deltas)]

               if (length(abs_deltas_clean) > 0 && max(abs_deltas_clean) > min(abs_deltas_clean)) {
                 # Normalize, filling NAs with 0 (minimum)
                 normalized <- (abs_deltas - min(abs_deltas_clean)) / (max(abs_deltas_clean) - min(abs_deltas_clean))
                 normalized[is.na(normalized)] <- 0
                 normalized
               } else if (length(abs_deltas_clean) > 0) {
                 # All same value, use middle
                 result <- rep(0.5, length(abs_deltas))
                 result[is.na(abs_deltas)] <- 0
                 result
               } else {
                 # All NA, fall back to strength
                 strengths
               }
             } else {
               try(showNotification(
                 "BIC/AIC score data (score_delta) not available for this analysis. Edge width falling back to Bootstrap Strength.",
                 type = "warning", duration = 8, id = "ic_fallback_warn"
               ), silent = TRUE)
               strengths
             }
           },
           strengths  # default fallback
    )
  } else {
    strengths  # default if not specified
  }

  if (length(strengths) == 0) {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, paste("No edges meet the threshold criteria", title_suffix, "\nTry lowering the thresholds or minimum edge size."),
         cex = 1.5, col = "red")
    return()
  }

  # FIX: Use common layout if provided, otherwise calculate
  if (!is.null(common_layout)) {
    layout_param <- common_layout[vars, , drop = FALSE]
  } else {
    # Handle layout types directly
    if (params$layout == "circle") {
      layout_param <- "circle"
    } else if (params$layout == "cascade") {
      tryCatch({
        if (!requireNamespace("igraph", quietly = TRUE)) {
          warning("igraph package required for cascade layout. Falling back to spring layout.")
          layout_param <- "spring"
        } else if (nrow(boot_selected) > 0) {
          # Build adjacency matrix for cascade calculation
          adj_binary <- (adj > 0) * 1
          g_DAG <- igraph::graph_from_adjacency_matrix(adj_binary, mode = "directed", diag = FALSE)

          # Calculate centrality
          tryCatch({
            katz_vals <- igraph::alpha_centrality(g_DAG, alpha = 0.1)
          }, error = function(e) {
            # Fallback to degree centrality
            katz_vals <- igraph::degree(g_DAG, mode = "in")
          })

          if(length(unique(katz_vals)) > 1) {
            katz_normalized <- (katz_vals - min(katz_vals)) / (max(katz_vals) - min(katz_vals))
          } else {
            katz_normalized <- rep(0.5, n)
          }

          custom_coords <- matrix(nrow = n, ncol = 2)
          y_coords <- 1 - katz_normalized
          # Use 8 discrete levels for finer granularity
          y_levels <- round(y_coords * 8) / 8

          for (level in unique(y_levels)) {
            nodes_in_level <- which(abs(y_levels - level) < 0.07)
            if (length(nodes_in_level) > 1) {
              x_positions <- seq(-1, 1, length.out = length(nodes_in_level))
              custom_coords[nodes_in_level, 1] <- x_positions
            } else {
              custom_coords[nodes_in_level, 1] <- 0
            }
            # Snap all nodes in the same level to exactly the same y coordinate
            custom_coords[nodes_in_level, 2] <- level
          }

          layout_param <- custom_coords
        } else {
          # No edges, use circle layout
          layout_param <- "circle"
        }
      }, error = function(e) {
        warning(paste("Error calculating cascade layout:", e$message, "Falling back to spring layout."))
        layout_param <- "spring"
      })
    } else {
      # Spring layout (default)
      layout_param <- "spring"
    }
  }

  # Calculate colors and node properties
  plot_data <- calculate_plot_properties(boot_selected, vars, continuous_vars, discrete_vars,
                                        data, strengths, params, group_name)

  # Apply edge transparency
  if (!is.null(params$edge_transparency)) {
    if (is.matrix(plot_data$edge_colors)) {
      plot_data$edge_colors <- apply(plot_data$edge_colors, c(1,2), function(x) {
        if (x != "#6b8eb1") {
          col_rgb <- col2rgb(x)
          rgb(col_rgb[1], col_rgb[2], col_rgb[3],
              alpha = params$edge_transparency * 255, maxColorValue = 255)
        } else {
          x
        }
      })
    } else {
      col_rgb <- col2rgb(plot_data$edge_colors)
      plot_data$edge_colors <- rgb(col_rgb[1,], col_rgb[2,], col_rgb[3,],
                                  alpha = params$edge_transparency * 255, maxColorValue = 255)
    }
  }

  edge_label_values <- if (params$show_labels) {
    # Follow the edge_display_type selection for labels
    if (!is.null(params$edge_display_type)) {
      switch(params$edge_display_type,
             "strength" = strengths,
             "direction" = directions,
             "combined" = strengths * directions,
             "ic" = {
               if (length(score_deltas) > 0 && sum(!is.na(score_deltas)) > 0) {
                 score_deltas  # keep NAs as NA — qgraph shows blank for NA labels
               } else {
                 strengths
               }
             },
             strengths  # default
      )
    } else {
      strengths  # default if not specified
    }
  } else {
    strengths  # fallback
  }

  # Define the title only when group_name is a string (split analysis)
  plot_title <- if (is.character(group_name)) paste("Group:", group_name) else NULL

  # For cascade layout, compute a per-edge curve matrix so edges bend around
  # intermediate nodes instead of passing through them.
  is_cascade <- !is.null(params$layout) && params$layout == "cascade"
  actual_layout_coords <- if (!is.null(common_layout)) common_layout[vars, , drop = FALSE] else layout_param
  cascade_curve_mat <- if (is_cascade && is.matrix(actual_layout_coords)) {
    compute_cascade_curves(adj, actual_layout_coords)
  } else {
    NULL
  }

  # Create the plot: include title.cex only when plot_title is not NULL
  qgraph(
    adj,
    directed = TRUE,
    layout = actual_layout_coords,
    edge.color = plot_data$edge_colors,
    vsize = plot_data$node_sizes,
    color = plot_data$node_colors,
    edge.width = {
      ev_min <- min(edge_width_values, na.rm = TRUE)
      ev_max <- max(edge_width_values, na.rm = TRUE)
      if (ev_max > ev_min) {
        norm_w <- (edge_width_values - ev_min) / (ev_max - ev_min)
      } else {
        norm_w <- rep(0.5, length(edge_width_values))
      }
      pmax(0.5, 0.5 + norm_w * params$edge_width * 4)
    },
    labels = plot_data$labels,
    label.cex = params$label_size,
    label.scale = FALSE,
    label.font = if(params$bold_labels) 2 else 1,
    edge.labels = if(params$show_labels) round(edge_label_values, 3) else FALSE,
    edge.label.cex = if(!is.null(params$edge_label_size)) params$edge_label_size else 0.7,
    edge.label.color = if (!is.null(params$edge_display_type)) {
      switch(params$edge_display_type,
             "strength" = "darkblue",
             "direction" = "darkred",
             "combined" = "darkorchid",
             "ic" = "darkgreen",
             "darkblue")
    } else {
      "darkblue"
    },
    edge.label.font = if(params$bold_edge_labels) 2 else 1,
    asize = params$arrow_size,
    curve = if (!is.null(cascade_curve_mat)) cascade_curve_mat else 0.2,
    curveAll = if (!is.null(cascade_curve_mat)) TRUE else FALSE,
    arrowAngle = pi/4,
    open = FALSE,
    title = plot_title,
    title.cex = if (!is.null(plot_title)) 0.9 else NULL,
    border.width = params$node_border_width,
    border.color = params$node_border_color,
    rescale = TRUE,
    normalize = FALSE,
    repulsion = if(!is.null(params$node_spacing)) params$node_spacing else 1
  )

  # Add legends with cleaner styling
  add_plot_legends(strengths, directions, params, vars, continuous_vars, discrete_vars, data)
}

# Calculate plot properties (colors, sizes, etc.)
calculate_plot_properties <- function(boot_selected, vars, continuous_vars, discrete_vars, data, strengths, params, group_name = NULL) {
  n <- length(vars)

  # Edge colors based on strength with more subtle palette
  pal <- brewer.pal(min(9, max(3, length(unique(strengths)))), params$palette)
  strength_breaks <- seq(min(strengths), max(strengths), length.out = length(pal))
  edge_colors <- pal[findInterval(strengths, strength_breaks, all.inside = TRUE)]

  # Node labels - always use full variable names
  labels <- vars

  # Calculate node means
  node_means <- numeric(n)
  names(node_means) <- vars

  if (length(continuous_vars) > 0) {
    cont_means <- colMeans(data[continuous_vars], na.rm = TRUE)
    node_means[continuous_vars] <- cont_means
  }

  # Set discrete variables to neutral value
  if (length(discrete_vars) > 0) {
    if (length(continuous_vars) > 0) {
      node_means[discrete_vars] <- mean(node_means[continuous_vars], na.rm = TRUE)
    } else {
      node_means[discrete_vars] <- 0.5
    }
  }

  # Node sizes and colors
  if (params$scale_nodes && length(continuous_vars) > 0) {
    # Normalize to smaller range for less dramatic differences
    if (length(unique(node_means)) > 1) {
      normalized_means <- (node_means - min(node_means)) / (max(node_means) - min(node_means))
    } else {
      normalized_means <- rep(0.5, length(node_means))
    }
    node_sizes <- normalized_means * (params$max_node_size - params$min_node_size) + params$min_node_size

    # More subtle alpha variation with data type differentiation using user-selected colors
    alpha_vals <- normalized_means * 0.4 + 0.6
    node_colors <- character(n)
    # Convert hex to RGB with alpha
    cont_col_rgb <- col2rgb(params$node_continuous_color) / 255
    disc_col_rgb <- col2rgb(params$node_discrete_color) / 255
    node_colors[vars %in% continuous_vars] <- rgb(cont_col_rgb[1], cont_col_rgb[2], cont_col_rgb[3],
                                                   alpha = alpha_vals[vars %in% continuous_vars])
    node_colors[vars %in% discrete_vars] <- rgb(disc_col_rgb[1], disc_col_rgb[2], disc_col_rgb[3],
                                                 alpha = alpha_vals[vars %in% discrete_vars])

  } else {
    node_sizes <- rep(params$node_size, n)
    # Different colors for different data types using user-selected colors
    node_colors <- character(n)
    node_colors[vars %in% continuous_vars] <- params$node_continuous_color
    node_colors[vars %in% discrete_vars] <- params$node_discrete_color
  }

  return(list(
    edge_colors = edge_colors,
    node_sizes = node_sizes,
    node_colors = node_colors,
    labels = labels
  ))
}

# Add legends to plot
add_plot_legends <- function(strengths, directions, params, vars, continuous_vars, discrete_vars, data) {
  legend_y_pos <- "topright"
  legend_items <- list()

  if (params$show_labels && !is.null(params$edge_display_type)) {
    # Follow edge_display_type selection
    legend_info <- switch(params$edge_display_type,
      "strength" = list(
        text = c("Edge Labels:", "Bootstrap Strength", "(0 = never appears)", "(1 = always appears)"),
        color = "darkblue"
      ),
      "direction" = list(
        text = c("Edge Labels:", "Direction Probability", "(0.5 = undirected)", "(1.0 = strong direction)"),
        color = "darkred"
      ),
      "combined" = list(
        text = c("Edge Labels:", "Strength x Direction", "Combined measure", "Higher = Stronger and/or more directed"),
        color = "darkorchid"
      ),
      "ic" = {
        ic_lbl <- if (!is.null(params$ic_score_type) && params$ic_score_type == "relative")
          c("Edge Labels:", "IC Score (Relative)", "% of total network score", "Higher = More important edge")
        else
          c("Edge Labels:", "IC Score (Absolute)", "BIC/AIC delta from arc removal", "Higher = More important edge")
        list(text = ic_lbl, color = "darkgreen")
      },
      NULL
    )

    if (!is.null(legend_info)) {
      legend("topright",
             legend = legend_info$text,
             col = c(NA, legend_info$color, legend_info$color, legend_info$color),
             lty = c(NA, 1, 1, 1),
             cex = 0.6,
             bty = "n",
             text.col = c("black", legend_info$color, legend_info$color, legend_info$color))
      legend_y_pos <- "topleft"
    }
  }

  # Add data type legend if mixed
  if (length(continuous_vars) > 0 && length(discrete_vars) > 0) {
    legend(legend_y_pos,
           legend = c("Node Types:", "Continuous", "Discrete"),
           col = c(NA, params$node_continuous_color, params$node_discrete_color),
           pch = c(NA, 19, 19),
           pt.cex = 1.5,
           cex = 0.6,
           bty = "n",
           text.col = c("black", "#1976d2", "#f57900"))
  }
}

# ============================================================================
# FOLDED TEMPORAL GRAPH FUNCTIONS
# ============================================================================

# Detect whether selected variables have a temporal prefix/suffix structure.
# Returns list(type, time_slices, base_names) or NULL if not detected.
detect_temporal_structure <- function(vars) {
  if (length(vars) < 2) return(NULL)

  has_underscore <- grepl("_", vars)

  prefix_result <- NULL
  suffix_result <- NULL

  if (sum(has_underscore) >= 2) {
    # --- Prefix candidate (NORM_DIS, LD_DIS): split on first underscore ---
    prefix_parts <- ifelse(has_underscore, sub("^([^_]+)_.*$", "\\1", vars), NA)
    base_parts   <- ifelse(has_underscore, sub("^[^_]+_(.+)$", "\\1",   vars), NA)
    valid <- which(!is.na(prefix_parts) & !is.na(base_parts))
    if (length(valid) >= 2) {
      unique_prefixes   <- unique(prefix_parts[valid])
      base_prefix_count <- table(base_parts[valid])
      shared_bases      <- names(base_prefix_count[base_prefix_count >= 2])
      if (length(unique_prefixes) >= 2 && length(shared_bases) >= 1)
        prefix_result <- list(type = "prefix", time_slices = unique_prefixes, base_names = shared_bases)
    }

    # --- Suffix candidate (VAR_T0, VAR_T1): split on last underscore ---
    suffix_parts     <- ifelse(has_underscore, sub("^.*_([^_]+)$", "\\1", vars), NA)
    base_from_suffix <- ifelse(has_underscore, sub("_[^_]+$",      "",    vars), NA)
    valid2 <- which(!is.na(suffix_parts) & !is.na(base_from_suffix))
    if (length(valid2) >= 2) {
      suffix_counts <- table(suffix_parts[valid2])
      shared_suf    <- names(suffix_counts[suffix_counts >= 2])
      base_sf_count <- table(base_from_suffix[valid2])
      shared_bases2 <- names(base_sf_count[base_sf_count >= 2])
      if (length(shared_suf) >= 2 && length(shared_bases2) >= 1)
        suffix_result <- list(type = "suffix", time_slices = shared_suf, base_names = shared_bases2)
    }

    # --- Choose: fewer time slices = more likely the true temporal labels.
    #     Ties go to prefix (prefix tokens are usually the time labels). ---
    chosen <- NULL
    if (!is.null(prefix_result) && !is.null(suffix_result)) {
      chosen <- if (length(suffix_result$time_slices) < length(prefix_result$time_slices))
                  suffix_result else prefix_result
    } else {
      chosen <- if (!is.null(prefix_result)) prefix_result else suffix_result
    }

    if (!is.null(chosen)) {
      message("[TemporalDetect] ", chosen$type, " structure found: slices=",
              paste(chosen$time_slices, collapse=","), "  shared_bases=", length(chosen$base_names))
      return(chosen)
    }
  }

  # --- Digit-suffix detection (SQ0/SQ1, DIS0/DIS1) - no underscore needed ---
  has_digit_sfx <- grepl("\\d+$", vars) & nchar(sub("\\d+$", "", vars)) > 0
  if (sum(has_digit_sfx) >= 2) {
    digit_slices <- ifelse(has_digit_sfx, regmatches(vars, regexpr("\\d+$", vars)), NA)
    digit_bases  <- ifelse(has_digit_sfx, sub("\\d+$", "", vars), NA)
    valid3 <- which(!is.na(digit_slices) & !is.na(digit_bases))
    if (length(valid3) >= 2) {
      unique_digits <- unique(digit_slices[valid3])
      base_cnt      <- table(digit_bases[valid3])
      shared_bases3 <- names(base_cnt[base_cnt >= 2])
      if (length(unique_digits) >= 2 && length(shared_bases3) >= 1) {
        message("[TemporalDetect] digit-suffix structure found: slices=",
                paste(unique_digits, collapse=","), "  shared_bases=", length(shared_bases3))
        return(list(type = "digit_suffix", time_slices = unique_digits, base_names = shared_bases3))
      }
    }
  }

  message("[TemporalDetect] no temporal structure found in vars: ", paste(head(vars,10), collapse=","))
  NULL
}

# Extract the base variable name given the temporal structure.
get_temporal_base <- function(var, ts) {
  if      (ts$type == "prefix")       sub("^[^_]+_(.+)$", "\\1", var)
  else if (ts$type == "suffix")       sub("_[^_]+$",      "",    var)
  else                                 sub("\\d+$",         "",    var)   # digit_suffix
}

# Extract the time-slice label from a variable name.
get_temporal_slice <- function(var, ts) {
  if      (ts$type == "prefix")       sub("^([^_]+)_.*$", "\\1", var)
  else if (ts$type == "suffix")       sub("^.*_([^_]+)$", "\\1", var)
  else                                 regmatches(var, regexpr("\\d+$", var))  # digit_suffix
}

# Test whether a variable belongs to a given time slice.
var_in_slice <- function(var, slice, ts) {
  switch(ts$type,
    prefix       = startsWith(var, paste0(slice, "_")),
    suffix       = endsWith(var,   paste0("_", slice)),
    digit_suffix = endsWith(var,   as.character(slice)) && grepl("\\d+$", var)
  )
}

# Infer temporal ordering from blacklisted edges.
infer_temporal_order <- function(ts, blacklist) {
  slices <- ts$time_slices
  if (length(slices) < 2 || is.null(blacklist) || nrow(blacklist) == 0) return(slices)

  # Build precedence tally: earlier_than[A, B] = how often A precedes B
  earlier_than <- matrix(0, length(slices), length(slices), dimnames = list(slices, slices))

  for (i in seq_len(nrow(blacklist))) {
    fv <- blacklist$from[i]; tv <- blacklist$to[i]
    fs <- NA; ts_sl <- NA
    for (sl in slices) {
      matched_f <- var_in_slice(fv, sl, ts)
      matched_t <- var_in_slice(tv, sl, ts)
      if (matched_f) fs    <- sl
      if (matched_t) ts_sl <- sl
    }
    if (!is.na(fs) && !is.na(ts_sl) && fs != ts_sl) {
      # from_slice -> to_slice is blacklisted  =>  to_slice is EARLIER than from_slice
      earlier_than[ts_sl, fs] <- earlier_than[ts_sl, fs] + 1
    }
  }

  earlier_counts <- rowSums(earlier_than > 0)
  message("[TemporalOrder] earlier_counts: ",
          paste(names(earlier_counts), earlier_counts, sep = "=", collapse = ", "))
  if (length(unique(earlier_counts)) == 1) {
    message("[TemporalOrder] Could not determine order from blacklist; returning: ",
            paste(slices, collapse = ", "))
    return(slices)
  }
  # Slices "earlier than" more others should come first => decreasing = TRUE
  ordered <- slices[order(earlier_counts, decreasing = TRUE)]
  message("[TemporalOrder] Inferred order: ", paste(ordered, collapse = " -> "))
  ordered
}

# Build folded edge data frame from bootstrap-selected edges.
build_folded_edges <- function(boot_selected, ts, ordered_slices, weight_type = "combined") {
  if (is.null(ts) || nrow(boot_selected) == 0) return(NULL)
  rows <- list()
  for (i in seq_len(nrow(boot_selected))) {
    fv <- boot_selected$from[i]; tv <- boot_selected$to[i]
    fs <- NA; tsl <- NA
    for (sl in ordered_slices) {
      mf <- var_in_slice(fv, sl, ts)
      mt <- var_in_slice(tv, sl, ts)
      if (mf) fs  <- sl
      if (mt) tsl <- sl
    }
    if (is.na(fs) || is.na(tsl)) next
    fb <- get_temporal_base(fv, ts); tb <- get_temporal_base(tv, ts)
    str_val <- as.numeric(boot_selected$strength[i])
    dir_val <- as.numeric(boot_selected$direction[i])
    w <- switch(weight_type,
      strength  = str_val,
      direction = dir_val,
      combined  = str_val * dir_val,
      str_val * dir_val
    )
    etype <- if (fs == tsl) paste0("contemp_", fs) else if (fb == tb) "self_loop" else "cross_lag"
    rows[[length(rows) + 1]] <- data.frame(
      from = fb, to = tb, from_slice = fs, to_slice = tsl,
      edge_type = etype, weight = w, strength = str_val, direction = dir_val,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  aggregate_folded_edges(do.call(rbind, rows))
}

# Collapse parallel edges sharing the same (from, to, edge_type) across
# multiple transitions by averaging weight, strength, and direction.
aggregate_folded_edges <- function(folded) {
  if (is.null(folded) || nrow(folded) == 0) return(folded)
  keys <- paste(folded$from, folded$to, folded$edge_type, sep = "|||")
  result <- do.call(rbind, lapply(split(seq_len(nrow(folded)), keys), function(idx) {
    grp <- folded[idx, , drop = FALSE]
    data.frame(
      from          = grp$from[1],
      to            = grp$to[1],
      from_slice    = paste(sort(unique(grp$from_slice)), collapse = "/"),
      to_slice      = paste(sort(unique(grp$to_slice)),   collapse = "/"),
      edge_type     = grp$edge_type[1],
      weight        = mean(grp$weight,    na.rm = TRUE),
      strength      = mean(grp$strength,  na.rm = TRUE),
      direction     = mean(grp$direction, na.rm = TRUE),
      n_transitions = nrow(grp),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}

# Main function: plot the folded temporal graph using qgraph.
plot_folded_temporal_graph <- function(boot_selected, vars, blacklist, params) {

  ts <- detect_temporal_structure(vars)
  if (is.null(ts)) {
    plot.new(); par(mar = c(2,2,2,2))
    text(0.5, 0.5, "No temporal structure detected.\nVariables need shared prefixes or suffixes\n(e.g. NORM_VAR / LD_VAR or VAR_1 / VAR_2)",
         cex = 1.1, col = "firebrick", adj = 0.5); return(invisible(NULL))
  }

  ordered_slices <- infer_temporal_order(ts, blacklist)
  wt <- if (!is.null(params$temporal_weight_type)) params$temporal_weight_type else "combined"
  folded <- build_folded_edges(boot_selected, ts, ordered_slices, wt)

  # All base nodes always shown, regardless of edge connectivity
  base_nodes <- unique(ts$base_names)

  # --- Apply edge visibility toggles ---
  if (!is.null(folded)) {
    if (!isTRUE(params$show_crosslag_edges)) {
      folded <- folded[!folded$edge_type %in% c("cross_lag", "self_loop"), ]
    }
    show_ct <- if (!is.null(params$show_contemp_slices)) params$show_contemp_slices else list()
    for (sl in names(show_ct)) {
      if (!isTRUE(show_ct[[sl]])) folded <- folded[folded$edge_type != paste0("contemp_", sl), ]
    }
    if (nrow(folded) == 0) folded <- NULL

    # --- Apply minimum edge weight filter ---
    if (!is.null(folded)) {
      min_w <- if (!is.null(params$min_edge_weight) && is.numeric(params$min_edge_weight))
                  params$min_edge_weight else 0
      if (min_w > 0) {
        folded <- folded[abs(folded$weight) >= min_w, ]
        if (nrow(folded) == 0) folded <- NULL
      }
    }
  }
  n          <- length(base_nodes)

  # Circular layout (start at top, go counter-clockwise)
  angles     <- seq(pi / 2, pi / 2 + 2 * pi * (1 - 1 / n), length.out = n)
  qg_layout  <- cbind(cos(angles), sin(angles))

  # --- Edge weights -> widths ---
  e_type   <- if (!is.null(folded)) folded$edge_type  else character(0)
  e_weight <- if (!is.null(folded)) abs(folded$weight) else numeric(0)
  if (length(e_weight) > 0 && max(e_weight, na.rm = TRUE) > min(e_weight, na.rm = TRUE)) {
    e_width <- (e_weight - min(e_weight, na.rm = TRUE)) /
               (max(e_weight, na.rm = TRUE) - min(e_weight, na.rm = TRUE)) * 4.5 + 0.5
  } else {
    e_width <- rep(2, length(e_weight))
  }
  mult    <- if (!is.null(params$edge_multiplier) && is.numeric(params$edge_multiplier) && params$edge_multiplier > 0)
               params$edge_multiplier else 1
  e_width <- pmax(0.1, e_width * mult)

  # --- Colours ---
  crosslag_color <- if (!is.null(params$temporal_crosslag_color) && nzchar(params$temporal_crosslag_color))
                      params$temporal_crosslag_color else "#ED7474"
  default_contemp  <- c("#FFBF6B", "#C6DB72", "#9B59B6", "#F39C12", "#1ABC9C")
  contemp_colors_p <- if (!is.null(params$temporal_contemp_colors)) params$temporal_contemp_colors else list()
  slice_colors <- setNames(
    lapply(seq_along(ordered_slices), function(i) {
      sl  <- ordered_slices[i]
      col <- contemp_colors_p[[sl]]
      if (!is.null(col) && nzchar(col)) col else default_contemp[min(i, length(default_contemp))]
    }),
    ordered_slices
  )

  e_color <- vapply(e_type, function(et) {
    if (et %in% c("cross_lag", "self_loop")) return(crosslag_color)
    sl <- sub("^contemp_", "", et)
    if (sl %in% names(slice_colors)) slice_colors[[sl]] else "#888888"
  }, character(1))

  # Dashed or solid for contemporaneous edges based on user toggle
  contemp_lty <- if (isTRUE(params$contemp_dotted)) 2L else 1L
  e_lty <- ifelse(grepl("^contemp_", e_type), contemp_lty, 1L)

  # --- qgraph edge list (integer from/to indices) ---
  if (!is.null(folded) && nrow(folded) > 0) {
    from_idx <- match(folded$from, base_nodes)
    to_idx   <- match(folded$to,   base_nodes)
    qg_edges <- cbind(from_idx, to_idx, e_weight)
  } else {
    qg_edges <- matrix(nrow = 0, ncol = 3)
  }

  # Self-loop rotation: per-node, point outward (= same angle as node on circle)
  loop_rotation <- angles   # angles[i] is the outward direction for node i

  lbl_cex <- if (!is.null(params$label_size)) params$label_size else 0.85

  message("[FoldedPlot] Using qgraph. nodes=", paste(base_nodes, collapse=","),
          "  edges=", if (!is.null(folded)) nrow(folded) else 0,
          "  loop_rotation=", paste(round(loop_rotation, 2), collapse=","))

  node_fill   <- if (!is.null(params$temporal_node_color)  && nzchar(params$temporal_node_color))
                   params$temporal_node_color  else "#E5F2FF"
  node_border <- if (!is.null(params$temporal_node_border) && nzchar(params$temporal_node_border))
                   params$temporal_node_border else "#B3CBE6"
  node_size   <- if (!is.null(params$temporal_node_size) && is.numeric(params$temporal_node_size) &&
                     params$temporal_node_size > 0) params$temporal_node_size else 10

  qgraph::qgraph(
    qg_edges,
    directed      = TRUE,
    layout        = qg_layout,
    vsize         = node_size,
    color         = node_fill,
    border.color  = node_border,
    border.width  = 2,
    labels        = base_nodes,
    label.scale   = FALSE,
    label.cex     = lbl_cex,
    edge.color    = e_color,
    lty           = e_lty,
    esize         = e_width,
    asize         = ifelse(grepl("^contemp_", e_type) & isTRUE(params$hide_contemp_arrows), 0,
                          if (!is.null(params$arrow_size)) params$arrow_size else 4),
    edge.labels   = if (isTRUE(params$show_edge_labels)) round(e_weight, 2) else FALSE,
    edge.label.cex = if (!is.null(params$edge_label_size) && is.numeric(params$edge_label_size))
                       params$edge_label_size else 0.7,
    loopRotation  = loop_rotation,
    mar           = c(8, 8, 8, 8)
  )

  # --- Legend (drawn on top of qgraph output) ---
  leg_labels <- c("Self-loop / Cross-lag (solid)")
  leg_cols   <- c(crosslag_color)
  leg_lty    <- c(1)
  contemp_used <- unique(grep("^contemp_", e_type, value = TRUE))
  for (sl in ordered_slices) {
    ct <- paste0("contemp_", sl)
    if (ct %in% contemp_used) {
      col <- if (sl %in% names(slice_colors)) slice_colors[[sl]] else "#888888"
      leg_labels <- c(leg_labels, paste0("Contemporaneous: ", sl, " (dashed)"))
      leg_cols   <- c(leg_cols, col)
      leg_lty    <- c(leg_lty, 2)
    }
  }
  legend("topleft", legend = leg_labels, col = leg_cols, lty = leg_lty,
         lwd = 2.5, bty = "n", cex = 0.75)

  mtext(paste("Time order:", paste(ordered_slices, collapse = " -> ")),
        side = 1, cex = 0.8, col = "gray40")

  invisible(NULL)
}

# Run Bayesian network analysis
run_analysis <- function(data, params) {

  incProgress(0.05, detail = 'Preparing data...')

  # Store original data for reference
  original_data <- data
  analysis_data <- data

  # Store the missing data method from params
  missing_method <- if(!is.null(params$missing_method)) params$missing_method else "listwise"
  missing_warnings <- character()

  # Check if we have any missing data
  has_missing <- any(is.na(data))
  n_missing <- sum(!complete.cases(data))
  pct_missing <- round((n_missing / nrow(data)) * 100, 1)

  # Initialize net variable
  net <- NULL

  if (missing_method == "pairwise") {
    # Pairwise deletion (limited support)
    warning("Pairwise deletion has limited support in bnlearn. Results may be approximate.")
    analysis_data <- data

  } else if (missing_method == "none") {
    # Keep all data as is
    analysis_data <- data

  } else {
    # Default: listwise deletion
    if (has_missing) {
      analysis_data <- na.omit(data)
      incProgress(0.02, detail = paste('Removed', nrow(data) - nrow(analysis_data), 'rows with missing data'))
    } else {
      analysis_data <- data
    }
  }

  # Now learn the network
  if (is.null(net)) {
    incProgress(0.1, detail = paste('Learning structure with', params$algorithm, 'algorithm...'))

    # Prepare arguments
    base_args <- list(x = analysis_data)

    # Add blacklist and whitelist if provided
    if (!is.null(params$blacklist)) {
      base_args$blacklist <- params$blacklist
      incProgress(0.02, detail = paste('Applied', nrow(params$blacklist), 'blacklist constraints'))
    }
    if (!is.null(params$whitelist)) {
      base_args$whitelist <- params$whitelist
      incProgress(0.02, detail = paste('Applied', nrow(params$whitelist), 'whitelist constraints'))
    }

    # Learn the network based on algorithm type
    tryCatch({
      if (params$algorithm %in% c("hc", "tabu")) {
        # Pure score-based algorithms
        args <- c(base_args, list(
          score = params$score,
          restart = params$restarts,
          perturb = params$perturb
        ))
        net <- do.call(params$algorithm, args)
      } else if (params$algorithm == "mmhc") {
        # Hybrid algorithm
        args <- c(base_args, list(
          maximize.args = list(score = params$score)
        ))
        net <- do.call("mmhc", args)
      } else {
        # Constraint-based algorithms
        net <- do.call(params$algorithm, base_args)
      }

      # Verify the network is valid
      if (!inherits(net, "bn")) {
        stop("Network learning failed to produce valid bn object")
      }

    }, error = function(e) {
      stop(paste("Failed to learn network structure:", e$message))
    })
  }

  incProgress(0.15, detail = paste('Initial network learned with', narcs(net), 'arcs'))

  # Bootstrap analysis - different parameter sets by algorithm type
  boot_base_args <- list()
  if (!is.null(params$blacklist)) {
    boot_base_args$blacklist <- params$blacklist
  }
  if (!is.null(params$whitelist)) {
    boot_base_args$whitelist <- params$whitelist
  }

  if (params$algorithm %in% c("hc", "tabu")) {
    # Pure score-based: include score, restart, and perturb
    boot_args <- c(boot_base_args, list(
      score = params$score,
      restart = params$boot_restarts,
      perturb = params$boot_perturb
    ))
  } else if (params$algorithm == "mmhc") {
    # Hybrid: pass score through maximize.args
    boot_args <- c(boot_base_args, list(
      maximize.args = list(score = params$score)
    ))
  } else {
    # Constraint-based: no additional parameters beyond constraints
    boot_args <- boot_base_args
  }

  # Detailed bootstrap progress
  incProgress(0.05, detail = paste('Starting bootstrap with', params$boot_r, 'iterations...'))

  # Estimate time and show progress
  start_time <- Sys.time()

  # Run a small sample first to estimate time
  if (params$boot_r >= 50) {
    incProgress(0.05, detail = 'Running sample bootstrap to estimate time...')

    sample_boot <- boot.strength(analysis_data, R = min(10, max(2, params$boot_r %/% 20)),
                                algorithm = params$algorithm,
                                algorithm.args = boot_args)

    elapsed_sample <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    estimated_total_time <- (elapsed_sample / min(10, max(2, params$boot_r %/% 20))) * params$boot_r

    incProgress(0.05, detail = paste('Estimated time:', round(estimated_total_time, 1), 'seconds'))
  }

  # Run full bootstrap
  incProgress(0.1, detail = paste('Running', params$boot_r, 'bootstrap iterations...'))

  boot <- boot.strength(analysis_data, R = params$boot_r,
                       algorithm = params$algorithm,
                       algorithm.args = boot_args)

  bootstrap_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  incProgress(0.4, detail = paste('Bootstrap completed in', round(bootstrap_time, 1), 'seconds'))

  # Create averaged network
  incProgress(0.1, detail = 'Creating bootstrap-averaged network...')
  avg <- averaged.network(boot, threshold = params$threshold)

  # Force full directionality if requested OR if threshold = 0 (to avoid undirected edges)
  if (!is.null(params$force_directionality) && params$force_directionality) {
    incProgress(0.03, detail = 'Forcing full directionality with cextend...')
    avg <- cextend(avg)
  } else if (params$threshold == 0 || params$direction == 0) {
    # Automatically extend when using zero thresholds to avoid undirected edges
    incProgress(0.03, detail = 'Auto-extending network (threshold = 0 detected)...')
    avg <- cextend(avg)
    message("Note: Automatically applied cextend() because threshold or direction = 0. This orients all edges to create a fully directed network.")
  }

  # Compute arc strengths on the finalized averaged network (correct reference)
  incProgress(0.05, detail = 'Computing arc strengths on averaged network...')
  astr <- tryCatch(
    arc.strength(avg, analysis_data, criterion = params$score),
    error = function(e) {
      warning("Could not compute arc strengths on averaged network: ", e$message)
      NULL
    }
  )

  incProgress(0.05, detail = 'Filtering significant edges...')
  sel_boot <- subset(boot, strength >= params$threshold & direction >= params$direction)

  # Merge arc strength score deltas with boot_selected
  if (!is.null(sel_boot) && is.data.frame(sel_boot) && nrow(sel_boot) > 0 && !is.null(astr) && is.data.frame(astr)) {
    sel_boot <- merge(sel_boot, astr[, c("from", "to", "strength")],
                      by = c("from", "to"), all.x = TRUE, suffixes = c("", "_score"))
    if ("strength_score" %in% names(sel_boot)) {
      names(sel_boot)[names(sel_boot) == "strength_score"] <- "score_delta"
    }
  }

  incProgress(0.05, detail = paste('Found', nrow(sel_boot), 'significant edges'))

  # Compute statistics
  incProgress(0.05, detail = 'Computing final statistics...')

  # Safely compute scores (handle potential errors from undirected edges)
  score_orig <- tryCatch(
    score(net, analysis_data, type = params$score),
    error = function(e) {
      warning("Could not compute score for original network: ", e$message)
      NA
    }
  )

  score_avg <- tryCatch(
    score(avg, analysis_data, type = params$score),
    error = function(e) {
      warning("Could not compute score for averaged network: ", e$message,
              "\nThis may occur with undirected edges. Try enabling 'Force Full Directionality'.")
      NA
    }
  )

  # Add relative score_delta (each arc's contribution as % of total network score)
  if (!is.null(sel_boot) && "score_delta" %in% names(sel_boot) &&
      !is.null(score_avg) && !is.na(score_avg) && score_avg != 0) {
    sel_boot$score_delta_relative <- sel_boot$score_delta / abs(score_avg) * 100
  } else if (!is.null(sel_boot)) {
    sel_boot$score_delta_relative <- NA_real_
  }

  stats <- list(
    nodes = nnodes(net),
    arcs = narcs(net),
    avg_nodes = nnodes(avg),
    avg_arcs = narcs(avg),
    score_orig = score_orig,
    score_avg = score_avg,
    blacklist_count = if(!is.null(params$blacklist)) nrow(params$blacklist) else 0,
    whitelist_count = if(!is.null(params$whitelist)) nrow(params$whitelist) else 0,
    bootstrap_time = bootstrap_time,
    continuous_vars = params$continuous_vars,
    discrete_vars = params$discrete_vars,
    missing_method = missing_method,
    missing_method_requested = params$missing_method,
    missing_method_warnings = missing_warnings,
    n_missing = sum(is.na(original_data)),
    n_complete_cases = nrow(na.omit(original_data)),
    n_after_processing = nrow(analysis_data),
    n_original = nrow(original_data)
  )

  incProgress(0.05, detail = 'Analysis complete!')

  return(list(
    net = net,
    astr = astr,
    boot = boot,
    boot_selected = sel_boot,
    avg = avg,
    data = analysis_data,
    original_data = original_data,
    stats = stats,
    blacklist = params$blacklist,
    whitelist = params$whitelist
  ))
}

# MODIFIED: Complete generate_r_code function with full R script generation
generate_r_code <- function(params, continuous_vars, discrete_vars, data_filename = "your_data.csv", split_var = NULL) {

  # Ensure params exists and has required fields
  if (is.null(params)) {
    stop("params argument is required")
  }

  # Set defaults for missing params
  if (is.null(params$missing_method)) params$missing_method <- "listwise"

  all_vars <- c(continuous_vars, discrete_vars)

  # Build the R script
  code <- paste0(
    "# Bayesian Network Analysis - Standalone R Code\n",
    "# Generated from bnlearn DAG Analysis App\n",
    "# ", Sys.time(), "\n\n",

    "# Load required libraries\n",
    "library(bnlearn)\n",
    "library(qgraph)\n",
    "library(RColorBrewer)\n",
    if (params$layout == "cascade") "library(igraph)  # Required for cascade layout (Katz centrality)\n" else "",
    "# library(visNetwork)  # Uncomment for interactive plots\n",
    "\n",

    "# Load and prepare data\n",
    "data <- read.csv('", data_filename, "')\n",
    "all_vars <- c(", paste0("'", all_vars, "'", collapse = ", "), ")\n",
    "continuous_vars <- c(", paste0("'", continuous_vars, "'", collapse = ", "), ")\n",
    "discrete_vars <- c(", paste0("'", discrete_vars, "'", collapse = ", "), ")\n"
  )

  # Add split variable handling if specified
  if (!is.null(split_var) && split_var != "") {
    code <- paste0(code,
      "\n# Split variable for subgroup analysis\n",
      "split_var <- '", split_var, "'\n",
      "split_levels <- unique(data[[split_var]])\n",
      "split_results <- list()\n\n",
      "# Perform analysis for each subgroup\n",
      "for (level in split_levels) {\n",
      "  cat('\\n========================================\\n')\n",
      "  cat('Analyzing subgroup:', split_var, '=', level, '\\n')\n",
      "  cat('========================================\\n')\n",
      "  \n",
      "  # Filter data for this subgroup\n",
      "  subgroup_data <- data[data[[split_var]] == level, ]\n",
      "  analysis_data <- subgroup_data[, all_vars]\n"
    )
  } else {
    code <- paste0(code,
      "\nanalysis_data <- data[, all_vars]\n"
    )
  }

  # Add missing data handling
  missing_method <- if(!is.null(params$missing_method)) params$missing_method else "listwise"

  code <- paste0(code, "\n",
    "# Missing data handling method: ", missing_method, "\n"
  )

  if (missing_method == "pairwise") {
    code <- paste0(code,
      "# Using pairwise deletion (limited support in bnlearn)\n",
      "cat('Note: Pairwise deletion has limited support. Results may be approximate.\\n')\n",
      "# bnlearn will handle missing data internally where possible\n\n"
    )

  } else if (missing_method == "none") {
    code <- paste0(code,
      "# Keep all data (no removal of missing values)\n",
      "cat('Keeping all data including missing values.\\n')\n\n"
    )

  } else {
    # Default: listwise deletion
    code <- paste0(code,
      "# Remove missing values (listwise deletion)\n",
      "analysis_data <- na.omit(analysis_data)\n",
      "cat('Removed', nrow(",
      if (!is.null(split_var) && split_var != "") "subgroup_data[, all_vars]" else "data[, all_vars]",
      ") - nrow(analysis_data), 'rows with missing data\\n')\n\n"
    )
  }

  # Add data preparation
  code <- paste0(code,
    "# Prepare mixed data types\n",
    "if (length(continuous_vars) > 0) {\n",
    "  # Convert continuous variables to numeric\n",
    "  cont_data <- as.matrix(analysis_data[, continuous_vars, drop = FALSE])\n",
    "  storage.mode(cont_data) <- 'numeric'\n",
    "  analysis_data[, continuous_vars] <- as.data.frame(cont_data)\n",
    "}\n\n",

    "if (length(discrete_vars) > 0) {\n",
    "  # Convert discrete variables to factors\n",
    "  for (var in discrete_vars) {\n",
    "    analysis_data[[var]] <- as.factor(analysis_data[[var]])\n",
    "  }\n",
    "}\n\n"
  )

  # Add constraints if any
  if (!is.null(params$blacklist) || !is.null(params$whitelist)) {
    code <- paste0(code, "# Network constraints\n")

    if (!is.null(params$blacklist) && nrow(params$blacklist) > 0) {
      bl_code <- paste0(
        "blacklist <- data.frame(\n",
        "  from = c(", paste0("'", params$blacklist$from, "'", collapse = ", "), "),\n",
        "  to = c(", paste0("'", params$blacklist$to, "'", collapse = ", "), ")\n",
        ")\n\n"
      )
      code <- paste0(code, bl_code)
    } else {
      code <- paste0(code, "blacklist <- NULL\n\n")
    }

    if (!is.null(params$whitelist) && nrow(params$whitelist) > 0) {
      wl_code <- paste0(
        "whitelist <- data.frame(\n",
        "  from = c(", paste0("'", params$whitelist$from, "'", collapse = ", "), "),\n",
        "  to = c(", paste0("'", params$whitelist$to, "'", collapse = ", "), ")\n",
        ")\n\n"
      )
      code <- paste0(code, wl_code)
    } else {
      code <- paste0(code, "whitelist <- NULL\n\n")
    }
  } else {
    code <- paste0(code, "blacklist <- NULL\nwhitelist <- NULL\n\n")
  }

  # Network learning section
  code <- paste0(code,
    "# Learn network structure\n",
    "if (!exists('net') || is.null(net)) {\n",
    "  cat('Learning network structure...\\n')\n"
  )

  # Add algorithm-specific learning code
  if (params$algorithm %in% c("hc", "tabu")) {
    code <- paste0(code,
      "  # Score-based algorithm with restart and perturbation\n",
      "  net <- ", params$algorithm, "(analysis_data,\n",
      "    score = '", params$score, "',\n",
      "    restart = ", params$restarts, ",\n",
      "    perturb = ", params$perturb,
      if (!is.null(params$blacklist) || !is.null(params$whitelist)) ",\n    blacklist = blacklist,\n    whitelist = whitelist" else "",
      ")\n"
    )
  } else if (params$algorithm == "mmhc") {
    code <- paste0(code,
      "  # Hybrid algorithm (MMHC)\n",
      "  net <- mmhc(analysis_data,\n",
      "    maximize.args = list(score = '", params$score, "')",
      if (!is.null(params$blacklist) || !is.null(params$whitelist)) ",\n    blacklist = blacklist,\n    whitelist = whitelist" else "",
      ")\n"
    )
  } else {
    code <- paste0(code,
      "  # Constraint-based algorithm\n",
      "  net <- ", params$algorithm, "(analysis_data",
      if (!is.null(params$blacklist) || !is.null(params$whitelist)) ",\n    blacklist = blacklist,\n    whitelist = whitelist" else "",
      ")\n"
    )
  }

  code <- paste0(code, "}\n\n")  # Close the if (!exists('net')) block

  # Add arc strength calculation (on averaged network — same reference as the app)
  code <- paste0(code,
    "# Calculate arc strengths relative to the averaged network\n",
    "arc_strengths <- arc.strength(avg_net, analysis_data, criterion = '", params$score, "')\n",
    "total_score <- bnlearn::score(avg_net, analysis_data, type = '", params$score, "')\n",
    "arc_strengths$score_delta_relative <- arc_strengths$strength / abs(total_score) * 100\n",
    "cat('Averaged network has', narcs(avg_net), 'arcs\\n')\n\n"
  )

  # Bootstrap analysis
  code <- paste0(code,
    "# Bootstrap analysis\n",
    "cat('Running bootstrap with ", params$boot_r, " iterations...\\n')\n"
  )

  # Build bootstrap arguments based on algorithm type
  if (params$algorithm %in% c("hc", "tabu")) {
    boot_args <- paste0(
      "  algorithm.args = list(\n",
      "    score = '", params$score, "',\n",
      "    restart = ", params$boot_restarts, ",\n",
      "    perturb = ", params$boot_perturb,
      if (!is.null(params$blacklist) || !is.null(params$whitelist)) ",\n    blacklist = blacklist,\n    whitelist = whitelist" else "",
      "\n  )"
    )
  } else if (params$algorithm == "mmhc") {
    boot_args <- paste0(
      "  algorithm.args = list(\n",
      "    maximize.args = list(score = '", params$score, "')",
      if (!is.null(params$blacklist) || !is.null(params$whitelist)) ",\n    blacklist = blacklist,\n    whitelist = whitelist" else "",
      "\n  )"
    )
  } else {
    if (!is.null(params$blacklist) || !is.null(params$whitelist)) {
      boot_args <- paste0(
        "  algorithm.args = list(\n",
        "    blacklist = blacklist,\n",
        "    whitelist = whitelist\n",
        "  )"
      )
    } else {
      boot_args <- "  algorithm.args = list()"
    }
  }

  code <- paste0(code,
    "boot_strength <- boot.strength(analysis_data,\n",
    "  R = ", params$boot_r, ",\n",
    "  algorithm = '", params$algorithm, "',\n",
    boot_args, "\n",
    ")\n\n"
  )

  # Create averaged network
  code <- paste0(code,
    "# Create bootstrap-averaged network\n",
    "avg_net <- averaged.network(boot_strength, threshold = ", params$threshold, ")\n"
  )

  if (!is.null(params$force_directionality) && params$force_directionality) {
    code <- paste0(code,
      "# Force full directionality\n",
      "avg_net <- cextend(avg_net)\n"
    )
  }

  code <- paste0(code,
    "\n# Select significant edges\n",
    "significant_edges <- subset(boot_strength,\n",
    "  strength >= ", params$threshold, " & direction >= ", params$direction, ")\n",
    "cat('Found', nrow(significant_edges), 'significant edges\\n')\n\n"
  )

  # Model evaluation
  code <- paste0(code,
    "# Model evaluation\n",
    "score_original <- score(net, analysis_data, type = '", params$score, "')\n",
    "score_averaged <- score(avg_net, analysis_data, type = '", params$score, "')\n",
    "cat('Original network score:', score_original, '\\n')\n",
    "cat('Averaged network score:', score_averaged, '\\n')\n\n"
  )

  # Visualization code
  code <- paste0(code,
    "# Visualization\n",
    "# Using qgraph for static visualization\n",
    "if (nrow(significant_edges) > 0) {\n",
    "  # Prepare adjacency matrix\n",
    "  nodes <- unique(c(continuous_vars, discrete_vars))\n",
    "  n <- length(nodes)\n",
    "  adj_matrix <- matrix(0, n, n, dimnames = list(nodes, nodes))\n",
    "  \n",
    "  for (i in 1:nrow(significant_edges)) {\n",
    "    adj_matrix[significant_edges$from[i], significant_edges$to[i]] <- significant_edges$strength[i]\n",
    "  }\n",
    "  \n",
    "  # Create network plot\n",
    "  qgraph(adj_matrix,\n",
    "    directed = TRUE,\n",
    "    layout = '", params$layout, "',\n",
    "    edge.width = ", params$edge_width, ",\n",
    "    vsize = ", params$node_size, ",\n",
    "    asize = ", params$arrow_size, ",\n",
    "    edge.color = '", params$palette, "',\n",
    "    labels = nodes,\n",
    "    label.cex = ", params$label_size, ",\n",
    "    label.scale = FALSE,\n",
    "    edge.labels = ", if(params$show_labels) "round(significant_edges$strength, 3)" else "FALSE", ",\n",
    "    title = 'Bayesian Network DAG'\n",
    "  )\n",
    "} else {\n",
    "  cat('No edges meet the threshold criteria.\\n')\n",
    "}\n"
  )

  # Add results summary
  code <- paste0(code,
    "\n# Results summary\n",
    "cat('\\n========================================\\n')\n",
    "cat('ANALYSIS SUMMARY\\n')\n",
    "cat('========================================\\n')\n",
    "cat('Algorithm:', '", params$algorithm, "'\\n')\n",
    "cat('Score type:', '", params$score, "'\\n')\n",
    "cat('Missing data method:', '", missing_method, "'\\n')\n",
    "cat('Bootstrap iterations:', ", params$boot_r, "\\n')\n",
    "cat('Strength threshold:', ", params$threshold, "\\n')\n",
    "cat('Direction threshold:', ", params$direction, "\\n')\n",
    "cat('Nodes:', nnodes(net), '\\n')\n",
    "cat('Original arcs:', narcs(net), '\\n')\n",
    "cat('Averaged arcs:', narcs(avg_net), '\\n')\n",
    "cat('Significant edges:', nrow(significant_edges), '\\n')\n"
  )

  # Close the loop if split analysis
  if (!is.null(split_var) && split_var != "") {
    code <- paste0(code,
      "\n  # Store results for this group\n",
      "  split_results[[as.character(level)]] <- list(\n",
      "    net = net,\n",
      "    avg_net = avg_net,\n",
      "    boot_strength = boot_strength,\n",
      "    significant_edges = significant_edges,\n",
      "    score_original = score_original,\n",
      "    score_averaged = score_averaged\n",
      "  )\n",
      "}\n\n",
      "# Summary of split analysis\n",
      "cat('\\n========================================\\n')\n",
      "cat('SPLIT ANALYSIS SUMMARY\\n')\n",
      "cat('========================================\\n')\n",
      "for (level in names(split_results)) {\n",
      "  res <- split_results[[level]]\n",
      "  cat('\\nGroup:', level, '\\n')\n",
      "  cat('  Original arcs:', narcs(res$net), '\\n')\n",
      "  cat('  Averaged arcs:', narcs(res$avg_net), '\\n')\n",
      "  cat('  Significant edges:', nrow(res$significant_edges), '\\n')\n",
      "  cat('  Original score:', round(res$score_original, 2), '\\n')\n",
      "  cat('  Averaged score:', round(res$score_averaged, 2), '\\n')\n",
      "}\n"
    )
  }

  # Add save results section
  code <- paste0(code,
    "\n# Save results\n",
    "# save(list = ls(), file = 'dag_analysis_results.RData')\n",
    "# write.csv(significant_edges, 'significant_edges.csv', row.names = FALSE)\n",
    "# write.csv(arc_strengths, 'arc_strengths.csv', row.names = FALSE)\n"
  )

  return(code)
}

# MODIFIED: Create constraints panel with multi-select capability for both FROM and TO
create_constraints_panel <- function() {
  wellPanel(
    h4("🚫 Network Constraints", style = "color: #1976d2;"),

    div(class = "info-box",
        "Control edge learning: Blacklist prevents edges, Whitelist forces edges.",
        tags$br(),
        tags$strong("NEW: You can now select multiple variables in both 'From:' and 'To:' dropdowns!")),

    # Blacklist section
    h5("🚫 Blacklisted Edges", style = "color: #d32f2f;"),
    p("Prevent these edges from being learned:", style = "font-size: 0.9em; color: #666;"),

    div(id = "blacklist-container",
        div(id = "blacklist-row-1", class = "row blacklist-row", style = "margin-bottom: 5px;",
            div(class = "col-sm-5",
                selectInput("bl_from_1", "From (multi-select):",
                           choices = NULL,
                           width = "100%",
                           multiple = TRUE)),
            div(class = "col-sm-5",
                selectInput("bl_to_1", "To (multi-select):",
                           choices = NULL,
                           width = "100%",
                           multiple = TRUE)),
            div(class = "col-sm-2",
                actionButton("remove_bl_1", "x", class = "btn-sm btn-outline-danger", style = "margin-top: 25px;"))
        )
    ),

    div(style = "margin: 10px 0;",
        actionButton("add_blacklist", "+ Add Blacklist", class = "btn-sm btn-outline-secondary")),

    hr(),

    # Whitelist section
    h5("Whitelisted Edges", style = "color: #388e3c;"),
    p("Force these edges to be included:", style = "font-size: 0.9em; color: #666;"),

    div(id = "whitelist-container",
        div(id = "whitelist-row-1", class = "row whitelist-row", style = "margin-bottom: 5px;",
            div(class = "col-sm-5",
                selectInput("wl_from_1", "From (multi-select):",
                           choices = NULL,
                           width = "100%",
                           multiple = TRUE)),
            div(class = "col-sm-5",
                selectInput("wl_to_1", "To (multi-select):",
                           choices = NULL,
                           width = "100%",
                           multiple = TRUE)),
            div(class = "col-sm-2",
                actionButton("remove_wl_1", "x", class = "btn-sm btn-outline-danger", style = "margin-top: 25px;"))
        )
    ),

    div(style = "margin: 10px 0;",
        actionButton("add_whitelist", "+ Add Whitelist", class = "btn-sm btn-outline-secondary")),

    br(),
    div(class = "info-box", style = "background-color: #fff3cd;",
        "Tip: Use blacklists to prevent impossible relationships and whitelists to ensure known causal relationships are included.",
        tags$br(),
        "Multi-select: Hold Ctrl/Cmd while clicking to select multiple variables in both dropdowns.",
        tags$br(),
        "This creates all combinations: e.g., From: [A,B] x To: [C,D] creates A->C, A->D, B->C, B->D")
  )
}

# Create split analysis panel
create_split_panel <- function() {
  wellPanel(
    h4("Split Analysis", style = "color: #1976d2;"),

    div(class = "info-box",
        "Perform separate DAG analyses for subgroups defined by a splitting variable."),

    checkboxInput("enableSplit", "Enable Split Analysis", FALSE),

    conditionalPanel(
      condition = "input.enableSplit",
      selectInput("splitVar", "Split Variable:", choices = NULL),
      div(class = "info-box", style = "background-color: #e8f5e8;",
          "The split variable should be categorical with 2-10 groups. Networks will be learned separately for each group."),

      h5("Comparison Options:", style = "color: #1976d2; margin-top: 15px;"),
      checkboxInput("compareSideBySide", "Show networks side-by-side", TRUE),
      checkboxInput("showGroupStats", "Show group statistics", TRUE),
      checkboxInput("synchronizeLayouts", "Synchronize node positions across groups", TRUE),
      div(class = "info-box", style = "background-color: #e1f5fe; font-size: 0.85em;",
          "When synchronized, nodes start at the same positions but can still be dragged. Uncheck to allow independent layouts per group."),

      conditionalPanel(
        condition = "input.synchronizeLayouts && input.plotType == 'interactive'",
        checkboxInput("enablePhysicsForFixed", "Enable physics for circle/cascade layouts", FALSE),
        div(class = "info-box", style = "background-color: #fff9c4; font-size: 0.85em;",
            "By default, physics is disabled for circle/cascade layouts to keep nodes stable. Enable this to allow spring-like movement.")
      ),

      div(class = "info-box", style = "background-color: #fff3cd; margin-top: 10px;",
          "Note: Each subgroup needs sufficient data (>30 obs) for reliable network learning."),

      hr(),
      h5("DAG Difference Graph:", style = "color: #1976d2; margin-top: 10px;"),
      checkboxInput("showDiffGraph", "Show difference graph between two groups", FALSE),
      conditionalPanel(
        condition = "input.showDiffGraph",
        selectInput("diffGroup1", "Group A (red = unique to A):", choices = NULL),
        selectInput("diffGroup2", "Group B (green = unique to B):", choices = NULL),
        div(class = "info-box", style = "font-size: 0.85em; margin-bottom: 6px;",
            "Toggle edge categories:"),
        checkboxInput("diffShowShared", HTML("<span style='color:#5C8DB8; font-weight:bold;'>&#9632;</span> Shared by both"), TRUE),
        checkboxInput("diffShowOnly1",  HTML("<span style='color:#E53935; font-weight:bold;'>&#9632;</span> Only in Group A"), TRUE),
        checkboxInput("diffShowOnly2",  HTML("<span style='color:#43A047; font-weight:bold;'>&#9632;</span> Only in Group B"), TRUE)
      )
    )
  )
}

# ============================================================================
# UI HELPER FUNCTIONS
# ============================================================================

# Create data input panel
create_data_panel <- function() {
  wellPanel(
    h4("Data Input", style = "color: #1976d2;"),

    # Option to start fresh or load previous work
    radioButtons("dataSource", "Start with:",
                choices = list(
                  "New Data (CSV/Excel/SPSS/Text)" = "new",
                  "Load Previous Analysis" = "analysis",
                  "Load Settings Only" = "settings"
                ),
                selected = "new",
                inline = FALSE),

    hr(),

    # New data upload
    conditionalPanel(
      condition = "input.dataSource == 'new'",
      fileInput("datafile", "Upload Data File:",
                accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls", ".sav")),
      div(class = "info-box",
          "Supported formats: CSV (.csv), Text (.txt, .tsv), Excel (.xlsx, .xls), SPSS (.sav). Max 30MB."),
      checkboxInput("dataIsLong", "Data is in long format (long-by-time)", value = FALSE),
      conditionalPanel(
        condition = "input.dataIsLong",
        div(class = "info-box", style = "font-size:0.85em;",
            "Select the subject ID and time columns.",
            tags$br(),
            "All other columns become variables renamed ", tags$code("var_<time>"),
            " (e.g. anx_0, anx_1)."),
        uiOutput("longColSelectorsUI")
      )
    ),

    # Load previous analysis
    conditionalPanel(
      condition = "input.dataSource == 'analysis'",
      fileInput("loadAnalysis", "Load Analysis (.RDS):", accept = ".rds"),
      div(class = "info-box", style = "font-size: 0.85em;",
          "Load complete analysis with results and data",
          tags$br(),
          tags$strong("To load in R console:"),
          tags$code("readRDS('file.rds')"),
          tags$br(),
          tags$small("NOT load('file.rds')")
      )
    ),

    # Load settings only
    conditionalPanel(
      condition = "input.dataSource == 'settings'",
      fileInput("loadSettings", "Load Settings (.JSON):", accept = ".json"),
      div(class = "info-box", style = "font-size: 0.85em;",
          "Load parameters only (requires data upload)"),
      hr(),
      fileInput("datafile2", "Upload Data File:",
                accept = c(".csv", ".txt", ".tsv", ".xlsx", ".xls", ".sav"))
    ),

    hr(),
    checkboxInput("showSample", "Show data preview", FALSE),
    radioButtons("missingDataMethod", "Missing Data Handling:",
      choices = list(
        "Listwise Deletion (Complete Cases)" = "listwise",
        "Available Cases (Pairwise)" = "pairwise",
        "Keep All Data (No Removal)" = "none"
      ),
      selected = "listwise"
    ),

    conditionalPanel(
      condition = "input.showSample",
      numericInput("sampleRows", "Preview rows:", value = 6, min = 1, max = 50)
    )
  )
}

# Create algorithm panel
create_algorithm_panel <- function() {
  wellPanel(
    h4("Structure Learning", style = "color: #1976d2;"),

    # Algorithm Selection
    selectInput("algorithm", "Algorithm:", choices = ALGORITHMS),
    div(class = "info-box", textOutput("algorithmInfo")),

    # Show restart/perturb only for pure score-based algorithms
    conditionalPanel(
      condition = "input.algorithm == 'hc' || input.algorithm == 'tabu'",
      div(class = "info-box", style = "background-color: #fff3cd;",
          "Score-based algorithms: Support restart and perturbation parameters"),
      numericInput("restarts", "Restarts:", 50, min = 0, max = 100),
      numericInput("perturb", "Perturbations:", 100, min = 0, max = 200)
    ),

    # Show info for hybrid algorithms
    conditionalPanel(
      condition = "input.algorithm == 'mmhc'",
      div(class = "info-box", style = "background-color: #e1f5fe;",
          "Hybrid algorithm: Combines constraint-based and score-based approaches")
    ),

    # Show info for PC Stable
    conditionalPanel(
      condition = "input.algorithm == 'pc.stable'",
      div(class = "info-box", style = "background-color: #e8f5e8;",
          "PC Stable: The gold standard for causal discovery. Order-independent results and robust to variable ordering. Excellent for identifying causal relationships.")
    ),

    # Show info for constraint-based algorithms
    conditionalPanel(
      condition = "input.algorithm == 'pc.stable' || input.algorithm == 'gs' || input.algorithm == 'iamb' || input.algorithm == 'rsmax2'",
      div(class = "info-box", style = "background-color: #f3e5f5;",
          "Constraint-based algorithms: Use statistical independence tests")
    ),

    hr(),

    # Score Selection
    h5("Scoring Criterion:", style = "color: #1976d2;"),
    selectInput("score", NULL, choices = SCORE_TYPES, selected = "bic-g"),
    div(class = "info-box", style = "font-size: 0.85em;", textOutput("scoreInfo")),

    # Data Type Information
    conditionalPanel(
      condition = "output.dataTypesDetected",
      div(class = "info-box", style = "background-color: #e8f5e8;", textOutput("dataTypeInfo"))
    ),

    # Quick Score Recommendations
    div(class = "info-box", style = "background-color: #e8f5e8;",
        h6("Quick Guide:", style = "margin: 0 0 5px 0; color: #2e7d32;"),
        tags$ul(style = "margin: 5px 0; padding-left: 15px; font-size: 0.8em;",
          tags$li("Continuous data: Use '-g' scores (bic-g, aic-g, bge)"),
          tags$li("Discrete data: Use regular scores (bic, aic, bde)"),
          tags$li("Mixed data: Use '-cg' scores (bic-cg, aic-cg)"),
          tags$li("Large datasets: BIC variants (more conservative)"),
          tags$li("Small datasets: AIC variants (more permissive)"),
          tags$li("Sparse networks: Extended BIC (ebic) variants")
        )
    ),

    # Algorithm-Score Compatibility
    conditionalPanel(
      condition = "input.algorithm == 'pc.stable' || input.algorithm == 'gs' || input.algorithm == 'iamb' || input.algorithm == 'rsmax2'",
      div(class = "info-box", style = "background-color: #fff3e0;",
          "Note: Constraint-based algorithms don't use scores during learning but scores are used for bootstrap analysis and final network evaluation.")
    ),

    # Score Equivalence Information
    conditionalPanel(
      condition = "input.score == 'bde' || input.score == 'bge'",
      div(class = "info-box", style = "background-color: #e3f2fd;",
          "Score Equivalent: This score treats some differently-oriented arcs as equivalent when they don't change the network's conditional independencies.")
    ),

    conditionalPanel(
      condition = "input.score == 'k2'",
      div(class = "info-box", style = "background-color: #fce4ec;",
          "Direction Sensitive: This score is NOT score equivalent - arc directions matter and affect the score.")
    )
  )
}

# Create bootstrap panel
create_bootstrap_panel <- function() {
  wellPanel(
    h4("Bootstrap Analysis", style = "color: #1976d2;"),
    div(class = "info-box", "Bootstrap sampling assesses the stability of learned edges."),
    numericInput("bootR", "Bootstrap Samples:", 500, min = 50, max = 10000, step = 50),
    numericInput("bootRestarts", "Bootstrap Restarts:", 5, min = 0, max = 100),
    numericInput("bootPerturb", "Bootstrap Perturbations:", 10, min = 0, max = 100),
    hr(),
    h5("Averaging Thresholds:"),
    numericInput("threshold", "Strength Threshold:", 0.85, min = 0, max = 1, step = 0.05),
    div(class = "info-box", style = "padding: 5px; margin: 5px 0; font-size: 0.8em;",
        "Setting threshold to 0 shows ALL edges. Note: Zero thresholds automatically apply cextend() to orient all edges."),
    numericInput("direction", "Direction Threshold:", 0.60, min = 0, max = 1, step = 0.05),
    div(class = "info-box", "Higher thresholds = more conservative network with fewer edges. Start with higher values for cleaner networks. Zero thresholds may create undirected edges that are automatically oriented."),
    hr(),
    h5("DAG Extension:"),
    checkboxInput("forceDirectionality", "Force Full Directionality (cextend)", FALSE),
    div(class = "info-box", style = "padding: 5px; margin: 5px 0; font-size: 0.8em;",
        "Uses cextend() to orient all undirected edges, converting PDAG to a fully directed DAG. This ensures all edges have clear directionality but the orientations may be arbitrary for some edges.")
  )
}

# Create visualization panel
create_viz_panel <- function() {
  wellPanel(
    h4("Visualization", style = "color: #1976d2;"),

    # Plot Type Selection
    radioButtons("plotType", "Plot Type:",
                choices = list(
                  "Interactive (Drag Nodes)" = "interactive",
                  "Static (qgraph)" = "static",
                  "Folded Temporal" = "folded"
                ),
                selected = "static"),
    div(class = "info-box", style = "font-size: 0.85em;",
        "Interactive: drag nodes | Static: fixed layout | Folded Temporal: collapses time slices into base variables"),

    selectInput("layout", "Layout:", choices = LAYOUTS, selected = "spring"),
    div(class = "info-box", style = "font-size: 0.85em;",
        "Spring: Force-directed layout with node physics",
        tags$br(),
        "Circle: Fixed circular arrangement (no physics in split mode)",
        tags$br(),
        "Cascade: Fixed hierarchical layout based on centrality (no physics in split mode)"),
    conditionalPanel(
      condition = "input.layout == 'cascade'",
      div(class = "info-box", style = "background-color: #fff3e0;",
          "Cascade layout requires igraph package and may conflict with bnlearn. Parent nodes (higher Katz centrality) appear at the top.")
    ),
    conditionalPanel(
      condition = "input.plotType == 'static'",
      numericInput("nodeSpacing", "Node Spacing (Repulsion):", 1, min = 0.01, max = 2, step = 0.01),
      div(class = "info-box", style = "padding: 5px; margin: 5px 0; font-size: 0.8em;",
          "Controls the repulsive force between nodes in static plots")
    ),

    hr(),
    h5("Node Styling:", style = "color: #1976d2;"),

    # --- Checkboxes first ---
    checkboxInput("scaleNodes", "Scale nodes by data means (continuous only)", TRUE),
    div(class = "info-box", style = "font-size: 0.85em;",
        "Node scaling only uses continuous variables. Discrete variables get neutral size."),
    conditionalPanel(
      condition = "input.scaleNodes",
      numericInput("minNodeSize", "Min Node Size:", 5, min = 1, max = 15),
      numericInput("maxNodeSize", "Max Node Size:", 15, min = 8, max = 30)
    ),
    checkboxInput("boldLabels", "Bold Node Labels", FALSE),

    # --- Numeric/color inputs after ---
    numericInput("nodeSize", "Base Node Size:", 10, min = 1, max = 40),
    numericInput("labelSize", "Node Label Size:", 1.0, min = 0.3, max = 2, step = 0.1),

    h5("Node Colors:", style = "color: #1976d2; margin-top: 15px;"),
    colourpicker::colourInput("nodeContinuousColor", "Continuous variable color",
                              value = "#e8f4f8",
                              showColour = "both",
                              palette = "square",
                              returnName = FALSE),
    colourpicker::colourInput("nodeDiscreteColor", "Discrete variable color",
                              value = "#fff3e0",
                              showColour = "both",
                              palette = "square",
                              returnName = FALSE),
    colourpicker::colourInput("nodeBorderColor", "Node border color",
                              value = "#000000",
                              showColour = "both",
                              palette = "square",
                              returnName = FALSE),
    numericInput("nodeBorderWidth", "Node border width:", 1.5, min = 0, max = 5, step = 0.5),
    div(class = "info-box", style = "font-size: 0.85em;",
        "Choose colors for continuous and discrete variable nodes. Alpha transparency from node scaling is preserved."),

    hr(),
    h5("Edge Styling:", style = "color: #1976d2;"),

    # --- Checkboxes first ---
    checkboxInput("showLabels", "Show Edge Labels", FALSE),
    div(class = "info-box", style = "padding: 5px; margin: 5px 0; font-size: 0.8em;",
        "Edge labels show bootstrap strengths by default. When strength threshold = 0, labels show directional probabilities instead."),
    conditionalPanel(
      condition = "input.showLabels",
      numericInput("edgeLabelSize", "Edge Label Size:", 0.7, min = 0.2, max = 1.5, step = 0.1),
      checkboxInput("boldEdgeLabels", "Bold Edge Labels", FALSE)
    ),

    # --- Numeric/select inputs after ---
    selectInput("edgeDisplayType", "Edge Width Represents:",
                choices = c("Bootstrap Strength" = "strength",
                           "Information Criterion (BIC/AIC)" = "ic",
                           "Direction Probability" = "direction",
                           "Strength x Direction" = "combined"),
                selected = "strength"),
    conditionalPanel(
      condition = "input.edgeDisplayType == 'ic'",
      radioButtons("icScoreType", "IC Score Display:",
                   choices = c("Absolute (raw BIC/AIC delta)" = "absolute",
                                "Relative (% of total network score)" = "relative"),
                   selected = "absolute", inline = FALSE)
    ),
    div(class = "info-box", style = "padding: 5px; margin: 5px 0; font-size: 0.8em;",
        "Controls what edge thickness represents:",
        tags$br(),
        "Strength: Bootstrap stability",
        tags$br(),
        "IC Score: How much the network score drops when this arc is removed (relative to the bootstrap-averaged network). Absolute = raw BIC/AIC delta; Relative = % of total network score.",
        tags$br(),
        "Direction: Directional confidence",
        tags$br(),
        "Strength x Direction: Combined measure"),
    numericInput("edgeWidth", "Edge Width Multiplier:", 0.5, min = 0.1, max = 3, step = 0.1),
    numericInput("minEdgeSize", "Minimum Edge Size to Display:", 0, min = 0, max = 1, step = 0.01),
    div(class = "info-box", style = "padding: 5px; margin: 5px 0; font-size: 0.8em;",
        "Hide edges below this strength threshold for cleaner visualization"),
    numericInput("edgeTransparency", "Edge Transparency:", 0.8, min = 0.1, max = 1, step = 0.1),
    selectInput("palette", "Color Palette:",
               choices = rownames(brewer.pal.info[brewer.pal.info$category == "seq",]),
               selected = "Blues"),
    numericInput("arrowSize", "Arrow Size:", 5, min = 0.3, max = 10, step = 0.1),

    hr(),
    h5("Folded Temporal Graph:", style = "color: #1976d2; margin-top: 15px;"),
    div(class = "info-box", style = "font-size: 0.85em; background-color: #f3e5f5;",
        "Variables need shared temporal prefixes or suffixes (e.g. NORM_VAR / LD_VAR or VAR_1 / VAR_2).",
        tags$br(), "Temporal order is inferred automatically from blacklisted edges."),
    selectInput("temporalWeightType", "Temporal Edge Width Represents:",
                choices = c("Strength x Direction (default)" = "combined",
                           "Bootstrap Strength" = "strength",
                           "Direction Probability" = "direction"),
                selected = "combined"),
    colourpicker::colourInput("temporalCrosslagColor",
                              "Cross-lag / Self-loop color",
                              value = "#ED7474",
                              showColour = "both", palette = "square", returnName = FALSE),
    uiOutput("temporalContempColorsUI"),
    colourpicker::colourInput("temporalNodeColor",
                              "Node fill color",
                              value = "#E5F2FF",
                              showColour = "both", palette = "square", returnName = FALSE),
    colourpicker::colourInput("temporalNodeBorderColor",
                              "Node border color",
                              value = "#B3CBE6",
                              showColour = "both", palette = "square", returnName = FALSE),
    checkboxInput("temporalHideContempArrows", "Hide arrowheads on contemporaneous edges", value = FALSE),
    checkboxInput("temporalContempDotted", "Dotted lines for contemporaneous edges", value = TRUE),
    hr(),
    h6("Show / hide edge types:", style = "margin-bottom:4px;"),
    checkboxInput("temporalShowCrosslagEdges", "Cross-lag / Self-loop edges", value = TRUE),
    uiOutput("temporalContempToggleUI")
  )
}

# ============================================================================
# USER INTERFACE
# ============================================================================

ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .content-wrapper { padding: 20px; }
      .well { background-color: #f8f9fa; }
      .error-msg { color: #d32f2f; font-weight: bold; }
      .success-msg { color: #388e3c; font-weight: bold; }
      .info-box { background-color: #e3f2fd; padding: 10px; margin: 10px 0; border-radius: 5px; }
      body { overflow-x: hidden; }

      /* Sidebar panel styling */
      .sidebar-panel {
        height: calc(100vh - 120px);
        overflow-y: auto;
        overflow-x: hidden;
        padding-right: 15px;
        position: fixed;
        width: 25%;
        min-width: 250px;
        max-width: 50%;
        top: 120px;
        scroll-behavior: smooth;
        border-right: 3px solid #e0e0e0;
        background: white;
        z-index: 10;
        padding-bottom: 120px !important;
      }

      .main-panel {
        margin-left: 25%;
        width: 75%;
        padding-top: 20px;
        min-height: calc(100vh - 100px);
        transition: margin-left 0.1s ease, width 0.1s ease;
      }

      .resize-handle {
        position: fixed;
        top: 120px;
        left: 25%;
        width: 6px;
        height: calc(100vh - 120px);
        background: #e0e0e0;
        cursor: col-resize;
        z-index: 20;
        transition: background-color 0.2s ease;
      }

      .resize-handle:hover {
        background: #2196f3;
      }

      .resize-handle.dragging {
        background: #1976d2;
      }

      /* Custom scrollbar styling for sidebar */
      .sidebar-panel::-webkit-scrollbar {
        width: 8px;
      }

      .sidebar-panel::-webkit-scrollbar-track {
        background: #f1f1f1;
        border-radius: 4px;
      }

      .sidebar-panel::-webkit-scrollbar-thumb {
        background: #c1c1c1;
        border-radius: 4px;
      }

      .sidebar-panel::-webkit-scrollbar-thumb:hover {
        background: #a8a8a8;
      }

      /* Collapsible sections */
      .collapse-toggle {
        cursor: pointer;
        padding: 8px 12px;
        background: linear-gradient(135deg, #f5f5f5, #e8e8e8);
        border-radius: 5px;
        margin-bottom: 10px;
        display: flex;
        justify-content: space-between;
        align-items: center;
        transition: all 0.3s ease;
      }

      .collapse-toggle:hover {
        background: linear-gradient(135deg, #e8e8e8, #d8d8d8);
      }

      .collapse-toggle.active {
        background: linear-gradient(135deg, #2196f3, #1976d2);
        color: white;
      }

      .collapse-arrow {
        transition: transform 0.3s ease;
      }

      .collapse-arrow.rotated {
        transform: rotate(90deg);
      }

      .collapsible-content {
        max-height: 0;
        overflow: hidden;
        transition: max-height 0.5s ease-out;
      }

      .collapsible-content.expanded {
        max-height: 9999px;
        transition: max-height 0.5s ease-in;
      }

      /* FIXED: Variable container layout */
      .var-selector-container {
        display: flex !important;
        flex-direction: row !important;
        gap: 10px !important;
        align-items: stretch !important;
        justify-content: stretch !important;
        width: 100% !important;
        height: 250px !important;
      }

      .rank-list-container {
        flex: 1 !important;
        height: 100% !important;
        position: relative !important;
        border: 2px solid #e0e0e0 !important;
        border-radius: 8px !important;
        padding: 0px !important;
        background: white !important;
        overflow: hidden !important;
        display: flex !important;
        flex-direction: column !important;
      }

      /* Header styling (fixed at top) */
      .var-header {
        font-size: 0.85em !important;
        font-weight: bold !important;
        color: white !important;
        text-align: center !important;
        padding: 6px !important;
        margin: 0 !important;
        border-bottom: 1px solid rgba(255,255,255,0.3) !important;
        flex-shrink: 0 !important;
        z-index: 1 !important;
      }

      .var-header-available {
        background: linear-gradient(135deg, #2196f3, #1976d2) !important;
      }

      .var-header-continuous {
        background: linear-gradient(135deg, #4caf50, #388e3c) !important;
      }

      .var-header-discrete {
        background: linear-gradient(135deg, #ff9800, #f57c00) !important;
      }

      /* CRITICAL FIX: Direct targeting of rank list structure */
      .rank-list-container > div.rank-list {
        flex: 1 !important;
        overflow-y: auto !important;
        overflow-x: hidden !important;
        padding: 6px !important;
        margin: 0 !important;
        min-height: 0 !important;
        max-height: calc(100% - 35px) !important;
        scrollbar-width: thin !important;
        scrollbar-color: #bdbdbd #f8f9fa !important;
      }

      .rank-list-container > div[class*='rank-list'] {
        flex: 1 !important;
        overflow-y: auto !important;
        overflow-x: hidden !important;
        padding: 6px !important;
        margin: 0 !important;
        min-height: 0 !important;
        max-height: calc(100% - 35px) !important;
      }

      /* Custom scrollbar for sortable areas */
      .rank-list-container > div.rank-list::-webkit-scrollbar,
      .rank-list-container > div[class*='rank-list']::-webkit-scrollbar {
        width: 6px !important;
      }

      .rank-list-container > div.rank-list::-webkit-scrollbar-track,
      .rank-list-container > div[class*='rank-list']::-webkit-scrollbar-track {
        background: #f8f9fa !important;
        border-radius: 3px !important;
      }

      .rank-list-container > div.rank-list::-webkit-scrollbar-thumb,
      .rank-list-container > div[class*='rank-list']::-webkit-scrollbar-thumb {
        background: #bdbdbd !important;
        border-radius: 3px !important;
      }

      /* Variable item styling */
      .rank-list-item {
        padding: 6px 10px !important;
        margin: 2px 0 !important;
        font-size: 0.8em !important;
        cursor: pointer !important;
        transition: all 0.2s ease !important;
        border-radius: 4px !important;
        border: 1px solid #e0e0e0 !important;
        background: #f8f9fa !important;
        white-space: nowrap !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        width: 100% !important;
        box-sizing: border-box !important;
        display: block !important;
        line-height: 1.2 !important;
        word-wrap: break-word !important;
        max-width: 100% !important;
      }

      /* Hover and selection effects */
      .rank-list-item:hover {
        background-color: #e3f2fd !important;
        border-color: #2196f3 !important;
        transform: translateY(-1px) !important;
        box-shadow: 0 2px 4px rgba(33, 150, 243, 0.3) !important;
      }

      .selected-item {
        background-color: #2196f3 !important;
        color: white !important;
        border-color: #1976d2 !important;
        box-shadow: 0 2px 6px rgba(33, 150, 243, 0.4) !important;
      }

      .ghost-item {
        opacity: 0.6 !important;
        background-color: #bbdefb !important;
        transform: rotate(1deg) scale(0.98) !important;
        box-shadow: 0 3px 6px rgba(0,0,0,0.2) !important;
      }

      /* Multi-select mode styles */
      .multi-select-mode {
        cursor: crosshair !important;
      }

      .multi-select-active {
        background: rgba(255, 152, 0, 0.05) !important;
        border-radius: 5px !important;
      }

      /* Data type indicators */
      .data-type-continuous {
        border-left: 3px solid #4caf50 !important;
      }
      .data-type-discrete {
        border-left: 3px solid #ff9800 !important;
      }

      /* Responsive design for smaller screens */
      @media (max-width: 768px) {
        .sidebar-panel {
          position: relative !important;
          height: auto !important;
          width: 100% !important;
          min-width: auto !important;
          max-width: none !important;
          top: auto !important;
          border-right: none !important;
          padding-bottom: 20px;
        }

        .main-panel {
          margin-left: 0 !important;
          width: 100% !important;
        }

        .resize-handle {
          display: none !important;
        }

        /* Stack variable containers vertically on mobile */
        .var-selector-container {
          flex-direction: column !important;
          height: auto !important;
          gap: 8px !important;
        }

        .rank-list-container {
          height: 150px !important;
          min-height: 120px !important;
        }

        .rank-list-item {
          font-size: 0.75em !important;
          padding: 5px 8px !important;
        }
      }

      /* Run button container */
      .run-button-container {
        background-color: #f8f9fa;
        border-radius: 8px;
        margin: 20px 0;
        padding: 15px;
        border: 2px solid #e9ecef;
      }

      /* Button styling */
      .btn-primary.btn-lg {
        box-shadow: 0 4px 8px rgba(0,123,255,0.3);
        transition: all 0.2s ease;
      }

      .btn-primary.btn-lg:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0,123,255,0.4);
      }

      /* Split comparison styling */
      .split-comparison-container {
        display: flex;
        flex-wrap: wrap;
        gap: 20px;
        justify-content: space-around;
      }

      .split-network-panel {
        flex: 1;
        min-width: 400px;
        border: 2px solid #e0e0e0;
        border-radius: 8px;
        padding: 10px;
        background: white;
      }

      .split-network-title {
        font-weight: bold;
        color: #1976d2;
        text-align: center;
        margin-bottom: 10px;
        padding: 5px;
        background: #f0f7ff;
        border-radius: 4px;
      }

      /* Color picker styling */
      .colourpicker-input {
        display: flex;
        align-items: center;
        gap: 10px;
        margin-bottom: 10px;
      }

      .colourpicker-input input[type='text'] {
        font-family: 'Consolas', 'Monaco', monospace;
        font-size: 13px;
        border: 1px solid #ddd;
        border-radius: 4px;
        padding: 5px;
      }

      .colourpicker-input .input-color-container {
        display: inline-block;
      }
    ")),

    # JavaScript for resizing and multi-select
    tags$script(HTML("
      $(document).ready(function() {
        // Create resize handle
        $('body').append('<div class=\"resize-handle\"></div>');

        // Panel resizing functionality
        let isResizing = false;
        let startX = 0;
        let startWidth = 0;

        // Mouse events for resizing
        $('.resize-handle').on('mousedown', function(e) {
          isResizing = true;
          startX = e.clientX;
          startWidth = $('.sidebar-panel').width();
          $(this).addClass('dragging');
          $('body').addClass('resizing');
          e.preventDefault();
        });

        $(document).on('mousemove', function(e) {
          if (!isResizing) return;

          const deltaX = e.clientX - startX;
          const newWidth = startWidth + deltaX;
          const windowWidth = $(window).width();
          const minWidth = 250;
          const maxWidth = windowWidth * 0.5;

          if (newWidth >= minWidth && newWidth <= maxWidth) {
            const widthPercent = (newWidth / windowWidth) * 100;
            const remainingPercent = 100 - widthPercent;

            $('.sidebar-panel').css('width', widthPercent + '%');
            $('.main-panel').css({
              'margin-left': widthPercent + '%',
              'width': remainingPercent + '%'
            });
            $('.resize-handle').css('left', widthPercent + '%');
          }
        });

        $(document).on('mouseup', function() {
          if (isResizing) {
            isResizing = false;
            $('.resize-handle').removeClass('dragging');
            $('body').removeClass('resizing');
          }
        });

        // Window resize handler
        $(window).on('resize', function() {
          const windowWidth = $(window).width();
          const currentWidth = $('.sidebar-panel').width();
          const maxWidth = windowWidth * 0.5;

          if (currentWidth > maxWidth) {
            $('.sidebar-panel').css('width', '50%');
            $('.main-panel').css({
              'margin-left': '50%',
              'width': '50%'
            });
            $('.resize-handle').css('left', '50%');
          }
        });

        // Keyboard shortcuts for multi-select
        $(document).on('keydown', function(e) {
          if (e.ctrlKey || e.metaKey) {
            $('.rank-list-item').addClass('multi-select-mode');
            $('.rank-list-container').addClass('multi-select-active');
          }
        });

        $(document).on('keyup', function(e) {
          if (!e.ctrlKey && !e.metaKey) {
            $('.rank-list-item').removeClass('multi-select-mode');
            $('.rank-list-container').removeClass('multi-select-active');
          }
        });
      });
    "))
  ),

  titlePanel(
    div(
      h2("The DAGger", style = "color: #1976d2; font-size: 1em; margin-bottom: 5px;"),
      p("Bayesian network learning and inference wrapper based on bnlearn (Scutari, 2007)",
        style = "color: #666; margin-top: 5px; font-size: 0.5em;")
    )
  ),

  # Use custom layout instead of sidebarLayout
  div(
    div(class = "sidebar-panel",

      create_data_panel(),

      conditionalPanel(
        condition = "output.dataLoaded",
        wellPanel(
          h4("Variables", style = "color: #1976d2;"),
          div(class = "info-box", style = "padding: 8px; margin: 5px 0 10px 0; font-size: 0.85em;",
              "Classify your variables as discrete or continuous. Drag from Available to the appropriate type container. Hold Ctrl/Cmd to select multiple items at once."),
          div(
            style = "position: relative; margin-bottom: 45px;",
            uiOutput("var_selector"),
            # Add instruction message below
            div(
              style = "position: absolute; bottom: -40px; left: 0; right: 0; text-align: center; font-size: 0.7em; color: #666; font-style: italic; background: rgba(255,255,255,0.95); padding: 6px; border-radius: 4px; border: 1px solid #e0e0e0;",
              "Hold Ctrl/Cmd + click to select multiple variables"
            )
          ),
          # Auto-detect button
          div(style = "text-align: center; margin: 10px 0;",
              actionButton("autoDetect", "Auto-Detect Data Types", class = "btn-outline-info btn-sm"))
        )
      ),

      # Collapsible sections for analysis parameters
      div(id = "algorithmSection", style = "display: none;",
          div(class = "collapse-toggle", id = "algorithmToggle",
              span("Structure Learning"),
              span(class = "collapse-arrow", "▶")),
          div(id = "algorithmContent", class = "collapsible-content",
              create_algorithm_panel())),

      div(id = "constraintsSection", style = "display: none;",
          div(class = "collapse-toggle", id = "constraintsToggle",
              span("Network Constraints"),
              span(class = "collapse-arrow", "▶")),
          div(id = "constraintsContent", class = "collapsible-content",
              create_constraints_panel())),

      div(id = "splitSection", style = "display: none;",
          div(class = "collapse-toggle", id = "splitToggle",
              span("Split Analysis"),
              span(class = "collapse-arrow", "▶")),
          div(id = "splitContent", class = "collapsible-content",
              create_split_panel())),

      div(id = "bootstrapSection", style = "display: none;",
          div(class = "collapse-toggle", id = "bootstrapToggle",
              span("Bootstrap Analysis"),
              span(class = "collapse-arrow", "▶")),
          div(id = "bootstrapContent", class = "collapsible-content",
              create_bootstrap_panel())),

      div(id = "vizSection", style = "display: none;",
          div(class = "collapse-toggle", id = "vizToggle",
              span("Visualization"),
              span(class = "collapse-arrow", "▶")),
          div(id = "vizContent", class = "collapsible-content",
              create_viz_panel())),

      div(id = "runSection", style = "display: none;",
          div(class = "run-button-container", style = "text-align: center;",
              actionButton("run", "Run Analysis", class = "btn-primary btn-lg",
                          style = "width: 100%; font-weight: bold; padding: 12px;")))
    ),

    div(class = "main-panel",
      uiOutput("statusMessages"),

      conditionalPanel(
        condition = "input.showSample && output.dataLoaded",
        wellPanel(h4("Data Preview"), withSpinner(DT::dataTableOutput("dataPreview")))
      ),

      conditionalPanel(
        condition = "output.analysisComplete",
        tabsetPanel(
          tabPanel("Network Visualization",
            wellPanel(
              h4("Learned Bayesian Network"),
              p("Edges represent learned dependencies. Edge thickness and color indicate bootstrap strength."),
              conditionalPanel(
                condition = "input.plotType == 'interactive'",
                div(class = "info-box", style = "background-color: #e8f5e8;",
                    "Interactive Mode: Click and drag nodes to reposition them. Use mouse wheel to zoom. Hover over nodes and edges for details.")
              ),
              conditionalPanel(
                condition = "input.plotType == 'interactive' && input.enableSplit && input.synchronizeLayouts",
                div(class = "info-box", style = "background-color: #fff9c4; margin-top: -10px;",
                    "Node positions are synchronized across groups for easier comparison."),
                conditionalPanel(
                  condition = "(input.layout == 'circle' || input.layout == 'cascade') && !input.enablePhysicsForFixed",
                  tags$p(style = "margin: 0; font-weight: bold; color: #e65100;",
                         "Physics disabled for stable layout. Enable in Split Analysis settings if needed.")
                )
              ),
              div(style = "text-align: center;", downloadButton("downloadPlot", "Download Plot", class = "btn-outline-primary")),
              br(), br(),
              # Show split comparison or single network based on split analysis
              conditionalPanel(
                condition = "input.enableSplit && output.splitAnalysisComplete",
                h4("Split Analysis Results", style = "color: #1976d2;"),
                uiOutput("splitComparisonUI"),
                conditionalPanel(
                  condition = "input.showDiffGraph",
                  hr(),
                  h4("DAG Difference Graph", style = "color: #1976d2;"),
                  withSpinner(visNetworkOutput("dagDiffPlot", height = "600px"))
                )
              ),
              conditionalPanel(
                condition = "!input.enableSplit || !output.splitAnalysisComplete",
                # Conditional plot outputs based on plot type
                conditionalPanel(
                  condition = "input.plotType == 'interactive'",
                  withSpinner(visNetworkOutput("interactivePlot", height = "700px"))
                ),
                conditionalPanel(
                  condition = "input.plotType == 'static'",
                  withSpinner(plotOutput("dagPlot", height = "700px"))
                ),
                conditionalPanel(
                  condition = "input.plotType == 'folded'",
                  uiOutput("temporalStructureInfo"),
                  br(),
                  withSpinner(plotOutput("foldedTemporalPlotMain", height = "700px"))
                )
              )
            )
          ),

          tabPanel("Network Statistics",
            conditionalPanel(
              condition = "input.enableSplit && output.splitAnalysisComplete",
              wellPanel(
                h4("Split Analysis Statistics"),
                div(style = "margin-bottom: 10px;",
                    downloadButton("downloadSplitStats", "Download Statistics Table", class = "btn-outline-primary")
                ),
                withSpinner(DT::dataTableOutput("splitStatistics"))
              )
            ),
            conditionalPanel(
              condition = "!input.enableSplit || !output.splitAnalysisComplete",
              fluidRow(
                column(6, wellPanel(h4("Network Properties"), withSpinner(verbatimTextOutput("networkStats")))),
                column(6, wellPanel(h4("Model Selection Scores"), withSpinner(verbatimTextOutput("modelScores"))))
              ),
              wellPanel(
                h4("Missing Data Report"),
                p("Summary of missing data and handling method used."),
                withSpinner(verbatimTextOutput("missingDataReport"))
              ),
              wellPanel(
                h4("Arc Strengths (Original Network)"),
                p("Strength scores for edges in the initially learned network."),
                div(style = "margin-bottom: 10px;",
                    downloadButton("downloadArcStrengths", "Download Arc Strengths", class = "btn-outline-primary")
                ),
                withSpinner(DT::dataTableOutput("arcStrengths"))
              )
            )
          ),

          tabPanel("Bootstrap Results",
            conditionalPanel(
              condition = "input.enableSplit && output.splitAnalysisComplete",
              wellPanel(
                h4("Split Bootstrap Analysis"),
                div(style = "margin-bottom: 15px;",
                    downloadButton("downloadAllBootstrap", "Download All Groups Combined", class = "btn-primary"),
                    " ",
                    downloadButton("downloadBootstrapSummary", "Download Summary Statistics", class = "btn-outline-primary")
                ),
                uiOutput("splitBootstrapUI")
              )
            ),
            conditionalPanel(
              condition = "!input.enableSplit || !output.splitAnalysisComplete",
              wellPanel(
                h4("Bootstrap Edge Analysis"),
                p("Edges that meet the strength and direction thresholds from bootstrap analysis."),
                div(style = "margin-bottom: 10px;",
                    downloadButton("downloadBootstrap", "Download Significant Edges", class = "btn-outline-primary"),
                    " ",
                    downloadButton("downloadFullBootstrap", "Download All Edge Strengths", class = "btn-outline-secondary")
                ),
                br(),
                withSpinner(DT::dataTableOutput("bootstrapResults"))
              ),
              wellPanel(
                h4("Bootstrap Strength Distribution"),
                p("Distribution of bootstrap strengths for all possible edges."),
                withSpinner(plotOutput("strengthHist", height = "400px"))
              )
            )
          ),

          tabPanel("Model Comparison",
            conditionalPanel(
              condition = "input.enableSplit && output.splitAnalysisComplete",
              wellPanel(
                h4("Edge Comparison Across Groups"),
                div(style = "margin-bottom: 10px;",
                    downloadButton("downloadEdgeComparison", "Download Edge Comparison", class = "btn-outline-primary")
                ),
                withSpinner(DT::dataTableOutput("edgeComparison"))
              )
            ),
            conditionalPanel(
              condition = "!input.enableSplit || !output.splitAnalysisComplete",
              wellPanel(
                h4("Network Comparison"),
                p("Compare the original learned network vs. bootstrap-averaged network."),
                div(style = "margin-bottom: 10px;",
                    downloadButton("downloadNetworkComparison", "Download Network Comparison", class = "btn-outline-primary")
                ),
                fluidRow(
                  column(6, h5("Original Network"), verbatimTextOutput("originalArcs")),
                  column(6, h5("Bootstrap-Averaged Network"), verbatimTextOutput("averagedArcs"))
                )
              )
            )
          ),

          tabPanel("Folded Temporal Graph",
            wellPanel(
              h4("Folded Temporal Network"),
              p("Temporal variables (sharing a prefix or suffix, e.g. NORM_ / LD_) are folded into their base names. ",
                "Temporal order is inferred from blacklisted edges. ",
                tags$br(),
                tags$strong("Solid edges:"), " cross-lagged (different bases, different time slices) and self-loops (same base, different time slices). ",
                tags$br(),
                tags$strong("Dotted edges:"), " contemporaneous (same time slice, different bases) - one colour per time slice."),
              uiOutput("temporalStructureInfo"),
              br(),
              div(style = "text-align: center;",
                  downloadButton("downloadTemporalPlot", "Download Temporal Plot", class = "btn-outline-primary")),
              br(),
              withSpinner(plotOutput("foldedTemporalPlot", height = "700px"))
            )
          ),

          tabPanel("R Code Export",
            wellPanel(
              h4("Standalone R Code"),
              p("Copy this code to reproduce the analysis in R without the Shiny app."),
              div(style = "margin-bottom: 15px;",
                  downloadButton("downloadRCodeTab", "Download R Script", class = "btn-outline-success")),
              div(
                style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px; padding: 10px;",
                withSpinner(verbatimTextOutput("rCodeDisplay"))
              ),
              br(),
              div(class = "info-box",
                  "This code includes all your current settings: algorithm, parameters, constraints, and visualization options. For interactive drag-and-drop visualization, consider using the visNetwork package in R.")
            )
          ),

          tabPanel("Save/Export",
            fluidRow(
              column(6,
                wellPanel(
                  h4("Save Analysis", style = "color: #1976d2;"),
                  p("Save your complete analysis results and settings."),
                  downloadButton("saveAnalysis", "Save Complete Analysis (.RDS)", class = "btn-primary btn-block"),
                  br(),
                  p("Includes all results, data, and settings", style = "font-size: 0.85em; color: #666;"),
                  div(class = "info-box", style = "background-color: #e3f2fd; margin-top: 10px;",
                      tags$strong("How to load in R:"),
                      tags$br(),
                      tags$code("data <- readRDS('your_file.rds')"),
                      tags$br(),
                      tags$small("Note: Use readRDS(), not load(), for .rds files")
                  ),
                  br(),
                  downloadButton("saveSettings", "Save Settings Only (.JSON)", class = "btn-outline-primary btn-block"),
                  br(),
                  p("Parameters only, no data/results", style = "font-size: 0.85em; color: #666;"),
                  br(),
                  div(class = "info-box", style = "font-size: 0.85em;",
                      "Share .JSON settings with colleagues to reproduce analyses with different data")
                )
              ),
              column(6,
                wellPanel(
                  h4("Export All Results", style = "color: #1976d2;"),
                  p("Export all analysis results to CSV files in a ZIP archive."),
                  downloadButton("exportAllResults", "Export All Results (.ZIP)", class = "btn-success btn-block"),
                  br(),
                  p("Includes:", style = "font-size: 0.85em; color: #666;"),
                  tags$ul(style = "font-size: 0.85em; color: #666;",
                    tags$li("Bootstrap results (significant edges)"),
                    tags$li("Full bootstrap strength matrix"),
                    tags$li("Arc strengths from original network"),
                    tags$li("Network statistics"),
                    tags$li("Edge comparison (if split analysis)"),
                    tags$li("All group-specific results (if split analysis)")
                  ),
                  hr(),
                  h5("Session Information", style = "color: #1976d2;"),
                  verbatimTextOutput("sessionInfo")
                )
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Increase maximum upload size to 30MB
  options(shiny.maxRequestSize = 30*1024^2)
  
  # Reactive values
  values <- reactiveValues(
    raw_data = NULL,
    data = NULL,
    analysis_results = NULL, 
    split_results = NULL,
    common_layout = NULL,
    error_msg = NULL, 
    success_msg = NULL,
    blacklist_count = 1,
    whitelist_count = 1,
    auto_detected = NULL,
    force_update = 0,
    loaded_analysis = NULL,
    panels_expanded = FALSE
  )
  
  # JavaScript to handle collapsible panels
  observe({
    runjs("
      // Toggle function for collapsible sections
      function toggleSection(sectionId) {
        var toggle = $('#' + sectionId + 'Toggle');
        var content = $('#' + sectionId + 'Content');
        var arrow = toggle.find('.collapse-arrow');
        
        toggle.toggleClass('active');
        content.toggleClass('expanded');
        arrow.toggleClass('rotated');
      }
      
      // Attach click handlers
      $('#algorithmToggle').off('click').on('click', function() { toggleSection('algorithm'); });
      $('#constraintsToggle').off('click').on('click', function() { toggleSection('constraints'); });
      $('#splitToggle').off('click').on('click', function() { toggleSection('split'); });
      $('#bootstrapToggle').off('click').on('click', function() { toggleSection('bootstrap'); });
      $('#vizToggle').off('click').on('click', function() { toggleSection('viz'); });
    ")
  })
  
  # Show/hide analysis panels based on variable selection
  observe({
    cont_vars <- input$continuous_vars
    disc_vars <- input$discrete_vars
    
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    
    total_vars <- length(cont_vars) + length(disc_vars)
    
    if (total_vars >= 2 || values$panels_expanded) {
      shinyjs::show("algorithmSection")
      shinyjs::show("constraintsSection")
      shinyjs::show("splitSection")
      shinyjs::show("bootstrapSection")
      shinyjs::show("vizSection")
      shinyjs::show("runSection")
      
      if (!values$panels_expanded && total_vars >= 2) {
        runjs("
          if (!$('#algorithmContent').hasClass('expanded')) {
            $('#algorithmToggle .collapse-arrow').addClass('rotated');
          }
        ")
        values$panels_expanded <- TRUE
      }
    } else {
      shinyjs::hide("algorithmSection")
      shinyjs::hide("constraintsSection")
      shinyjs::hide("splitSection")
      shinyjs::hide("bootstrapSection")
      shinyjs::hide("vizSection")
      shinyjs::hide("runSection")
      values$panels_expanded <- FALSE
    }
  })
  
  # Update split variable choices
  observe({
    req(values$data)
    
    all_vars <- get_analyzable_vars(values$data)
    
    discrete_candidates <- names(values$data)[sapply(values$data, function(x) {
      if(is.factor(x) || is.character(x)) {
        n_levels <- length(unique(x[!is.na(x)]))
        return(n_levels >= 2 && n_levels <= 10)
      } else if(is.numeric(x)) {
        n_unique <- length(unique(x[!is.na(x)]))
        return(n_unique >= 2 && n_unique <= 10)
      }
      return(FALSE)
    })]
    
    if(length(discrete_candidates) > 0) {
      updateSelectInput(session, "splitVar", choices = c("", discrete_candidates))
    } else {
      updateSelectInput(session, "splitVar", choices = "")
    }
  })
  
  # Update common layout when synchronization setting changes
  observeEvent(input$synchronizeLayouts, {
    if (!is.null(values$split_results)) {
      if (input$synchronizeLayouts) {
        cont_vars <- input$continuous_vars
        disc_vars <- input$discrete_vars
        if (is.null(cont_vars)) cont_vars <- character(0)
        if (is.null(disc_vars)) disc_vars <- character(0)
        all_vars <- unique(c(cont_vars, disc_vars))
        values$common_layout <- calculate_common_layout(all_vars, values$split_results, input$layout)
      } else {
        values$common_layout <- NULL
      }
      values$force_update <- values$force_update + 1
    }
  })
  
  # Update common layout when layout type changes
  observeEvent(input$layout, {
    if (!is.null(values$split_results) && input$synchronizeLayouts) {
      cont_vars <- input$continuous_vars
      disc_vars <- input$discrete_vars
      if (is.null(cont_vars)) cont_vars <- character(0)
      if (is.null(disc_vars)) disc_vars <- character(0)
      all_vars <- unique(c(cont_vars, disc_vars))
      values$common_layout <- calculate_common_layout(all_vars, values$split_results, input$layout)
      values$force_update <- values$force_update + 1
    }
    
    if (!is.null(values$analysis_results)) {
      values$force_update <- values$force_update + 1
    }
  })
  
  # Force layout refresh when toggle changes
  observeEvent(input$forceLayoutRefresh, {
    values$force_update <- values$force_update + 1
    
    if (!is.null(values$split_results) && input$synchronizeLayouts) {
      cont_vars <- input$continuous_vars
      disc_vars <- input$discrete_vars
      if (is.null(cont_vars)) cont_vars <- character(0)
      if (is.null(disc_vars)) disc_vars <- character(0)
      all_vars <- unique(c(cont_vars, disc_vars))
      values$common_layout <- calculate_common_layout(all_vars, values$split_results, input$layout)
    }
  })
  
  # Update physics when enable physics checkbox changes
  observeEvent(input$enablePhysicsForFixed, {
    if (!is.null(values$split_results) || !is.null(values$analysis_results)) {
      values$force_update <- values$force_update + 1
    }
  })
  
  # Algorithm and Score information
  output$algorithmInfo <- renderText({
    if (!is.null(ALGORITHM_INFO[[input$algorithm]])) {
      ALGORITHM_INFO[[input$algorithm]]
    } else {
      "Choose an algorithm to see details"
    }
  })
  
  output$scoreInfo <- renderText({
    if (!is.null(SCORE_INFO[[input$score]])) {
      SCORE_INFO[[input$score]]
    } else {
      "Choose a score to see details"
    }
  })
  
  # Data loading - handle both datafile and datafile2
  observeEvent(list(input$datafile, input$datafile2), {
    file_input <- NULL
    if (!is.null(input$datafile) && input$dataSource == "new") {
      file_input <- input$datafile
    } else if (!is.null(input$datafile2) && input$dataSource == "settings") {
      file_input <- input$datafile2
    }
    
    if (is.null(file_input)) return()
    
    tryCatch({
      df <- read_data_file(file_input$datapath, file_input$name)
      values$raw_data <- df
      values$error_msg <- NULL
    }, error = function(e) {
      values$error_msg <- paste("Error loading data:", e$message)
      values$raw_data <- NULL
      values$data     <- NULL
    })
  })
  
  # Dynamic column selectors populated from raw_data column names
  output$longColSelectorsUI <- renderUI({
    req(values$raw_data, isTRUE(input$dataIsLong))
    cols <- names(values$raw_data)
    tagList(
      selectInput("longIdCol",   "Subject / ID column:",    choices = cols, selected = cols[1]),
      selectInput("longTimeCol", "Time / occasion column:", choices = cols,
                  selected = if (length(cols) > 1) cols[2] else cols[1])
    )
  })

  # Transformation observer: converts raw_data -> data (wide or pass-through)
  observe({
    req(values$raw_data)
    df <- values$raw_data

    if (isTRUE(input$dataIsLong)) {
      id_col   <- input$longIdCol
      time_col <- input$longTimeCol
      if (is.null(id_col) || is.null(time_col) ||
          !id_col %in% names(df) || !time_col %in% names(df) ||
          id_col == time_col) return()
      tryCatch({
        wide <- long_to_wide(df, id_col, time_col)
        errors <- validate_data(wide)
        if (length(errors) > 0) {
          values$error_msg <- paste(errors, collapse = " ")
          values$data <- NULL
          return()
        }
        values$data <- wide
        values$error_msg <- NULL
        values$success_msg <- paste0(
          "Reshaped to wide: ", nrow(wide), " subjects × ", ncol(wide), " columns",
          " (", length(unique(df[[time_col]])), " time points × ",
          length(setdiff(names(df), c(id_col, time_col))), " variables)."
        )
      }, error = function(e) {
        values$error_msg <- paste("Reshape error:", e$message)
        values$data <- NULL
      })
    } else {
      errors <- validate_data(df)
      if (length(errors) > 0) {
        values$error_msg <- paste(errors, collapse = " ")
        values$data <- NULL
        return()
      }
      values$data <- df
      values$error_msg <- NULL
      values$success_msg <- paste("Successfully loaded", nrow(df), "observations with", ncol(df), "variables.")
    }
  })

  # Load previous analysis (.RDS file)
  observeEvent(input$loadAnalysis, {
    req(input$loadAnalysis)
    
    tryCatch({
      loaded_obj <- readRDS(input$loadAnalysis$datapath)
      
      if (!is.list(loaded_obj)) {
        values$error_msg <- "Invalid analysis file. Expected a saved analysis object."
        return()
      }
      
      if (is.null(loaded_obj$data)) {
        values$error_msg <- "Analysis file does not contain data."
        return()
      }
      
      values$data <- loaded_obj$data
      values$analysis_results <- loaded_obj$analysis_results
      values$split_results <- loaded_obj$split_results
      values$common_layout <- loaded_obj$common_layout
      
      if (!is.null(loaded_obj$parameters)) {
        params <- loaded_obj$parameters
        
        if (!is.null(params$algorithm)) updateSelectInput(session, "algorithm", selected = params$algorithm)
        if (!is.null(params$score)) updateSelectInput(session, "score", selected = params$score)
        if (!is.null(params$restarts)) updateNumericInput(session, "restarts", value = params$restarts)
        if (!is.null(params$perturb)) updateNumericInput(session, "perturb", value = params$perturb)
        
        if (!is.null(params$boot_r)) updateNumericInput(session, "bootR", value = params$boot_r)
        if (!is.null(params$boot_restarts)) updateNumericInput(session, "bootRestarts", value = params$boot_restarts)
        if (!is.null(params$boot_perturb)) updateNumericInput(session, "bootPerturb", value = params$boot_perturb)
        if (!is.null(params$threshold)) updateNumericInput(session, "threshold", value = params$threshold)
        if (!is.null(params$direction)) updateNumericInput(session, "direction", value = params$direction)
        if (!is.null(params$force_directionality)) updateCheckboxInput(session, "forceDirectionality", value = params$force_directionality)
        
        if (!is.null(params$continuous_vars)) {
          session$sendInputMessage("continuous_vars", list(value = params$continuous_vars))
        }
        if (!is.null(params$discrete_vars)) {
          session$sendInputMessage("discrete_vars", list(value = params$discrete_vars))
        }
        
        if (!is.null(params$split_var)) {
          updateCheckboxInput(session, "enableSplit", value = TRUE)
          updateSelectInput(session, "splitVar", selected = params$split_var)
        }
      }
      
      values$error_msg <- NULL
      values$success_msg <- "Successfully loaded previous analysis."
      values$panels_expanded <- TRUE
      
    }, error = function(e) {
      if (grepl("magic number", e$message) || grepl("bad restore file", e$message)) {
        values$error_msg <- paste0(
          "Error loading .RDS file. This file must be loaded with readRDS(), not load().\n\n",
          "Correct usage in R:\n",
          "data <- readRDS('", input$loadAnalysis$name, "')\n\n",
          "Common mistake (won't work):\n",
          "load('", input$loadAnalysis$name, "')"
        )
      } else {
        values$error_msg <- paste("Error loading analysis:", e$message)
      }
    })
  })
  
  # Load settings (.JSON file)
  observeEvent(input$loadSettings, {
    req(input$loadSettings)
    
    tryCatch({
      settings_text <- readLines(input$loadSettings$datapath)
      settings <- jsonlite::fromJSON(paste(settings_text, collapse = "\n"))
      
      if (!is.null(settings$algorithm)) updateSelectInput(session, "algorithm", selected = settings$algorithm)
      if (!is.null(settings$score)) updateSelectInput(session, "score", selected = settings$score)
      if (!is.null(settings$restarts)) updateNumericInput(session, "restarts", value = settings$restarts)
      if (!is.null(settings$perturb)) updateNumericInput(session, "perturb", value = settings$perturb)
      if (!is.null(settings$boot_r)) updateNumericInput(session, "bootR", value = settings$boot_r)
      if (!is.null(settings$boot_restarts)) updateNumericInput(session, "bootRestarts", value = settings$boot_restarts)
      if (!is.null(settings$boot_perturb)) updateNumericInput(session, "bootPerturb", value = settings$boot_perturb)
      if (!is.null(settings$threshold)) updateNumericInput(session, "threshold", value = settings$threshold)
      if (!is.null(settings$direction)) updateNumericInput(session, "direction", value = settings$direction)
      if (!is.null(settings$force_directionality)) updateCheckboxInput(session, "forceDirectionality", value = settings$force_directionality)
      
      if (!is.null(settings$visualization)) {
        viz <- settings$visualization
        if (!is.null(viz$plot_type)) updateRadioButtons(session, "plotType", selected = viz$plot_type)
        if (!is.null(viz$layout)) updateSelectInput(session, "layout", selected = viz$layout)
        if (!is.null(viz$node_size)) updateNumericInput(session, "nodeSize", value = viz$node_size)
        if (!is.null(viz$scale_nodes)) updateCheckboxInput(session, "scaleNodes", value = viz$scale_nodes)
        if (!is.null(viz$min_node_size)) updateNumericInput(session, "minNodeSize", value = viz$min_node_size)
        if (!is.null(viz$max_node_size)) updateNumericInput(session, "maxNodeSize", value = viz$max_node_size)
        if (!is.null(viz$bold_labels)) updateCheckboxInput(session, "boldLabels", value = viz$bold_labels)
        if (!is.null(viz$label_size)) updateNumericInput(session, "labelSize", value = viz$label_size)
        if (!is.null(viz$edge_width)) updateNumericInput(session, "edgeWidth", value = viz$edge_width)
        if (!is.null(viz$min_edge_size)) updateNumericInput(session, "minEdgeSize", value = viz$min_edge_size)
        if (!is.null(viz$palette)) updateSelectInput(session, "palette", selected = viz$palette)
        if (!is.null(viz$show_labels)) updateCheckboxInput(session, "showLabels", value = viz$show_labels)
if (!is.null(viz$edge_label_size)) updateNumericInput(session, "edgeLabelSize", value = viz$edge_label_size)
        if (!is.null(viz$bold_edge_labels)) updateCheckboxInput(session, "boldEdgeLabels", value = viz$bold_edge_labels)
        if (!is.null(viz$arrow_size)) updateNumericInput(session, "arrowSize", value = viz$arrow_size)
        if (!is.null(viz$node_spacing)) updateNumericInput(session, "nodeSpacing", value = viz$node_spacing)
        if (!is.null(viz$edge_transparency)) updateNumericInput(session, "edgeTransparency", value = viz$edge_transparency)
        if (!is.null(viz$clean_theme)) updateCheckboxInput(session, "cleanTheme", value = viz$clean_theme)
      }
      
      if (!is.null(settings$enable_split)) updateCheckboxInput(session, "enableSplit", value = settings$enable_split)
      if (!is.null(settings$split_var)) updateSelectInput(session, "splitVar", selected = settings$split_var)
      
      values$success_msg <- "Settings loaded successfully. Please upload your data file."
      
    }, error = function(e) {
      values$error_msg <- paste("Error loading settings:", e$message)
    })
  })
  
  # Status messages
  output$statusMessages <- renderUI({
    msgs <- div()
    if (!is.null(values$error_msg)) {
      msgs <- tagList(msgs, div(class = "alert alert-danger", role = "alert", 
                               icon("exclamation-triangle"), " ", 
                               values$error_msg))
    }
    if (!is.null(values$success_msg)) {
      msgs <- tagList(msgs, div(class = "alert alert-success", role = "alert", 
                               icon("check-circle"), " ", 
                               values$success_msg))
    }
    msgs
  })
  
  # Reactive flags
  output$dataLoaded <- reactive({ !is.null(values$data) })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  output$varsSelected <- reactive({ 
    cont_vars <- input$continuous_vars
    disc_vars <- input$discrete_vars
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    total_vars <- length(cont_vars) + length(disc_vars)
    return(total_vars >= 2)
  })
  outputOptions(output, "varsSelected", suspendWhenHidden = FALSE)
  
  output$dataTypesDetected <- reactive({
    cont_vars <- input$continuous_vars
    disc_vars <- input$discrete_vars
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    return(length(cont_vars) > 0 || length(disc_vars) > 0)
  })
  outputOptions(output, "dataTypesDetected", suspendWhenHidden = FALSE)
  
  output$analysisComplete <- reactive({ 
    !is.null(values$analysis_results) || !is.null(values$split_results) 
  })
  outputOptions(output, "analysisComplete", suspendWhenHidden = FALSE)
  
  output$splitAnalysisComplete <- reactive({ !is.null(values$split_results) })
  outputOptions(output, "splitAnalysisComplete", suspendWhenHidden = FALSE)
  
  # Data type information
  output$dataTypeInfo <- renderText({
    cont_vars <- input$continuous_vars
    disc_vars <- input$discrete_vars
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    cont_count <- length(cont_vars)
    disc_count <- length(disc_vars)
    if (cont_count > 0 && disc_count > 0) {
      paste("Mixed Data Detected:", cont_count, "continuous,", disc_count, "discrete variables. Use '-cg' scores.")
    } else if (cont_count > 0) {
      paste("Continuous Data:", cont_count, "variables. Use '-g' scores (Gaussian).")
    } else if (disc_count > 0) {
      paste("Discrete Data:", disc_count, "variables. Use regular scores (discrete).")
    } else {
      "No variables selected yet."
    }
  })
  
  # Variable selector with proper scrollable containers
  output$var_selector <- renderUI({
    req(values$data)
    analyzable_vars <- get_analyzable_vars(values$data)
    
    if (length(analyzable_vars) == 0) {
      return(div(class = "alert alert-warning", "No analyzable variables found in the data."))
    }
    
    div(
      class = "var-selector-container",
      
      div(
        class = "rank-list-container",
        div(
          class = "var-header var-header-available",
          "📊 Available Variables"
        ),
        rank_list(
          text = NULL,
          labels = analyzable_vars,
          input_id = "avail",
          class = "rank-list",
          options = sortable_options(
            multiDrag = TRUE,
            selectedClass = "selected-item",
            ghostClass = "ghost-item",
            animation = 150,
            group = list(name = "vars", put = TRUE, pull = TRUE)
          )
        )
      ),
      
      div(
        class = "rank-list-container",
        div(
          class = "var-header var-header-continuous",
          "📈 Continuous Variables"
        ),
        rank_list(
          text = NULL,
          labels = NULL,
          input_id = "continuous_vars",
          class = "rank-list",
          options = sortable_options(
            multiDrag = TRUE,
            selectedClass = "selected-item",
            ghostClass = "ghost-item",
            animation = 150,
            group = list(name = "vars", put = TRUE, pull = TRUE)
          )
        )
      ),
      
      div(
        class = "rank-list-container",
        div(
          class = "var-header var-header-discrete",
          "🔢 Discrete Variables"
        ),
        rank_list(
          text = NULL,
          labels = NULL,
          input_id = "discrete_vars",
          class = "rank-list",
          options = sortable_options(
            multiDrag = TRUE,
            selectedClass = "selected-item",
            ghostClass = "ghost-item",
            animation = 150,
            group = list(name = "vars", put = TRUE, pull = TRUE)
          )
        )
      )
    )
  })
  
  # Auto-detect button functionality
  observeEvent(input$autoDetect, {
    req(values$data)
    analyzable_vars <- get_analyzable_vars(values$data)
    if (length(analyzable_vars) == 0) return()
    detected <- auto_detect_data_types(values$data, analyzable_vars)
    session$sendInputMessage("avail", list(value = character(0)))
    session$sendInputMessage("continuous_vars", list(value = detected$continuous))
    session$sendInputMessage("discrete_vars", list(value = detected$discrete))
    showNotification(
      paste("Auto-detected:", length(detected$continuous), "continuous,", length(detected$discrete), "discrete variables"),
      type = "message",
      duration = 3
    )
  })
  
  # Data preview
  output$dataPreview <- DT::renderDataTable({
    req(values$data)
    n_rows <- min(input$sampleRows, nrow(values$data))
    DT::datatable(head(values$data, n_rows), options = list(scrollX = TRUE, pageLength = n_rows))
  })
  
  # Main analysis button
  observeEvent(input$run, {
    cont_vars <- input$continuous_vars
    disc_vars <- input$discrete_vars
    
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    
    if (length(cont_vars) + length(disc_vars) < 2) {
      values$error_msg <- "Please select at least 2 variables for analysis."
      return()
    }
    
    values$error_msg <- NULL
    values$analysis_results <- NULL
    values$split_results <- NULL
    values$common_layout <- NULL
    
    withProgress(message = 'Running Bayesian Network Analysis...', value = 0, {
      
      if (input$enableSplit && !is.null(input$splitVar) && input$splitVar != "") {
        # Perform split analysis
        split_var <- input$splitVar
        split_levels <- unique(values$data[[split_var]])
        split_levels <- split_levels[!is.na(split_levels)]
        
        if (length(split_levels) < 2) {
          values$error_msg <- "Split variable must have at least 2 levels."
          return()
        }
        
        if (length(split_levels) > 10) {
          values$error_msg <- "Split variable has too many levels (max 10)."
          return()
        }
        
        # COLLECT BLACKLIST AND WHITELIST (before the loop)
        blacklist <- NULL
        if (values$blacklist_count > 0) {
          bl_from <- character()
          bl_to <- character()
          for(j in 1:values$blacklist_count) {
            tryCatch({
              from_vals <- isolate(input[[paste0("bl_from_", j)]])
              to_vals <- isolate(input[[paste0("bl_to_", j)]])
              if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
                for (from_val in from_vals) {
                  for (to_val in to_vals) {
                    if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                      bl_from <- c(bl_from, from_val)
                      bl_to <- c(bl_to, to_val)
                    }
                  }
                }
              }
            }, error = function(e) { NULL })
          }
          if (length(bl_from) > 0) {
            blacklist <- data.frame(from = bl_from, to = bl_to, stringsAsFactors = FALSE)
          }
        }
        
        whitelist <- NULL
        if (values$whitelist_count > 0) {
          wl_from <- character()
          wl_to <- character()
          for(j in 1:values$whitelist_count) {
            tryCatch({
              from_vals <- isolate(input[[paste0("wl_from_", j)]])
              to_vals <- isolate(input[[paste0("wl_to_", j)]])
              if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
                for (from_val in from_vals) {
                  for (to_val in to_vals) {
                    if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                      wl_from <- c(wl_from, from_val)
                      wl_to <- c(wl_to, to_val)
                    }
                  }
                }
              }
            }, error = function(e) { NULL })
          }
          if (length(wl_from) > 0) {
            whitelist <- data.frame(from = wl_from, to = wl_to, stringsAsFactors = FALSE)
          }
        }
        
        split_results <- list()
        n_groups <- length(split_levels)
        
        for (i in seq_along(split_levels)) {
          level <- split_levels[i]
          incProgress(0, detail = paste('Analyzing group:', split_var, '=', level,
                                        '(', i, '/', n_groups, ')'))
          
          subgroup_data <- values$data[values$data[[split_var]] == level, ]
          
          if (nrow(subgroup_data) < 30) {
            showNotification(
              paste("Warning: Group", level, "has only", nrow(subgroup_data),
                    "observations. Results may be unreliable."),
              type = "warning",
              duration = 5
            )
          }
          
          tryCatch({
            # Prepare data for this subgroup based on missing data method
            if (!is.null(input$missingDataMethod)) {
              if (input$missingDataMethod == "listwise") {
                subgroup_analysis_data <- prepare_mixed_data(
                  subgroup_data, cont_vars, disc_vars, omit_na = TRUE
                )
              } else if (input$missingDataMethod == "pairwise") {
                subgroup_analysis_data <- prepare_mixed_data(
                  subgroup_data, cont_vars, disc_vars, omit_na = FALSE
                )
              } else {  # "none"
                subgroup_analysis_data <- prepare_mixed_data(
                  subgroup_data, cont_vars, disc_vars, omit_na = FALSE
                )
              }
            } else {
              subgroup_analysis_data <- prepare_mixed_data(
                subgroup_data, cont_vars, disc_vars, omit_na = TRUE
              )
            }
            
            # Build parameters for this subgroup
            params <- list(
              algorithm = input$algorithm,
              score = input$score,
              restarts = input$restarts,
              perturb = input$perturb,
              boot_r = input$bootR,
              boot_restarts = input$bootRestarts,
              boot_perturb = input$bootPerturb,
              threshold = input$threshold,
              direction = input$direction,
              force_directionality = input$forceDirectionality,
              blacklist = blacklist,
              whitelist = whitelist,
              continuous_vars = cont_vars,
              discrete_vars = disc_vars,
              missing_method = input$missingDataMethod
            )
            
            result <- run_analysis(subgroup_analysis_data, params)
            
            if (!is.null(result) && !is.null(result$net) && inherits(result$net, "bn")) {
              result$n_obs <- nrow(subgroup_analysis_data)
              result$group_name <- as.character(level)
              split_results[[as.character(level)]] <- result
            } else {
              showNotification(
                paste("Analysis failed for group", level, ": Invalid network produced"),
                type = "error",
                duration = 10
              )
            }
            
          }, error = function(e) {
            showNotification(
              paste("Error analyzing group", level, ":", e$message),
              type = "error",
              duration = 10
            )
          })
        }  # END OF FOR LOOP
        
        if (length(split_results) == 0) {
          values$error_msg <- "Split analysis failed for all groups. Check your data and settings."
          return()
        }
        
        if (length(split_results) < length(split_levels)) {
          failed_groups <- setdiff(as.character(split_levels), names(split_results))
          showNotification(
            paste("Analysis failed for group(s):", paste(failed_groups, collapse = ", ")),
            type = "warning",
            duration = 10
          )
        }
        
        values$split_results <- split_results
        
        if (length(split_results) > 0 && input$synchronizeLayouts) {
          all_vars <- unique(c(cont_vars, disc_vars))
          values$common_layout <- calculate_common_layout(all_vars, split_results, input$layout)
        } else {
          values$common_layout <- NULL
        }
        
        values$success_msg <- paste("Split analysis completed for", length(split_results), "groups.")
        
      } else {
        # Regular single analysis
        tryCatch({
          # Prepare data based on missing data method
          if (!is.null(input$missingDataMethod)) {
            if (input$missingDataMethod == "listwise") {
              analysis_data <- prepare_mixed_data(
                values$data, cont_vars, disc_vars, omit_na = TRUE
              )
            } else if (input$missingDataMethod == "pairwise") {
              analysis_data <- prepare_mixed_data(
                values$data, cont_vars, disc_vars, omit_na = FALSE
              )
            } else {  # "none"
              analysis_data <- prepare_mixed_data(
                values$data, cont_vars, disc_vars, omit_na = FALSE
              )
            }
          } else {
            analysis_data <- prepare_mixed_data(
              values$data, cont_vars, disc_vars, omit_na = TRUE
            )
          }
          
          # Collect blacklist constraints
          blacklist <- NULL
          if (values$blacklist_count > 0) {
            bl_from <- character()
            bl_to <- character()
            for(i in 1:values$blacklist_count) {
              tryCatch({
                from_vals <- isolate(input[[paste0("bl_from_", i)]])
                to_vals <- isolate(input[[paste0("bl_to_", i)]])
                if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
                  for (from_val in from_vals) {
                    for (to_val in to_vals) {
                      if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                        bl_from <- c(bl_from, from_val)
                        bl_to <- c(bl_to, to_val)
                      }
                    }
                  }
                }
              }, error = function(e) { NULL })
            }
            if (length(bl_from) > 0) {
              blacklist <- data.frame(from = bl_from, to = bl_to, stringsAsFactors = FALSE)
            }
          }
          
          # Collect whitelist constraints
          whitelist <- NULL
          if (values$whitelist_count > 0) {
            wl_from <- character()
            wl_to <- character()
            for(i in 1:values$whitelist_count) {
              tryCatch({
                from_vals <- isolate(input[[paste0("wl_from_", i)]])
                to_vals <- isolate(input[[paste0("wl_to_", i)]])
                if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
                  for (from_val in from_vals) {
                    for (to_val in to_vals) {
                      if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                        wl_from <- c(wl_from, from_val)
                        wl_to <- c(wl_to, to_val)
                      }
                    }
                  }
                }
              }, error = function(e) { NULL })
            }
            if (length(wl_from) > 0) {
              whitelist <- data.frame(from = wl_from, to = wl_to, stringsAsFactors = FALSE)
            }
          }
          
          # Build parameters
          params <- list(
            algorithm = input$algorithm,
            score = input$score,
            restarts = input$restarts,
            perturb = input$perturb,
            boot_r = input$bootR,
            boot_restarts = input$bootRestarts,
            boot_perturb = input$bootPerturb,
            threshold = input$threshold,
            direction = input$direction,
            force_directionality = input$forceDirectionality,
            blacklist = blacklist,
            whitelist = whitelist,
            continuous_vars = cont_vars,
            discrete_vars = disc_vars,
            missing_method = input$missingDataMethod,
            # Visualization parameters
            layout = input$layout,
            node_size = input$nodeSize,
            scale_nodes = input$scaleNodes,
            min_node_size = input$minNodeSize,
            max_node_size = input$maxNodeSize,
            bold_labels = input$boldLabels,
            label_size = input$labelSize,
            edge_width = input$edgeWidth,
            edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
            min_edge_size = input$minEdgeSize,
            palette = input$palette,
            node_continuous_color = input$nodeContinuousColor,
            node_discrete_color = input$nodeDiscreteColor,
            node_border_color = input$nodeBorderColor,
            show_labels = input$showLabels,
            edge_label_size = input$edgeLabelSize,
            bold_edge_labels = input$boldEdgeLabels,
            arrow_size = input$arrowSize,
            node_spacing = input$nodeSpacing,
            show_legends_toggle = input$showLegendsToggle,
            edge_transparency = input$edgeTransparency,
            clean_theme = input$cleanTheme,
            omit_na = (input$missingDataMethod == "listwise")
          )
          
          results <- run_analysis(analysis_data, params)
          values$analysis_results <- results
          values$success_msg <- "Analysis completed successfully!"
          
        }, error = function(e) {
          values$error_msg <- paste("Analysis error:", e$message)
        })
      }
    })
  })
  
  # Update constraint variable choices when selected variables change
  observe({
    cont_vars <- input$continuous_vars
    disc_vars <- input$discrete_vars
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    all_vars <- c(cont_vars, disc_vars)
    if (length(all_vars) == 0) return()
    choices <- c("", all_vars)
    
    for(i in 1:values$blacklist_count) {
      from_id <- paste0("bl_from_", i)
      to_id <- paste0("bl_to_", i)
      current_from <- isolate(input[[from_id]])
      current_to <- isolate(input[[to_id]])
      if (!is.null(current_from)) {
        valid_from_selections <- current_from[current_from %in% choices]
        updateSelectInput(session, from_id, choices = choices,
                         selected = if(length(valid_from_selections) > 0) valid_from_selections else NULL)
      } else {
        updateSelectInput(session, from_id, choices = choices, selected = NULL)
      }
      if (!is.null(current_to)) {
        valid_to_selections <- current_to[current_to %in% choices]
        updateSelectInput(session, to_id, choices = choices,
                         selected = if(length(valid_to_selections) > 0) valid_to_selections else NULL)
      } else {
        updateSelectInput(session, to_id, choices = choices, selected = NULL)
      }
    }
    
    for(i in 1:values$whitelist_count) {
      from_id <- paste0("wl_from_", i)
      to_id <- paste0("wl_to_", i)
      current_from <- isolate(input[[from_id]])
      current_to <- isolate(input[[to_id]])
      if (!is.null(current_from)) {
        valid_from_selections <- current_from[current_from %in% choices]
        updateSelectInput(session, from_id, choices = choices,
                         selected = if(length(valid_from_selections) > 0) valid_from_selections else NULL)
      } else {
        updateSelectInput(session, from_id, choices = choices, selected = NULL)
      }
      if (!is.null(current_to)) {
        valid_to_selections <- current_to[current_to %in% choices]
        updateSelectInput(session, to_id, choices = choices,
                         selected = if(length(valid_to_selections) > 0) valid_to_selections else NULL)
      } else {
        updateSelectInput(session, to_id, choices = choices, selected = NULL)
      }
    }
  })
  
  # Add blacklist row with multi-select
  observeEvent(input$add_blacklist, {
    values$blacklist_count <- values$blacklist_count + 1
    i <- values$blacklist_count
    cont_vars <- isolate(input$continuous_vars)
    disc_vars <- isolate(input$discrete_vars)
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    choices <- c("", c(cont_vars, disc_vars))
    insertUI(
      selector = "#blacklist-container",
      where = "beforeEnd",
      ui = div(id = paste0("blacklist-row-", i), class = "row blacklist-row", style = "margin-bottom: 5px;",
               div(class = "col-sm-5",
                   selectInput(paste0("bl_from_", i), "From (multi-select):",
                              choices = choices, width = "100%", multiple = TRUE, selected = character(0))),
               div(class = "col-sm-5",
                   selectInput(paste0("bl_to_", i), "To (multi-select):",
                              choices = choices, width = "100%", multiple = TRUE, selected = character(0))),
               div(class = "col-sm-2",
                   actionButton(paste0("remove_bl_", i), "x",
                              class = "btn-sm btn-outline-danger", style = "margin-top: 25px;"))
      ),
      immediate = TRUE
    )
  })
  
  # Add whitelist row with multi-select
  observeEvent(input$add_whitelist, {
    values$whitelist_count <- values$whitelist_count + 1
    i <- values$whitelist_count
    cont_vars <- isolate(input$continuous_vars)
    disc_vars <- isolate(input$discrete_vars)
    if (is.null(cont_vars)) cont_vars <- character(0)
    if (is.null(disc_vars)) disc_vars <- character(0)
    choices <- c("", c(cont_vars, disc_vars))
    insertUI(
      selector = "#whitelist-container",
      where = "beforeEnd",
      ui = div(id = paste0("whitelist-row-", i), class = "row whitelist-row", style = "margin-bottom: 5px;",
               div(class = "col-sm-5",
                   selectInput(paste0("wl_from_", i), "From (multi-select):",
                              choices = choices, width = "100%", multiple = TRUE, selected = character(0))),
               div(class = "col-sm-5",
                   selectInput(paste0("wl_to_", i), "To (multi-select):",
                              choices = choices, width = "100%", multiple = TRUE, selected = character(0))),
               div(class = "col-sm-2",
                   actionButton(paste0("remove_wl_", i), "x",
                              class = "btn-sm btn-outline-danger", style = "margin-top: 25px;"))
      ),
      immediate = TRUE
    )
  })
  
  # Remove blacklist rows
  lapply(1:100, function(i) {
    observeEvent(input[[paste0("remove_bl_", i)]], {
      removeUI(selector = paste0("#blacklist-row-", i))
    })
  })
  
  # Remove whitelist rows
  lapply(1:100, function(i) {
    observeEvent(input[[paste0("remove_wl_", i)]], {
      removeUI(selector = paste0("#whitelist-row-", i))
    })
  })
  
  # Visualization outputs for split analysis
  output$splitComparisonUI <- renderUI({
    req(values$split_results)
    input$layout
    values$force_update
    
    if (input$plotType == "folded") {
      tabs <- lapply(names(values$split_results), function(group) {
        tabPanel(
          title = paste(input$splitVar, "=", group),
          withSpinner(plotOutput(paste0("splitFoldedTemporal_", group), height = "700px"))
        )
      })
      do.call(tabsetPanel, tabs)

    } else if (input$compareSideBySide && input$plotType == "interactive") {
      plot_outputs <- lapply(names(values$split_results), function(group) {
        div(class = "split-network-panel",
            div(class = "split-network-title", paste(input$splitVar, "=", group)),
            visNetworkOutput(paste0("splitPlot_", group), height = "500px")
        )
      })
      div(class = "split-comparison-container", plot_outputs)
      
    } else if (input$compareSideBySide && input$plotType == "static") {
      plotOutput("splitPlotsStatic", height = paste0(400 * ceiling(length(values$split_results) / 2), "px"))
      
    } else {
      tabs <- lapply(names(values$split_results), function(group) {
        tabPanel(
          title = paste(input$splitVar, "=", group),
          if (input$plotType == "interactive") {
            visNetworkOutput(paste0("splitPlotTab_", group), height = "600px")
          } else {
            plotOutput(paste0("splitPlotTabStatic_", group), height = "600px")
          }
        )
      })
      do.call(tabsetPanel, tabs)
    }
  })
  
  # Generate interactive plots for each split group
  observe({
    req(values$split_results)
    values$force_update
    input$layout
    
    lapply(names(values$split_results), function(group) {
      output[[paste0("splitPlot_", group)]] <- renderVisNetwork({
        result <- values$split_results[[group]]
        params <- list(
          layout = input$layout,
          node_size = input$nodeSize,
          scale_nodes = input$scaleNodes,
          min_node_size = input$minNodeSize,
          max_node_size = input$maxNodeSize,
          bold_labels = input$boldLabels,
          label_size = input$labelSize,
          edge_width = input$edgeWidth,
          edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
          min_edge_size = input$minEdgeSize,
          palette = input$palette,
          node_continuous_color = input$nodeContinuousColor,
          node_discrete_color = input$nodeDiscreteColor,
          node_border_color = input$nodeBorderColor,
          show_labels = input$showLabels,
          bold_edge_labels = input$boldEdgeLabels,
          arrow_size = input$arrowSize,
          edge_transparency = input$edgeTransparency,
          clean_theme = input$cleanTheme,
          threshold = input$threshold,
          enable_physics_for_fixed = input$enablePhysicsForFixed
        )

        create_interactive_network_plot(
          result$boot_selected,
          result$stats$continuous_vars,
          result$stats$discrete_vars,
          result$data,
          params,
          title_suffix = paste("(", input$splitVar, "=", group, ")"),
          common_layout = values$common_layout,
          group_name = group, original_data = result$original_data
        )
      })

      output[[paste0("splitPlotTab_", group)]] <- renderVisNetwork({
        result <- values$split_results[[group]]
        params <- list(
          layout = input$layout,
          node_size = input$nodeSize,
          scale_nodes = input$scaleNodes,
          min_node_size = input$minNodeSize,
          max_node_size = input$maxNodeSize,
          bold_labels = input$boldLabels,
          label_size = input$labelSize,
          edge_width = input$edgeWidth,
          edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
          min_edge_size = input$minEdgeSize,
          palette = input$palette,
          node_continuous_color = input$nodeContinuousColor,
          node_discrete_color = input$nodeDiscreteColor,
          node_border_color = input$nodeBorderColor,
          show_labels = input$showLabels,
          bold_edge_labels = input$boldEdgeLabels,
          arrow_size = input$arrowSize,
          edge_transparency = input$edgeTransparency,
          clean_theme = input$cleanTheme,
          threshold = input$threshold,
          enable_physics_for_fixed = input$enablePhysicsForFixed
        )

        create_interactive_network_plot(
          result$boot_selected,
          result$stats$continuous_vars,
          result$stats$discrete_vars,
          result$data,
          params,
          title_suffix = paste("(", input$splitVar, "=", group, ")"),
          common_layout = values$common_layout,
          group_name = group
        )
      })
    })
  })

  # Generate static plots for tabbed view in split analysis
  observe({
    req(values$split_results)
    values$force_update
    input$layout
    
    lapply(names(values$split_results), function(group) {
      output[[paste0("splitPlotTabStatic_", group)]] <- renderPlot({
        result <- values$split_results[[group]]
        
        params <- list(
          layout = input$layout,
          node_size = input$nodeSize,
          scale_nodes = input$scaleNodes,
          min_node_size = input$minNodeSize,
          max_node_size = input$maxNodeSize,
          bold_labels = input$boldLabels,
          label_size = input$labelSize,
          edge_width = input$edgeWidth,
          edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
          min_edge_size = input$minEdgeSize,
          palette = input$palette,
          node_continuous_color = input$nodeContinuousColor,
          node_discrete_color = input$nodeDiscreteColor,
          node_border_color = input$nodeBorderColor,
          show_labels = input$showLabels,
          edge_label_size = input$edgeLabelSize,
          bold_edge_labels = input$boldEdgeLabels,
          arrow_size = input$arrowSize,
          node_spacing = input$nodeSpacing,
          show_legends_toggle = input$showLegendsToggle,
          edge_transparency = input$edgeTransparency,
          clean_theme = input$cleanTheme,
          threshold = input$threshold,
          direction = input$direction,
          force_directionality = input$forceDirectionality
        )

        create_network_plot(
          result$boot_selected,
          result$stats$continuous_vars,
          result$stats$discrete_vars,
          result$data,
          params,
          title_suffix = paste("(", input$splitVar, "=", group, ")"),
          common_layout = values$common_layout,
          group_name = group
        )
      })
    })
  })
  
  # Static plots for split analysis
  output$splitPlotsStatic <- renderPlot({
    req(values$split_results)
    values$force_update
    input$layout
    
    n_groups <- length(values$split_results)
    n_cols <- min(2, n_groups)
    n_rows <- ceiling(n_groups / n_cols)
    
    par(mfrow = c(n_rows, n_cols), mar = c(2, 2, 3, 2))
    
    for (group in names(values$split_results)) {
      result <- values$split_results[[group]]
      
      params <- list(
        layout = input$layout,
        node_size = input$nodeSize,
        scale_nodes = input$scaleNodes,
        min_node_size = input$minNodeSize,
        max_node_size = input$maxNodeSize,
        bold_labels = input$boldLabels,
        label_size = input$labelSize,
        edge_width = input$edgeWidth,
        edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
        min_edge_size = input$minEdgeSize,
        palette = input$palette,
        node_continuous_color = input$nodeContinuousColor,
        node_discrete_color = input$nodeDiscreteColor,
        node_border_color = input$nodeBorderColor,
        show_labels = input$showLabels,
        edge_label_size = input$edgeLabelSize,
        bold_edge_labels = input$boldEdgeLabels,
        arrow_size = input$arrowSize,
        node_spacing = input$nodeSpacing,
        show_legends_toggle = FALSE,
        edge_transparency = input$edgeTransparency,
        clean_theme = input$cleanTheme,
        threshold = input$threshold,
        direction = input$direction,
        force_directionality = input$forceDirectionality
      )

      create_network_plot(
        result$boot_selected,
        result$stats$continuous_vars,
        result$stats$discrete_vars,
        result$data,
        params,
        title_suffix = paste("(", input$splitVar, "=", group, ")"),
        common_layout = values$common_layout,
        group_name = group
      )
    }
  })

  # Generate folded temporal plots for each split group
  observe({
    req(values$split_results)
    lapply(names(values$split_results), function(group) {
      local({
        grp <- group
        output[[paste0("splitFoldedTemporal_", grp)]] <- renderPlot({
          result <- values$split_results[[grp]]
          info   <- temporal_info_for_result(result)
          if (is.null(info)) {
            plot.new(); par(mar = c(2,2,2,2))
            text(0.5, 0.5,
                 paste0("No temporal structure detected for group: ", grp,
                        "\nVariables need shared temporal prefixes or suffixes."),
                 cex = 1.1, col = "firebrick", adj = 0.5)
            return(invisible(NULL))
          }
          params <- make_temporal_params(info)
          plot_folded_temporal_graph(
            boot_selected = result$boot_selected,
            vars          = info$vars,
            blacklist     = info$blacklist,
            params        = params
          )
        })
      })
    })
  })

  # Split analysis statistics
  output$splitStatistics <- DT::renderDataTable({
    req(values$split_results)
    
    stats_df <- do.call(rbind, lapply(names(values$split_results), function(group) {
      result <- values$split_results[[group]]
      n_total <- if(!is.null(result$original_data)) nrow(result$original_data) else result$n_obs
      n_complete <- result$n_obs
      pct_complete <- round((n_complete / n_total) * 100, 1)
      data.frame(
        Group = group,
        `Complete Cases` = n_complete,
        `% Complete` = pct_complete,
        `Original Arcs` = result$stats$arcs,
        `Averaged Arcs` = result$stats$avg_arcs,
        `Significant Edges` = nrow(result$boot_selected),
        `Original Score` = round(result$stats$score_orig, 2),
        `Averaged Score` = round(result$stats$score_avg, 2),
        `Bootstrap Time (s)` = round(result$stats$bootstrap_time, 1),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }))
    
    DT::datatable(stats_df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Update diffGroup1/diffGroup2 choices when split results arrive
  observeEvent(values$split_results, {
    grps <- names(values$split_results)
    if (length(grps) >= 2) {
      updateSelectInput(session, "diffGroup1", choices = grps, selected = grps[1])
      updateSelectInput(session, "diffGroup2", choices = grps, selected = grps[2])
    }
  })

  # DAG difference graph
  output$dagDiffPlot <- renderVisNetwork({
    req(values$split_results, input$showDiffGraph)
    req(input$diffGroup1, input$diffGroup2)
    req(input$diffGroup1 != input$diffGroup2)
    req(input$diffGroup1 %in% names(values$split_results))
    req(input$diffGroup2 %in% names(values$split_results))
    input$layout
    values$force_update

    params <- list(
      layout                = input$layout,
      node_size             = input$nodeSize,
      node_continuous_color = input$nodeContinuousColor,
      node_discrete_color   = input$nodeDiscreteColor,
      node_border_color     = input$nodeBorderColor,
      node_border_width     = input$nodeBorderWidth,
      bold_labels           = input$boldLabels,
      label_size            = input$labelSize,
      edge_width            = input$edgeWidth,
      show_labels           = input$showLabels,
      diff_show_shared      = isTRUE(input$diffShowShared),
      diff_show_only1       = isTRUE(input$diffShowOnly1),
      diff_show_only2       = isTRUE(input$diffShowOnly2)
    )

    create_dag_diff_plot(
      result1       = values$split_results[[input$diffGroup1]],
      result2       = values$split_results[[input$diffGroup2]],
      group1_name   = input$diffGroup1,
      group2_name   = input$diffGroup2,
      params        = params,
      common_layout = values$common_layout
    )
  })

  # Regular single analysis outputs
  output$interactivePlot <- renderVisNetwork({
    req(values$analysis_results)
    if(!is.null(values$split_results)) return(NULL)
    
    input$layout
    values$force_update
    
    params <- list(
      layout = input$layout,
      node_size = input$nodeSize,
      scale_nodes = input$scaleNodes,
      min_node_size = input$minNodeSize,
      max_node_size = input$maxNodeSize,
      bold_labels = input$boldLabels,
      label_size = input$labelSize,
      edge_width = input$edgeWidth,
      edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
      min_edge_size = input$minEdgeSize,
      palette = input$palette,
      node_continuous_color = input$nodeContinuousColor,
      node_discrete_color = input$nodeDiscreteColor,
      node_border_color = input$nodeBorderColor,
      show_labels = input$showLabels,
      bold_edge_labels = input$boldEdgeLabels,
      arrow_size = input$arrowSize,
      edge_transparency = input$edgeTransparency,
      clean_theme = input$cleanTheme,
      threshold = input$threshold,
      enable_physics_for_fixed = input$enablePhysicsForFixed
    )

    create_interactive_network_plot(
      values$analysis_results$boot_selected,
      values$analysis_results$stats$continuous_vars,
      values$analysis_results$stats$discrete_vars,
      values$analysis_results$data,
      params,
      common_layout = values$common_layout
    )
  })

  output$dagPlot <- renderPlot({
    req(values$analysis_results)
    if(!is.null(values$split_results)) return(NULL)

    input$layout
    values$force_update

    params <- list(
      layout = input$layout,
      node_size = input$nodeSize,
      scale_nodes = input$scaleNodes,
      min_node_size = input$minNodeSize,
      max_node_size = input$maxNodeSize,
      bold_labels = input$boldLabels,
      label_size = input$labelSize,
      edge_width = input$edgeWidth,
      edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
      min_edge_size = input$minEdgeSize,
      palette = input$palette,
      node_continuous_color = input$nodeContinuousColor,
      node_discrete_color = input$nodeDiscreteColor,
      node_border_color = input$nodeBorderColor,
      show_labels = input$showLabels,
      edge_label_size = input$edgeLabelSize,
      bold_edge_labels = input$boldEdgeLabels,
      arrow_size = input$arrowSize,
      node_spacing = input$nodeSpacing,
      show_legends_toggle = input$showLegendsToggle,
      edge_transparency = input$edgeTransparency,
      clean_theme = input$cleanTheme,
      threshold = input$threshold,
      direction = input$direction,
      force_directionality = input$forceDirectionality
    )

    create_network_plot(
      values$analysis_results$boot_selected,
      values$analysis_results$stats$continuous_vars,
      values$analysis_results$stats$discrete_vars,
      values$analysis_results$data,
      params
    )
  })

  # Rest of the server outputs
  output$networkStats <- renderPrint({
    if(!is.null(values$analysis_results)) {
      stats <- values$analysis_results$stats
      cat("NETWORK STRUCTURE STATISTICS\n")
      cat("=============================\n")
      cat("Nodes:", stats$nodes, "\n")
      cat("Arcs (original):", stats$arcs, "\n")
      cat("Arcs (averaged):", stats$avg_arcs, "\n")
      cat("Significant edges:", nrow(values$analysis_results$boot_selected), "\n")
      cat("\nDATA TYPES\n")
      cat("Continuous variables:", length(stats$continuous_vars), "\n")
      cat("Discrete variables:", length(stats$discrete_vars), "\n")
      cat("\nCONSTRAINTS\n")
      cat("Blacklisted edges:", stats$blacklist_count, "\n")
      cat("Whitelisted edges:", stats$whitelist_count, "\n")
      cat("\nCOMPUTATION\n")
      cat("Bootstrap time:", round(stats$bootstrap_time, 1), "seconds\n")
      cat("Bootstrap iterations:", input$bootR, "\n")
      cat("Time per iteration:", round(stats$bootstrap_time / input$bootR, 3), "seconds\n")
    } else if(!is.null(values$split_results)) {
      cat("Split Analysis Summary\n")
      cat("=====================\n")
      cat("Split variable:", input$splitVar, "\n")
      cat("Number of groups:", length(values$split_results), "\n")
      cat("Groups:", paste(names(values$split_results), collapse = ", "), "\n")
    }
  })
  
  output$modelScores <- renderPrint({
    if(!is.null(values$analysis_results)) {
      stats <- values$analysis_results$stats
      cat("MODEL SELECTION SCORES\n")
      cat("=====================\n")
      cat("Score type:", input$score, "\n")
      cat("Original network:", round(stats$score_orig, 2), "\n")
      cat("Averaged network:", round(stats$score_avg, 2), "\n")
      cat("\nThresholds Used:\n")
      cat("Strength threshold:", input$threshold, "\n")
      cat("Direction threshold:", input$direction, "\n")
    } else if(!is.null(values$split_results)) {
      cat("Split Analysis Scores\n")
      cat("=====================\n")
      for(group in names(values$split_results)) {
        result <- values$split_results[[group]]
        cat("\nGroup:", group, "\n")
        cat("  Original score:", round(result$stats$score_orig, 2), "\n")
        cat("  Averaged score:", round(result$stats$score_avg, 2), "\n")
      }
    }
  })
  
  # Missing Data Report output
  output$missingDataReport <- renderPrint({
    if(!is.null(values$analysis_results)) {
      result <- values$analysis_results
      
      if (!is.null(result$original_data)) {
        orig_data <- result$original_data
      } else {
        orig_data <- result$data
      }
      
      cat("MISSING DATA REPORT\n")
      cat("===================\n")
      
      n_total <- nrow(orig_data)
      n_complete <- nrow(na.omit(orig_data))
      n_missing_cases <- n_total - n_complete
      pct_missing_cases <- round((n_missing_cases / n_total) * 100, 1)
      
      cat("Total observations:", n_total, "\n")
      cat("Complete cases:", n_complete, "\n")
      cat("Cases with any missing:", n_missing_cases, sprintf("(%.1f%%)\n", pct_missing_cases))
      cat("\n")
      
      if (!is.null(result$stats$missing_method)) {
        cat("Missing data method used:", result$stats$missing_method, "\n")
        if (!is.null(result$stats$missing_method_requested) && 
            result$stats$missing_method_requested != result$stats$missing_method) {
          cat("Originally requested:", result$stats$missing_method_requested, "\n")
          cat("** Method was changed due to data constraints **\n")
        }
      } else {
        cat("Missing data method: listwise deletion (default)\n")
      }
      
      cat("Observations used in analysis:", nrow(result$data), "\n")
      
      cat("\nMissing Values by Variable:\n")
      cat("---------------------------\n")
      
      selected_vars <- c(result$stats$continuous_vars, result$stats$discrete_vars)
      
      for (var in selected_vars) {
        if (var %in% names(orig_data)) {
          n_missing <- sum(is.na(orig_data[[var]]))
          n_available <- sum(!is.na(orig_data[[var]]))
          pct_missing <- round((n_missing / n_total) * 100, 1)
          
          if (n_missing > 0) {
            cat(sprintf("%-15s: %3d missing (%5.1f%%), %3d available\n",
                       var, n_missing, pct_missing, n_available))
          } else {
            cat(sprintf("%-15s: No missing values\n", var))
          }
        }
      }
      
      cat("\nMissingness Pattern:\n")
      cat("-------------------\n")
      missing_pattern <- apply(is.na(orig_data[selected_vars]), 1, sum)
      pattern_table <- table(missing_pattern)
      
      for (i in 1:length(pattern_table)) {
        n_vars_missing <- as.numeric(names(pattern_table)[i])
        n_cases <- pattern_table[i]
        pct_cases <- round((n_cases / n_total) * 100, 1)
        
        if (n_vars_missing == 0) {
          cat(sprintf("Complete cases: %d (%.1f%%)\n", n_cases, pct_cases))
        } else {
          cat(sprintf("%d variable(s) missing: %d cases (%.1f%%)\n",
                     n_vars_missing, n_cases, pct_cases))
        }
      }
      
    } else if(!is.null(values$split_results)) {
      cat("MISSING DATA REPORT - SPLIT ANALYSIS\n")
      cat("====================================\n")
      
      for(group in names(values$split_results)) {
        result <- values$split_results[[group]]
        cat("\nGroup:", group, "\n")
        cat("----------------\n")
        
        if (!is.null(result$original_data)) {
          orig_data <- result$original_data
        } else {
          orig_data <- result$data
        }
        
        n_total <- nrow(orig_data)
        n_complete <- nrow(na.omit(orig_data))
        n_used <- nrow(result$data)
        
        cat("  Total observations:", n_total, "\n")
        cat("  Complete cases:", n_complete, "\n")
        cat("  Used in analysis:", n_used, "\n")
        cat("  Percent complete:", round((n_complete/n_total)*100, 1), "%\n")
      }
    }
  })

  output$arcStrengths <- DT::renderDataTable({
    if(!is.null(values$analysis_results)) {
      df <- values$analysis_results$astr
      df$strength <- round(df$strength, 4)
      DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE))
    } else {
      DT::datatable(data.frame(Message = "Arc strengths available for single analysis only"))
    }
  })
  
  output$bootstrapResults <- DT::renderDataTable({
    if(!is.null(values$analysis_results)) {
      df <- values$analysis_results$boot_selected
      df$strength <- round(df$strength, 4)
      df$direction <- round(df$direction, 4)
      DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE))
    } else {
      DT::datatable(data.frame(Message = "Bootstrap results shown in Split Bootstrap Analysis tab"))
    }
  })
  
  output$strengthHist <- renderPlot({
    if(!is.null(values$analysis_results)) {
      boot_data <- values$analysis_results$boot
      hist(boot_data$strength, breaks = 30,
           main = "Distribution of Bootstrap Edge Strengths",
           xlab = "Bootstrap Strength",
           ylab = "Frequency",
           col = "lightblue", border = "darkblue")
      abline(v = input$threshold, col = "red", lwd = 2, lty = 2)
      legend("topright", legend = c(paste("Threshold =", input$threshold)),
             col = c("red"), lty = c(2), lwd = 2)
    }
  })
  
  output$originalArcs <- renderPrint({
    if(!is.null(values$analysis_results)) {
      arcs(values$analysis_results$net)
    } else if(!is.null(values$split_results)) {
      cat("Original arcs for each group:\n")
      cat("=============================\n")
      for(group in names(values$split_results)) {
        cat("\nGroup:", group, "\n")
        print(arcs(values$split_results[[group]]$net))
      }
    }
  })
  
  output$averagedArcs <- renderPrint({
    if(!is.null(values$analysis_results)) {
      arcs(values$analysis_results$avg)
    } else if(!is.null(values$split_results)) {
      cat("Averaged arcs for each group:\n")
      cat("=============================\n")
      for(group in names(values$split_results)) {
        cat("\nGroup:", group, "\n")
        print(arcs(values$split_results[[group]]$avg))
      }
    }
  })
  
  # Split bootstrap UI with download buttons for each group
  output$splitBootstrapUI <- renderUI({
    req(values$split_results)
    
    tabs <- lapply(names(values$split_results), function(group) {
      tabPanel(
        title = paste(input$splitVar, "=", group),
        div(
          div(style = "margin-bottom: 10px;",
              downloadButton(paste0("downloadBootstrap_", gsub("[^A-Za-z0-9]", "_", group)),
                           paste("Download", group, "Results"),
                           class = "btn-sm btn-outline-primary")
          ),
          DT::dataTableOutput(paste0("splitBootstrap_", group))
        )
      )
    })
    
    do.call(tabsetPanel, tabs)
  })
  
  # Generate bootstrap tables for each group
  observe({
    req(values$split_results)
    
    lapply(names(values$split_results), function(group) {
      output[[paste0("splitBootstrap_", group)]] <- DT::renderDataTable({
        result <- values$split_results[[group]]
        df <- result$boot_selected
        if(nrow(df) > 0) {
          df$strength <- round(df$strength, 4)
          df$direction <- round(df$direction, 4)
        }
        DT::datatable(df, options = list(pageLength = 15, scrollX = TRUE))
      })
      
      local({
        current_group <- group
        output[[paste0("downloadBootstrap_", gsub("[^A-Za-z0-9]", "_", current_group))]] <- downloadHandler(
          filename = function() {
            paste0("bootstrap_", input$splitVar, "_", current_group, "_", Sys.Date(), ".csv")
          },
          content = function(file) {
            result <- values$split_results[[current_group]]
            df <- result$boot_selected
            if(nrow(df) > 0) {
              df$group <- current_group
              df$split_variable <- input$splitVar
              write.csv(df, file, row.names = FALSE)
            } else {
              write.csv(data.frame(
                from = character(),
                to = character(),
                strength = numeric(),
                direction = numeric(),
                group = character(),
                split_variable = character()
              ), file, row.names = FALSE)
            }
          }
        )
      })
    })
  })
  
  # R Code generation
  output$rCodeDisplay <- renderText({
    if(!is.null(values$analysis_results) || !is.null(values$split_results)) {
      cont_vars <- input$continuous_vars
      disc_vars <- input$discrete_vars
      if (is.null(cont_vars)) cont_vars <- character(0)
      if (is.null(disc_vars)) disc_vars <- character(0)
      
      blacklist <- NULL
      if (values$blacklist_count > 0) {
        bl_from <- character()
        bl_to <- character()
        for(i in 1:values$blacklist_count) {
          tryCatch({
            from_vals <- isolate(input[[paste0("bl_from_", i)]])
            to_vals <- isolate(input[[paste0("bl_to_", i)]])
            if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
              for (from_val in from_vals) {
                for (to_val in to_vals) {
                  if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                    bl_from <- c(bl_from, from_val)
                    bl_to <- c(bl_to, to_val)
                  }
                }
              }
            }
          }, error = function(e) { NULL })
        }
        if (length(bl_from) > 0) {
          blacklist <- data.frame(from = bl_from, to = bl_to, stringsAsFactors = FALSE)
        }
      }
      
      whitelist <- NULL
      if (values$whitelist_count > 0) {
        wl_from <- character()
        wl_to <- character()
        for(i in 1:values$whitelist_count) {
          tryCatch({
            from_vals <- isolate(input[[paste0("wl_from_", i)]])
            to_vals <- isolate(input[[paste0("wl_to_", i)]])
            if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
              for (from_val in from_vals) {
                for (to_val in to_vals) {
                  if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                    wl_from <- c(wl_from, from_val)
                    wl_to <- c(wl_to, to_val)
                  }
                }
              }
            }
          }, error = function(e) { NULL })
        }
        if (length(wl_from) > 0) {
          whitelist <- data.frame(from = wl_from, to = wl_to, stringsAsFactors = FALSE)
        }
      }
      
      params <- list(
        algorithm = input$algorithm,
        score = input$score,
        restarts = input$restarts,
        perturb = input$perturb,
        boot_r = input$bootR,
        boot_restarts = input$bootRestarts,
        boot_perturb = input$bootPerturb,
        threshold = input$threshold,
        direction = input$direction,
        force_directionality = input$forceDirectionality,
        blacklist = blacklist,
        whitelist = whitelist,
        layout = input$layout,
        node_size = input$nodeSize,
        scale_nodes = input$scaleNodes,
        min_node_size = input$minNodeSize,
        max_node_size = input$maxNodeSize,
        bold_labels = input$boldLabels,
        label_size = input$labelSize,
        edge_width = input$edgeWidth,
        edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
        palette = input$palette,
        node_continuous_color = input$nodeContinuousColor,
        node_discrete_color = input$nodeDiscreteColor,
        node_border_color = input$nodeBorderColor,
        show_labels = input$showLabels,
        arrow_size = input$arrowSize,
        omit_na = (input$missingDataMethod == "listwise"),
        missing_method = input$missingDataMethod
      )
      
      split_var <- if(input$enableSplit && !is.null(input$splitVar)) input$splitVar else NULL
      
      generate_r_code(params, cont_vars, disc_vars, "your_data.csv", split_var)
    } else {
      "No analysis results available. Please run an analysis first."
    }
  })
  
  # Download handler for R code
  output$downloadRCodeTab <- downloadHandler(
    filename = function() {
      paste0("bayesian_network_analysis_", Sys.Date(), ".R")
    },
    content = function(file) {
      if(!is.null(values$analysis_results) || !is.null(values$split_results)) {
        cont_vars <- input$continuous_vars
        disc_vars <- input$discrete_vars
        if (is.null(cont_vars)) cont_vars <- character(0)
        if (is.null(disc_vars)) disc_vars <- character(0)
        
        blacklist <- NULL
        if (values$blacklist_count > 0) {
          bl_from <- character()
          bl_to <- character()
          for(i in 1:values$blacklist_count) {
            tryCatch({
              from_vals <- isolate(input[[paste0("bl_from_", i)]])
              to_vals <- isolate(input[[paste0("bl_to_", i)]])
              if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
                for (from_val in from_vals) {
                  for (to_val in to_vals) {
                    if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                      bl_from <- c(bl_from, from_val)
                      bl_to <- c(bl_to, to_val)
                    }
                  }
                }
              }
            }, error = function(e) { NULL })
          }
          if (length(bl_from) > 0) {
            blacklist <- data.frame(from = bl_from, to = bl_to, stringsAsFactors = FALSE)
          }
        }
        
        whitelist <- NULL
        if (values$whitelist_count > 0) {
          wl_from <- character()
          wl_to <- character()
          for(i in 1:values$whitelist_count) {
            tryCatch({
              from_vals <- isolate(input[[paste0("wl_from_", i)]])
              to_vals <- isolate(input[[paste0("wl_to_", i)]])
              if (!is.null(from_vals) && length(from_vals) > 0 && !is.null(to_vals) && length(to_vals) > 0) {
                for (from_val in from_vals) {
                  for (to_val in to_vals) {
                    if (!is.null(from_val) && from_val != "" && !is.null(to_val) && to_val != "") {
                      wl_from <- c(wl_from, from_val)
                      wl_to <- c(wl_to, to_val)
                    }
                  }
                }
              }
            }, error = function(e) { NULL })
          }
          if (length(wl_from) > 0) {
            whitelist <- data.frame(from = wl_from, to = wl_to, stringsAsFactors = FALSE)
          }
        }
        
        params <- list(
          algorithm = input$algorithm,
          score = input$score,
          restarts = input$restarts,
          perturb = input$perturb,
          boot_r = input$bootR,
          boot_restarts = input$bootRestarts,
          boot_perturb = input$bootPerturb,
          threshold = input$threshold,
          direction = input$direction,
          force_directionality = input$forceDirectionality,
          blacklist = blacklist,
          whitelist = whitelist,
          layout = input$layout,
          node_size = input$nodeSize,
          scale_nodes = input$scaleNodes,
          min_node_size = input$minNodeSize,
          max_node_size = input$maxNodeSize,
          bold_labels = input$boldLabels,
          label_size = input$labelSize,
          edge_width = input$edgeWidth,
          edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
          palette = input$palette,
          node_continuous_color = input$nodeContinuousColor,
          node_discrete_color = input$nodeDiscreteColor,
          node_border_color = input$nodeBorderColor,
          show_labels = input$showLabels,
          arrow_size = input$arrowSize,
          omit_na = (input$missingDataMethod == "listwise"),
          missing_method = input$missingDataMethod
        )
        
        split_var <- if(input$enableSplit && !is.null(input$splitVar)) input$splitVar else NULL
        
        code <- generate_r_code(params, cont_vars, disc_vars, "your_data.csv", split_var)
        writeLines(code, file)
      } else {
        writeLines("# No analysis results available. Please run an analysis first.", file)
      }
    }
  )
  
  # Export all results as ZIP
  output$exportAllResults <- downloadHandler(
    filename = function() {
      if (!is.null(values$split_results)) {
        paste0("dag_analysis_split_", input$splitVar, "_", Sys.Date(), ".zip")
      } else {
        paste0("dag_analysis_results_", Sys.Date(), ".zip")
      }
    },
    content = function(file) {
      temp_dir <- tempdir()
      temp_files <- c()
      
      if (!is.null(values$split_results)) {
        all_results <- do.call(rbind, lapply(names(values$split_results), function(group) {
          result <- values$split_results[[group]]
          df <- result$boot_selected
          if(nrow(df) > 0) {
            df$group <- group
            df$split_variable <- input$splitVar
            df$n_observations <- result$n_obs
            df
          } else { NULL }
        }))
        
        if (!is.null(all_results) && nrow(all_results) > 0) {
          f <- file.path(temp_dir, "bootstrap_all_groups.csv")
          write.csv(all_results, f, row.names = FALSE)
          temp_files <- c(temp_files, f)
        }
        
        for (group in names(values$split_results)) {
          result <- values$split_results[[group]]
          if (nrow(result$boot_selected) > 0) {
            f <- file.path(temp_dir, paste0("bootstrap_", gsub("[^A-Za-z0-9]", "_", group), ".csv"))
            write.csv(result$boot_selected, f, row.names = FALSE)
            temp_files <- c(temp_files, f)
          }
        }
        
        all_edges <- unique(do.call(rbind, lapply(values$split_results, function(result) {
          if (nrow(result$boot_selected) > 0) {
            data.frame(edge = paste(result$boot_selected$from, "->", result$boot_selected$to), stringsAsFactors = FALSE)
          } else { NULL }
        })))
        
        if (!is.null(all_edges) && nrow(all_edges) > 0) {
          summary_data <- do.call(rbind, lapply(all_edges$edge, function(edge) {
            edge_parts <- strsplit(edge, " -> ")[[1]]
            from_node <- edge_parts[1]; to_node <- edge_parts[2]
            edge_data <- do.call(rbind, lapply(names(values$split_results), function(group) {
              result <- values$split_results[[group]]
              if (nrow(result$boot_selected) > 0) {
                edge_rows <- which(result$boot_selected$from == from_node & result$boot_selected$to == to_node)
                if (length(edge_rows) > 0) {
                  data.frame(group = group, strength = result$boot_selected$strength[edge_rows],
                             direction = result$boot_selected$direction[edge_rows])
                } else { NULL }
              } else { NULL }
            }))
            if (!is.null(edge_data) && nrow(edge_data) > 0) {
              data.frame(from = from_node, to = to_node,
                         n_groups_present = nrow(edge_data), total_groups = length(values$split_results),
                         consistency = nrow(edge_data) / length(values$split_results),
                         mean_strength = mean(edge_data$strength), sd_strength = sd(edge_data$strength),
                         min_strength = min(edge_data$strength), max_strength = max(edge_data$strength),
                         mean_direction = mean(edge_data$direction), sd_direction = sd(edge_data$direction),
                         groups_present = paste(edge_data$group, collapse = "; "))
            } else { NULL }
          }))
          f <- file.path(temp_dir, "bootstrap_summary.csv")
          write.csv(summary_data, f, row.names = FALSE)
          temp_files <- c(temp_files, f)
        }
        
        stats_df <- do.call(rbind, lapply(names(values$split_results), function(group) {
          result <- values$split_results[[group]]
          data.frame(Group = group, Split_Variable = input$splitVar, Observations = result$n_obs,
                     Original_Arcs = result$stats$arcs, Averaged_Arcs = result$stats$avg_arcs,
                     Significant_Edges = nrow(result$boot_selected),
                     Original_Score = round(result$stats$score_orig, 4),
                     Averaged_Score = round(result$stats$score_avg, 4),
                     Bootstrap_Time_Seconds = round(result$stats$bootstrap_time, 2), stringsAsFactors = FALSE)
        }))
        f <- file.path(temp_dir, "split_statistics.csv")
        write.csv(stats_df, f, row.names = FALSE)
        temp_files <- c(temp_files, f)
        
        for (group in names(values$split_results)) {
          result <- values$split_results[[group]]
          f <- file.path(temp_dir, paste0("full_bootstrap_", gsub("[^A-Za-z0-9]", "_", group), ".csv"))
          df <- result$boot; df$group <- group
          write.csv(df, f, row.names = FALSE)
          temp_files <- c(temp_files, f)
        }
        
      } else if (!is.null(values$analysis_results)) {
        f <- file.path(temp_dir, "bootstrap_significant_edges.csv")
        write.csv(values$analysis_results$boot_selected, f, row.names = FALSE)
        temp_files <- c(temp_files, f)
        
        f <- file.path(temp_dir, "bootstrap_full_strength.csv")
        write.csv(values$analysis_results$boot, f, row.names = FALSE)
        temp_files <- c(temp_files, f)
        
        f <- file.path(temp_dir, "arc_strengths_original.csv")
        write.csv(values$analysis_results$astr, f, row.names = FALSE)
        temp_files <- c(temp_files, f)
        
        orig_arcs <- arcs(values$analysis_results$net)
        avg_arcs <- arcs(values$analysis_results$avg)
        all_arcs <- unique(rbind(
          data.frame(from = orig_arcs[,1], to = orig_arcs[,2], stringsAsFactors = FALSE),
          data.frame(from = avg_arcs[,1], to = avg_arcs[,2], stringsAsFactors = FALSE)
        ))
        if (nrow(all_arcs) > 0) {
          all_arcs$edge <- paste(all_arcs$from, "->", all_arcs$to)
          all_arcs$in_original <- paste(all_arcs$from, all_arcs$to) %in% paste(orig_arcs[,1], orig_arcs[,2])
          all_arcs$in_averaged <- paste(all_arcs$from, all_arcs$to) %in% paste(avg_arcs[,1], avg_arcs[,2])
          for (i in 1:nrow(all_arcs)) {
            boot_idx <- which(values$analysis_results$boot$from == all_arcs$from[i] &
                            values$analysis_results$boot$to == all_arcs$to[i])
            if (length(boot_idx) > 0) {
              all_arcs$bootstrap_strength[i] <- round(values$analysis_results$boot$strength[boot_idx], 4)
              all_arcs$bootstrap_direction[i] <- round(values$analysis_results$boot$direction[boot_idx], 4)
            }
          }
          f <- file.path(temp_dir, "network_comparison.csv")
          write.csv(all_arcs, f, row.names = FALSE)
          temp_files <- c(temp_files, f)
        }
        
        stats_data <- data.frame(
          Metric = c("Nodes", "Original Arcs", "Averaged Arcs", "Significant Edges",
                    "Original Score", "Averaged Score", "Bootstrap Time (s)",
                    "Continuous Variables", "Discrete Variables",
                    "Blacklisted Edges", "Whitelisted Edges"),
          Value = c(values$analysis_results$stats$nodes,
                   values$analysis_results$stats$arcs,
                   values$analysis_results$stats$avg_arcs,
                   nrow(values$analysis_results$boot_selected),
                   round(values$analysis_results$stats$score_orig, 4),
                   round(values$analysis_results$stats$score_avg, 4),
                   round(values$analysis_results$stats$bootstrap_time, 2),
                   length(values$analysis_results$stats$continuous_vars),
                   length(values$analysis_results$stats$discrete_vars),
                   values$analysis_results$stats$blacklist_count,
                   values$analysis_results$stats$whitelist_count)
        )
        f <- file.path(temp_dir, "network_statistics.csv")
        write.csv(stats_data, f, row.names = FALSE)
        temp_files <- c(temp_files, f)
      }
      
      params_data <- data.frame(
        Parameter = c("Algorithm", "Score", "Bootstrap Iterations",
                     "Strength Threshold", "Direction Threshold",
                     "Force Directionality", "Remove Missing Data"),
        Value = c(input$algorithm, input$score, input$bootR,
                 input$threshold, input$direction,
                 input$forceDirectionality, input$omitNA)
      )
      f <- file.path(temp_dir, "analysis_parameters.csv")
      write.csv(params_data, f, row.names = FALSE)
      temp_files <- c(temp_files, f)
      
      zip(file, temp_files, flags = "-j")
    }
  )
  
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })
  
  # Download handlers
  output$downloadPlot <- downloadHandler(
    filename = function() {
      ext <- if (!is.null(input$plotType) && input$plotType == "folded") ".pdf" else ".png"
      paste0("bayesian_network_", Sys.Date(), ext)
    },
    content = function(file) {
      if (!is.null(values$analysis_results)) {
        if (!is.null(input$plotType) && input$plotType == "folded") {
          info <- temporal_info()
          params <- make_temporal_params(info)
          pdf(file, width = 10, height = 10)
          plot_folded_temporal_graph(
            boot_selected = values$analysis_results$boot_selected,
            vars          = if (!is.null(info)) info$vars else
                              c(values$analysis_results$stats$continuous_vars,
                                values$analysis_results$stats$discrete_vars),
            blacklist     = if (!is.null(info)) info$blacklist else NULL,
            params        = params
          )
          dev.off()
        } else {
          png(file, width = 1200, height = 800)
          params <- list(
            layout = input$layout,
            node_size = input$nodeSize,
            scale_nodes = input$scaleNodes,
            min_node_size = input$minNodeSize,
            max_node_size = input$maxNodeSize,
            bold_labels = input$boldLabels,
            label_size = input$labelSize,
            edge_width = input$edgeWidth,
            edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
            min_edge_size = input$minEdgeSize,
            palette = input$palette,
            node_continuous_color = input$nodeContinuousColor,
            node_discrete_color = input$nodeDiscreteColor,
            node_border_color = input$nodeBorderColor,
            show_labels = input$showLabels,
            edge_label_size = input$edgeLabelSize,
            bold_edge_labels = input$boldEdgeLabels,
            arrow_size = input$arrowSize,
            node_spacing = input$nodeSpacing,
            show_legends_toggle = input$showLegendsToggle,
            edge_transparency = input$edgeTransparency,
            clean_theme = input$cleanTheme,
            threshold = input$threshold,
            direction = input$direction,
            force_directionality = input$forceDirectionality
          )
          create_network_plot(
            values$analysis_results$boot_selected,
            values$analysis_results$stats$continuous_vars,
            values$analysis_results$stats$discrete_vars,
            values$analysis_results$data,
            params
          )
          dev.off()
        }
      }
    }
  )
  
  output$downloadBootstrap <- downloadHandler(
    filename = function() {
      paste0("bootstrap_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$analysis_results)) {
        write.csv(values$analysis_results$boot_selected, file, row.names = FALSE)
      }
    }
  )
  
  # Download all bootstrap results combined (for split analysis)
  output$downloadAllBootstrap <- downloadHandler(
    filename = function() {
      paste0("bootstrap_all_groups_", input$splitVar, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$split_results)) {
        all_results <- do.call(rbind, lapply(names(values$split_results), function(group) {
          result <- values$split_results[[group]]
          df <- result$boot_selected
          if(nrow(df) > 0) {
            df$group <- group
            df$split_variable <- input$splitVar
            df$n_observations <- result$n_obs
            df
          } else { NULL }
        }))
        if (!is.null(all_results) && nrow(all_results) > 0) {
          write.csv(all_results, file, row.names = FALSE)
        } else {
          write.csv(data.frame(from = character(), to = character(), strength = numeric(),
                               direction = numeric(), group = character(),
                               split_variable = character(), n_observations = integer()),
                   file, row.names = FALSE)
        }
      }
    }
  )
  
  # Download bootstrap summary statistics (for split analysis)
  output$downloadBootstrapSummary <- downloadHandler(
    filename = function() {
      paste0("bootstrap_summary_", input$splitVar, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$split_results)) {
        all_edges <- unique(do.call(rbind, lapply(values$split_results, function(result) {
          if (nrow(result$boot_selected) > 0) {
            data.frame(edge = paste(result$boot_selected$from, "->", result$boot_selected$to), stringsAsFactors = FALSE)
          } else { NULL }
        })))
        
        if (!is.null(all_edges) && nrow(all_edges) > 0) {
          summary_data <- do.call(rbind, lapply(all_edges$edge, function(edge) {
            edge_parts <- strsplit(edge, " -> ")[[1]]
            from_node <- edge_parts[1]; to_node <- edge_parts[2]
            edge_data <- do.call(rbind, lapply(names(values$split_results), function(group) {
              result <- values$split_results[[group]]
              if (nrow(result$boot_selected) > 0) {
                edge_rows <- which(result$boot_selected$from == from_node & result$boot_selected$to == to_node)
                if (length(edge_rows) > 0) {
                  data.frame(group = group, strength = result$boot_selected$strength[edge_rows],
                             direction = result$boot_selected$direction[edge_rows])
                } else { NULL }
              } else { NULL }
            }))
            if (!is.null(edge_data) && nrow(edge_data) > 0) {
              data.frame(from = from_node, to = to_node,
                         n_groups_present = nrow(edge_data), total_groups = length(values$split_results),
                         consistency = nrow(edge_data) / length(values$split_results),
                         mean_strength = mean(edge_data$strength), sd_strength = sd(edge_data$strength),
                         min_strength = min(edge_data$strength), max_strength = max(edge_data$strength),
                         mean_direction = mean(edge_data$direction), sd_direction = sd(edge_data$direction),
                         groups_present = paste(edge_data$group, collapse = "; "))
            } else { NULL }
          }))
          summary_data <- summary_data[order(summary_data$consistency, summary_data$mean_strength, decreasing = TRUE), ]
          write.csv(summary_data, file, row.names = FALSE)
        } else {
          write.csv(data.frame(from = character(), to = character(), n_groups_present = integer(),
                               total_groups = integer(), consistency = numeric(), mean_strength = numeric(),
                               sd_strength = numeric(), min_strength = numeric(), max_strength = numeric(),
                               mean_direction = numeric(), sd_direction = numeric(), groups_present = character()),
                   file, row.names = FALSE)
        }
      }
    }
  )
  
  # Download full bootstrap strength data
  output$downloadFullBootstrap <- downloadHandler(
    filename = function() {
      paste0("bootstrap_full_strength_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$analysis_results)) {
        write.csv(values$analysis_results$boot, file, row.names = FALSE)
      } else if (!is.null(values$split_results)) {
        all_boot <- do.call(rbind, lapply(names(values$split_results), function(group) {
          result <- values$split_results[[group]]
          df <- result$boot
          df$group <- group
          df$split_variable <- input$splitVar
          df
        }))
        write.csv(all_boot, file, row.names = FALSE)
      }
    }
  )
  
  # Download arc strengths
  output$downloadArcStrengths <- downloadHandler(
    filename = function() {
      paste0("arc_strengths_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$analysis_results)) {
        write.csv(values$analysis_results$astr, file, row.names = FALSE)
      }
    }
  )
  
  # Download split statistics
  output$downloadSplitStats <- downloadHandler(
    filename = function() {
      paste0("split_statistics_", input$splitVar, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$split_results)) {
        stats_df <- do.call(rbind, lapply(names(values$split_results), function(group) {
          result <- values$split_results[[group]]
          data.frame(Group = group, Split_Variable = input$splitVar, Observations = result$n_obs,
                     Original_Arcs = result$stats$arcs, Averaged_Arcs = result$stats$avg_arcs,
                     Significant_Edges = nrow(result$boot_selected),
                     Original_Score = round(result$stats$score_orig, 4),
                     Averaged_Score = round(result$stats$score_avg, 4),
                     Bootstrap_Time_Seconds = round(result$stats$bootstrap_time, 2),
                     Algorithm = input$algorithm, Score_Type = input$score,
                     Bootstrap_Iterations = input$bootR, Strength_Threshold = input$threshold,
                     Direction_Threshold = input$direction, stringsAsFactors = FALSE)
        }))
        write.csv(stats_df, file, row.names = FALSE)
      }
    }
  )
  
  # Download edge comparison
  output$downloadEdgeComparison <- downloadHandler(
    filename = function() {
      paste0("edge_comparison_", input$splitVar, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$split_results)) {
        all_edges <- unique(do.call(rbind, lapply(values$split_results, function(result) {
          if (nrow(result$boot_selected) > 0) {
            data.frame(edge = paste(result$boot_selected$from, "->", result$boot_selected$to), stringsAsFactors = FALSE)
          } else {
            data.frame(edge = character(0), stringsAsFactors = FALSE)
          }
        })))
        
        if (nrow(all_edges) > 0) {
          edge_comparison <- data.frame(Edge = all_edges$edge, stringsAsFactors = FALSE)
          
          for (group in names(values$split_results)) {
            result <- values$split_results[[group]]
            group_edges <- character()
            if (nrow(result$boot_selected) > 0) {
              group_edges <- paste(result$boot_selected$from, "->", result$boot_selected$to)
              group_strengths <- result$boot_selected$strength
              group_directions <- result$boot_selected$direction
              edge_comparison[[paste0(group, "_Strength")]] <- sapply(all_edges$edge, function(e) {
                idx <- which(group_edges == e)
                if (length(idx) > 0) round(group_strengths[idx], 4) else NA
              })
              edge_comparison[[paste0(group, "_Direction")]] <- sapply(all_edges$edge, function(e) {
                idx <- which(group_edges == e)
                if (length(idx) > 0) round(group_directions[idx], 4) else NA
              })
              edge_comparison[[paste0(group, "_Present")]] <- all_edges$edge %in% group_edges
            } else {
              edge_comparison[[paste0(group, "_Strength")]] <- NA
              edge_comparison[[paste0(group, "_Direction")]] <- NA
              edge_comparison[[paste0(group, "_Present")]] <- FALSE
            }
          }
          
          presence_cols <- grep("_Present$", names(edge_comparison))
          edge_comparison$Groups_Present <- rowSums(edge_comparison[, presence_cols], na.rm = TRUE)
          edge_comparison$Consistency <- edge_comparison$Groups_Present / length(values$split_results)
          edge_comparison <- edge_comparison[order(edge_comparison$Consistency, decreasing = TRUE), ]
          
          write.csv(edge_comparison, file, row.names = FALSE)
        } else {
          write.csv(data.frame(Message = "No edges found in any group"), file, row.names = FALSE)
        }
      }
    }
  )
  
  # Download network comparison
  output$downloadNetworkComparison <- downloadHandler(
    filename = function() {
      paste0("network_comparison_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(values$analysis_results)) {
        orig_arcs <- arcs(values$analysis_results$net)
        avg_arcs <- arcs(values$analysis_results$avg)
        all_arcs <- unique(rbind(
          data.frame(from = orig_arcs[,1], to = orig_arcs[,2], stringsAsFactors = FALSE),
          data.frame(from = avg_arcs[,1], to = avg_arcs[,2], stringsAsFactors = FALSE)
        ))
        if (nrow(all_arcs) > 0) {
          all_arcs$edge <- paste(all_arcs$from, "->", all_arcs$to)
          all_arcs$in_original <- paste(all_arcs$from, all_arcs$to) %in% paste(orig_arcs[,1], orig_arcs[,2])
          all_arcs$in_averaged <- paste(all_arcs$from, all_arcs$to) %in% paste(avg_arcs[,1], avg_arcs[,2])
          all_arcs$bootstrap_strength <- NA
          all_arcs$bootstrap_direction <- NA
          for (i in 1:nrow(all_arcs)) {
            boot_idx <- which(values$analysis_results$boot$from == all_arcs$from[i] &
                            values$analysis_results$boot$to == all_arcs$to[i])
            if (length(boot_idx) > 0) {
              all_arcs$bootstrap_strength[i] <- round(values$analysis_results$boot$strength[boot_idx], 4)
              all_arcs$bootstrap_direction[i] <- round(values$analysis_results$boot$direction[boot_idx], 4)
            }
          }
          all_arcs$original_strength <- NA
          for (i in 1:nrow(all_arcs)) {
            orig_idx <- which(values$analysis_results$astr$from == all_arcs$from[i] &
                            values$analysis_results$astr$to == all_arcs$to[i])
            if (length(orig_idx) > 0) {
              all_arcs$original_strength[i] <- round(values$analysis_results$astr$strength[orig_idx], 4)
            }
          }
          all_arcs <- all_arcs[order(all_arcs$bootstrap_strength, decreasing = TRUE, na.last = TRUE), ]
          write.csv(all_arcs[, c("edge", "from", "to", "in_original", "in_averaged",
                                "original_strength", "bootstrap_strength", "bootstrap_direction")],
                   file, row.names = FALSE)
        } else {
          write.csv(data.frame(Message = "No arcs found in either network"), file, row.names = FALSE)
        }
      }
    }
  )
  
  output$saveAnalysis <- downloadHandler(
    filename = function() {
      paste0("dag_analysis_", Sys.Date(), ".rds")
    },
    content = function(file) {
      save_obj <- list(
        analysis_results = values$analysis_results,
        split_results = values$split_results,
        common_layout = values$common_layout,
        parameters = list(
          algorithm = input$algorithm,
          score = input$score,
          restarts = input$restarts,
          perturb = input$perturb,
          boot_r = input$bootR,
          boot_restarts = input$bootRestarts,
          boot_perturb = input$bootPerturb,
          threshold = input$threshold,
          direction = input$direction,
          force_directionality = input$forceDirectionality,
          continuous_vars = input$continuous_vars,
          discrete_vars = input$discrete_vars,
          split_var = if(input$enableSplit) input$splitVar else NULL
        ),
        data = values$data
      )
      saveRDS(save_obj, file)
    }
  )
  
  # Save settings handler
  output$saveSettings <- downloadHandler(
    filename = function() {
      paste0("dag_settings_", Sys.Date(), ".json")
    },
    content = function(file) {
      settings <- list(
        algorithm = input$algorithm,
        score = input$score,
        restarts = input$restarts,
        perturb = input$perturb,
        boot_r = input$bootR,
        boot_restarts = input$bootRestarts,
        boot_perturb = input$bootPerturb,
        threshold = input$threshold,
        direction = input$direction,
        force_directionality = input$forceDirectionality,
        continuous_vars = input$continuous_vars,
        discrete_vars = input$discrete_vars,
        enable_split = input$enableSplit,
        split_var = input$splitVar,
        visualization = list(
          plot_type = input$plotType,
          layout = input$layout,
          node_size = input$nodeSize,
          scale_nodes = input$scaleNodes,
          min_node_size = input$minNodeSize,
          max_node_size = input$maxNodeSize,
          bold_labels = input$boldLabels,
          label_size = input$labelSize,
          edge_width = input$edgeWidth,
          edge_display_type = input$edgeDisplayType,
            ic_score_type = input$icScoreType,
          min_edge_size = input$minEdgeSize,
          palette = input$palette,
          node_continuous_color = input$nodeContinuousColor,
          node_discrete_color = input$nodeDiscreteColor,
          node_border_color = input$nodeBorderColor,
          show_labels = input$showLabels,
          edge_label_size = input$edgeLabelSize,
          bold_edge_labels = input$boldEdgeLabels,
          arrow_size = input$arrowSize,
          node_spacing = input$nodeSpacing,
          show_legends_toggle = input$showLegendsToggle,
          edge_transparency = input$edgeTransparency,
          clean_theme = input$cleanTheme
        )
      )
      write(jsonlite::toJSON(settings, pretty = TRUE), file)
    }
  )

  # -----------------------------------------------------------------------
  # FOLDED TEMPORAL GRAPH outputs
  # -----------------------------------------------------------------------

  # Helper to extract blacklist from current UI inputs
  get_current_blacklist <- reactive({
    bl_count <- values$blacklist_count
    if (bl_count < 1) return(NULL)
    bl_from <- character(); bl_to <- character()
    for (j in seq_len(bl_count)) {
      fv <- isolate(input[[paste0("bl_from_", j)]])
      tv <- isolate(input[[paste0("bl_to_", j)]])
      if (!is.null(fv) && length(fv) > 0 && !is.null(tv) && length(tv) > 0) {
        for (f in fv) for (t in tv)
          if (nzchar(f) && nzchar(t)) { bl_from <- c(bl_from, f); bl_to <- c(bl_to, t) }
      }
    }
    if (length(bl_from) > 0) data.frame(from = bl_from, to = bl_to, stringsAsFactors = FALSE) else NULL
  })

  # Non-reactive helper: compute temporal info from any result object
  temporal_info_for_result <- function(result) {
    cont_vars <- result$stats$continuous_vars
    disc_vars <- result$stats$discrete_vars
    all_vars  <- c(cont_vars, disc_vars)
    ts <- detect_temporal_structure(all_vars)
    if (is.null(ts)) return(NULL)
    bl <- get_current_blacklist()
    ordered_slices <- infer_temporal_order(ts, bl)
    list(ts = ts, ordered_slices = ordered_slices, vars = all_vars, blacklist = bl)
  }

  # Reactive: detect temporal structure from selected variables
  temporal_info <- reactive({
    if (!is.null(values$analysis_results)) {
      return(temporal_info_for_result(values$analysis_results))
    }
    if (!is.null(values$split_results) && length(values$split_results) > 0) {
      return(temporal_info_for_result(values$split_results[[1]]))
    }
    NULL
  })

  # Helper: build temporal params list from current inputs
  make_temporal_params <- function(info) {
    list(
      temporal_weight_type    = input$temporalWeightType,
      temporal_crosslag_color = input$temporalCrosslagColor,
      temporal_contemp_colors = get_contemp_colors(info),
      temporal_node_color     = input$temporalNodeColor,
      temporal_node_border    = input$temporalNodeBorderColor,
      temporal_node_size      = input$nodeSize,
      label_size              = input$labelSize,
      edge_multiplier         = input$edgeWidth,
      hide_contemp_arrows     = isTRUE(input$temporalHideContempArrows),
      arrow_size              = if (!is.null(input$arrowSize)) input$arrowSize else 4,
      show_crosslag_edges     = isTRUE(input$temporalShowCrosslagEdges),
      show_contemp_slices     = get_contemp_toggles(info),
      contemp_dotted          = isTRUE(input$temporalContempDotted),
      min_edge_weight         = if (!is.null(input$minEdgeSize)) input$minEdgeSize else 0,
      show_edge_labels        = isTRUE(input$showLabels),
      edge_label_size         = input$edgeLabelSize
    )
  }

  output$temporalStructureInfo <- renderUI({
    info <- temporal_info()
    if (is.null(info)) {
      return(div(class = "info-box", style = "background-color: #fff3e0;",
                 tags$strong("No temporal structure detected."),
                 tags$br(),
                 "Select variables that share temporal prefixes (e.g. NORM_VAR, LD_VAR) or suffixes (e.g. VAR_1, VAR_2)."))
    }
    ts <- info$ts
    div(class = "info-box", style = "background-color: #e8f5e9;",
        tags$strong(paste0("Temporal structure detected (", ts$type, ")")),
        tags$br(),
        paste0("Time slices: ", paste(info$ordered_slices, collapse = " → ")),
        tags$br(),
        paste0("Shared base variables: ", paste(head(ts$base_names, 10), collapse = ", "),
               if (length(ts$base_names) > 10) paste0(" ... (+", length(ts$base_names) - 10, " more)") else "")
    )
  })

  # Dynamic color pickers - one per detected contemporaneous time slice
  output$temporalContempColorsUI <- renderUI({
    info <- temporal_info()
    if (is.null(info) || is.null(info$ordered_slices)) return(NULL)
    default_colors <- c("#FFBF6B", "#C6DB72", "#9B59B6", "#F39C12", "#1ABC9C")
    lapply(seq_along(info$ordered_slices), function(i) {
      sl  <- info$ordered_slices[i]
      col <- default_colors[min(i, length(default_colors))]
      colourpicker::colourInput(
        inputId    = paste0("temporalContempColor_", sl),
        label      = paste0("Contemporaneous color: ", sl),
        value      = col,
        showColour = "both", palette = "square", returnName = FALSE
      )
    })
  })

  # Dynamic show/hide checkboxes - one per contemporaneous time slice
  output$temporalContempToggleUI <- renderUI({
    info <- temporal_info()
    if (is.null(info) || is.null(info$ordered_slices)) return(NULL)
    lapply(info$ordered_slices, function(sl) {
      checkboxInput(
        inputId = paste0("temporalShowContemp_", sl),
        label   = paste0("Contemporaneous: ", sl),
        value   = TRUE
      )
    })
  })

  # Helper to read per-slice visibility into a named list
  get_contemp_toggles <- function(info) {
    if (is.null(info) || is.null(info$ordered_slices)) return(list())
    setNames(
      lapply(info$ordered_slices, function(sl) {
        v <- input[[paste0("temporalShowContemp_", sl)]]
        if (is.null(v)) TRUE else isTRUE(v)
      }),
      info$ordered_slices
    )
  }

  # Helper to read per-slice contemp colors into a named list
  get_contemp_colors <- function(info) {
    if (is.null(info) || is.null(info$ordered_slices)) return(list())
    setNames(
      lapply(info$ordered_slices, function(sl) {
        v <- input[[paste0("temporalContempColor_", sl)]]
        if (!is.null(v) && nzchar(v)) v else "#4A90D9"
      }),
      info$ordered_slices
    )
  }

  # Shared render helper to avoid duplicating logic
  render_folded_temporal <- function() {
    req(values$analysis_results)
    info <- temporal_info()
    message("[FoldedPlot] analysis_results present: TRUE")
    message("[FoldedPlot] temporal_info: ", if (is.null(info)) "NULL (no temporal structure)" else
              paste0("slices=", paste(info$ordered_slices, collapse="->"),
                     "  n_edges=", nrow(values$analysis_results$boot_selected)))
    if (is.null(info)) {
      plot.new(); par(mar = c(2,2,2,2))
      text(0.5, 0.5, "No temporal structure detected.\nSelect variables with shared temporal prefixes or suffixes.",
           cex = 1.1, col = "firebrick", adj = 0.5)
      return(invisible(NULL))
    }
    params <- make_temporal_params(info)
    folded <- build_folded_edges(
      values$analysis_results$boot_selected,
      info$ts, info$ordered_slices,
      if (!is.null(params$temporal_weight_type)) params$temporal_weight_type else "combined"
    )
    message("[FoldedPlot] folded edges: ", if (is.null(folded)) "NULL" else
              paste0(nrow(folded), " rows; types: ",
                     paste(sort(unique(folded$edge_type)), collapse=", ")))
    plot_folded_temporal_graph(
      boot_selected = values$analysis_results$boot_selected,
      vars          = info$vars,
      blacklist     = info$blacklist,
      params        = params
    )
  }

  # Used in the dedicated Folded Temporal Graph tab
  output$foldedTemporalPlot <- renderPlot({ render_folded_temporal() })

  # Used inline in the Network Visualization tab
  output$foldedTemporalPlotMain <- renderPlot({ render_folded_temporal() })

  output$downloadTemporalPlot <- downloadHandler(
    filename = function() paste0("folded_temporal_graph_", Sys.Date(), ".pdf"),
    content = function(file) {
      req(values$analysis_results)
      info <- temporal_info()
      params <- make_temporal_params(info)
      pdf(file, width = 10, height = 10)
      plot_folded_temporal_graph(
        boot_selected = values$analysis_results$boot_selected,
        vars          = info$vars,
        blacklist     = info$blacklist,
        params        = params
      )
      dev.off()
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
