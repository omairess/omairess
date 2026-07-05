# ------------------------------------------------------------   
# bootnet Shiny • VISUALLY ENHANCED VERSION
# Modern UI with improved styling, colors, and UX
# ------------------------------------------------------------

# --- Install & load packages ----
install_if_missing <- function(pkgs){
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing)) install.packages(missing, dependencies = TRUE)
}
core_pkgs <- c(
  "shiny","shinyWidgets","DT","shinydashboard","colourpicker",
  "dplyr","purrr","readr","vroom","readxl","tools",
  "bootnet","qgraph","mgm","relaimpo","igraph","Matrix","EGAnet","IsingFit","networktools","huge",
  "ggplot2","gridExtra","ggraph"
)
install_if_missing(core_pkgs)
invisible(lapply(core_pkgs, require, character.only = TRUE))

# Install GitHub packages
if (!requireNamespace("NetworkComparisonTest", quietly = TRUE)) {
  message("Installing NetworkComparisonTest from GitHub...")
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("cvborkulo/NetworkComparisonTest", ref = "development")
}
library(NetworkComparisonTest)

`%||%` <- function(a,b) if (is.null(a) || (is.character(a) && identical(a,""))) b else a
is_numeric_like <- function(x) is.numeric(x) || is.integer(x)
trim_ws <- function(x) gsub("^\\s+|\\s+$","", x)

# ---------- Helper functions (unchanged from original) ----------
is_dichotomous <- function(x) {
  unique_vals <- unique(na.omit(x))
  if (length(unique_vals) != 2) return(FALSE)
  if (is.numeric(x) && all(unique_vals %in% c(0, 1))) return(TRUE)
  if (is.logical(x)) return(TRUE)
  if (is.factor(x) && nlevels(x) == 2) return(TRUE)
  if (is.character(x) && length(unique_vals) == 2) return(TRUE)
  return(FALSE)
}

dichotomous_to_numeric <- function(x) {
  if (!is_dichotomous(x)) return(x)
  unique_vals <- sort(unique(na.omit(x)))
  if (is.numeric(x) && all(unique_vals %in% c(0, 1))) return(as.numeric(x))
  if (is.logical(x)) return(as.numeric(x))
  if (is.factor(x)) return(as.numeric(x) - 1)
  result <- numeric(length(x))
  result[is.na(x)] <- NA
  result[x == unique_vals[1]] <- 0
  result[x == unique_vals[2]] <- 1
  return(result)
}

predictability_R2 <- function(data_numeric){
  if (is.null(data_numeric)) return(numeric(0))
  Xall <- as.matrix(data_numeric); storage.mode(Xall) <- "double"
  if (ncol(Xall) == 0) return(numeric(0))
  if (ncol(Xall) == 1) return(0)
  p <- ncol(Xall); r2 <- numeric(p)
  for (i in seq_len(p)){
    y <- Xall[, i]; X <- Xall[, -i, drop = FALSE]
    Xdesign <- cbind(Intercept = 1, as.matrix(X))
    fit <- tryCatch(stats::lm.fit(x = Xdesign, y = y), error = function(e) NULL)
    if (is.null(fit) || is.null(fit$residuals)) r2[i] <- 0 else {
      rss <- sum(fit$residuals^2, na.rm = TRUE)
      tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
      r2[i] <- if (tss > 0) max(0, min(1, 1 - rss/tss)) else 0
    }
  }
  as.numeric(r2)
}

compute_centrality_measures <- function(W, measures = c("Strength", "Closeness", "Betweenness")) {
  W <- as.matrix(W)
  n <- ncol(W)
  if (n == 0) return(data.frame())
  results <- data.frame(Node = if(!is.null(colnames(W))) colnames(W) else paste0("V", 1:n))
  compute_strength <- function(W) rowSums(abs(W), na.rm = TRUE)
  compute_degree <- function(W) rowSums(W != 0, na.rm = TRUE)
  compute_expected_influence <- function(W) rowSums(W, na.rm = TRUE)
  if ("Strength" %in% measures) results$Strength <- compute_strength(W)
  if ("OutStrength" %in% measures) results$OutStrength <- compute_strength(W)
  if ("InStrength" %in% measures) results$InStrength <- compute_strength(t(W))
  if ("Degree" %in% measures) results$Degree <- compute_degree(W)
  if ("OutDegree" %in% measures) results$OutDegree <- compute_degree(W)
  if ("InDegree" %in% measures) results$InDegree <- compute_degree(t(W))
  if ("ExpectedInfluence" %in% measures) results$ExpectedInfluence <- compute_expected_influence(W)
  if ("ExpectedInfluence2" %in% measures) results$ExpectedInfluence2 <- compute_expected_influence(W)
  if (any(c("Closeness", "Betweenness", "Eigenvector") %in% measures)) {
    ig <- tryCatch({
      diag(W) <- 0
      igraph::graph_from_adjacency_matrix(abs(W), mode = "undirected", weighted = TRUE)
    }, error = function(e) NULL)
    if (!is.null(ig)) {
      if ("Closeness" %in% measures) results$Closeness <- tryCatch(igraph::closeness(ig), error = function(e) rep(0, n))
      if ("Betweenness" %in% measures) results$Betweenness <- tryCatch(igraph::betweenness(ig), error = function(e) rep(0, n))
      if ("Eigenvector" %in% measures) results$Eigenvector <- tryCatch(igraph::eigen_centrality(ig)$vector, error = function(e) rep(0, n))
    } else {
      if ("Closeness" %in% measures) results$Closeness <- rep(0, n)
      if ("Betweenness" %in% measures) results$Betweenness <- rep(0, n)
      if ("Eigenvector" %in% measures) results$Eigenvector <- rep(0, n)
    }
  }
  if (any(c("Hubs", "Authorities") %in% measures)) {
    ig_dir <- tryCatch({
      diag(W) <- 0
      igraph::graph_from_adjacency_matrix(abs(W), mode = "directed", weighted = TRUE)
    }, error = function(e) NULL)
    if (!is.null(ig_dir)) {
      hub_auth <- tryCatch(igraph::hub_score(ig_dir), error = function(e) list(vector = rep(0, n)))
      auth_score <- tryCatch(igraph::authority_score(ig_dir), error = function(e) list(vector = rep(0, n)))
      if ("Hubs" %in% measures) results$Hubs <- hub_auth$vector
      if ("Authorities" %in% measures) results$Authorities <- auth_score$vector
    } else {
      if ("Hubs" %in% measures) results$Hubs <- rep(0, n)
      if ("Authorities" %in% measures) results$Authorities <- rep(0, n)
    }
  }
  results
}

is_valid_color <- function(z){ 
  if (is.null(z) || length(z) == 0) return(FALSE)
  tryCatch({ 
    grDevices::col2rgb(z); 
    TRUE 
  }, error = function(e) FALSE) 
}
sanitize_colors <- function(cols, bg="white"){
  if (length(cols) == 0) return(cols)
  cols <- as.character(cols)
  # Handle "background" keyword
  cols[trimws(cols) == "background"] <- bg
  # Validate colors
  bad <- !vapply(cols, is_valid_color, logical(1))
  if (any(bad)) cols[bad] <- "#808080"  # default to gray for invalid colors
  cols
}
make_pie_color_arg <- function(pie_input, n_nodes, bg="white"){
  if (is.null(pie_input) || n_nodes <= 0) return(NULL)
  # pie_input is now a hex color from colourpicker
  txt <- trimws(as.character(pie_input))
  # If it's a valid color, repeat it for all nodes
  if (is_valid_color(txt)) {
    cols <- rep(txt, n_nodes)
  } else {
    # Fallback to lightblue if invalid
    cols <- rep("#ADD8E6", n_nodes)
  }
  sanitize_colors(cols, bg = bg)
}
scale_to_range <- function(x, rng = c(4,10), reverse = FALSE){
  if (!length(x)) return(numeric())
  
  a <- suppressWarnings(min(x, na.rm=TRUE))
  b <- suppressWarnings(max(x, na.rm=TRUE))
  
  if (!is.finite(a) || !is.finite(b) || a==b) return(rep(mean(rng), length(x)))
  
  if (reverse) {
    # Reverse formula: (max - X) / (max - min)
    # Smallest value gets 1.0, largest gets 0.0
    normalized <- (b - x) / (b - a)
  } else {
    # Standard formula: (X - min) / (max - min)
    # Smallest value gets 0.0, largest gets 1.0
    normalized <- (x - a) / (b - a)
  }
  
  # Scale to target range
  normalized * (rng[2] - rng[1]) + rng[1]
}
apply_relimp_threshold <- function(net, thr){
  if (is.null(net) || !is.numeric(thr) || thr <= 0) return(net)
  if (!is.null(net$graph)){
    W <- as.matrix(net$graph); storage.mode(W) <- "double"; W[is.na(W)] <- 0
    W[abs(W) < thr] <- 0; net$graph <- W
  }
  net
}

make_relimp_matrix <- function(df, smart_preprocess = TRUE, dummy = FALSE,
                               drop_nonfinite = TRUE, drop_zerovar = TRUE){
  
  # CRITICAL: Sanitize column names for relimp (formulas can't handle special chars)
  # Save original names for later restoration
  original_names <- colnames(df)
  safe_names <- make.names(colnames(df), unique = TRUE)
  colnames(df) <- safe_names
  
  if (smart_preprocess) {
    X <- df
    for (i in seq_len(ncol(X))) {
      if (is_dichotomous(X[[i]])) {
        X[[i]] <- dichotomous_to_numeric(X[[i]])
      } else if (!is_numeric_like(X[[i]])) {
        if (dummy) {
          next
        } else {
          num_attempt <- suppressWarnings(as.numeric(as.character(X[[i]])))
          if (all(is.na(num_attempt))) {
            X[[i]] <- NULL
          } else {
            X[[i]] <- num_attempt
          }
        }
      }
    }
    if (dummy && any(sapply(X, function(x) !is_numeric_like(x)))) {
      X <- stats::model.matrix(~ . - 1, data = X)
      X <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
    }
    X <- X[, sapply(X, is_numeric_like), drop = FALSE]
  } else {
    X <- df[, sapply(df, is_numeric_like), drop = FALSE]
  }
  if (ncol(X) < 2) {
    n_numeric <- sum(sapply(df, is_numeric_like))
    n_dichot <- sum(sapply(df, is_dichotomous))
    stop(sprintf("After preprocessing: %d numeric columns, %d dichotomous columns. Need at least 2 total columns for relimp.", 
                 n_numeric, n_dichot))
  }
  if (drop_zerovar){
    keep <- vapply(X, function(col){
      col <- as.numeric(col)
      v <- stats::var(col, na.rm = TRUE)
      is.finite(v) && v > 1e-10
    }, logical(1))
    if (any(!keep)) {
      dropped_names <- names(X)[!keep]
      X <- X[, keep, drop = FALSE]
      if (length(dropped_names) > 0) {
        message(sprintf("Dropped %d zero-variance columns: %s", 
                        length(dropped_names), 
                        paste(head(dropped_names, 5), collapse=", ")))
      }
    }
  }
  if (drop_nonfinite){
    good <- apply(as.matrix(X), 1, function(r) all(is.finite(r)))
    n_dropped <- sum(!good)
    X <- X[good, , drop = FALSE]
    if (n_dropped > 0) {
      message(sprintf("Dropped %d rows with non-finite values", n_dropped))
    }
  }
  if (ncol(X) < 2 || nrow(X) < 3) {
    stop(sprintf("Relimp preprocessing resulted in %d columns and %d rows (need ≥2 columns and ≥3 rows).", 
                 ncol(X), nrow(X)))
  }
  
  # Store mapping of safe names to original names as attribute
  attr(X, "original_names") <- original_names[match(colnames(X), safe_names)]
  
  as.data.frame(X)
}

CENT_ALL <- c("Strength","InStrength","OutStrength","Degree","InDegree","OutDegree",
              "ExpectedInfluence","ExpectedInfluence2","Closeness","Betweenness",
              "Eigenvector","Hubs","Authorities")

# ---------- Layout computation helpers ----------
compute_ega_layout <- function(W, data = NULL, return_communities = FALSE, separation = 3) {
  # Compute EGA layout using community detection
  # separation: factor for within-community attraction (higher = more separated communities)
  tryCatch({
    # Convert to absolute values for community detection
    W_abs <- abs(W)
    diag(W_abs) <- 0
    
    communities <- NULL
    n_communities <- 1
    
    # Try EGAnet with raw data if available
    if (!is.null(data) && requireNamespace("EGAnet", quietly = TRUE)) {
      # Get numeric data only
      data_numeric <- data[, sapply(data, is.numeric), drop = FALSE]
      
      if (ncol(data_numeric) >= 2 && nrow(data_numeric) >= ncol(data_numeric)) {
        ega_result <- tryCatch({
          EGAnet::EGA(data_numeric, plot.EGA = FALSE)
        }, error = function(e) {
          message("EGAnet::EGA with raw data failed: ", e$message)
          NULL
        })
        
        if (!is.null(ega_result) && !is.null(ega_result$wc)) {
          communities <- ega_result$wc
          n_communities <- length(unique(communities))
          
          # Create community-aware layout
          if (n_communities > 1) {
            # Use igraph's Fruchterman-Reingold with modified weights
            g <- igraph::graph_from_adjacency_matrix(W_abs, mode = "undirected", weighted = TRUE)
            
            # Modify weights based on community membership
            W_community <- W_abs
            within_factor <- separation  # Stronger attraction within communities
            between_factor <- 1 / separation  # Weaker attraction between communities
            
            for (i in 1:nrow(W_community)) {
              for (j in 1:ncol(W_community)) {
                if (communities[i] == communities[j]) {
                  # Increase weight for within-community edges (attract)
                  W_community[i, j] <- W_community[i, j] * within_factor
                } else {
                  # Decrease weight for between-community edges (repel)
                  W_community[i, j] <- W_community[i, j] * between_factor
                }
              }
            }
            
            g_comm <- igraph::graph_from_adjacency_matrix(W_community, mode = "undirected", weighted = TRUE)
            layout <- igraph::layout_with_fr(g_comm, weights = igraph::E(g_comm)$weight)
            
          } else {
            # Regular layout if only 1 community
            g <- igraph::graph_from_adjacency_matrix(W_abs, mode = "undirected", weighted = TRUE)
            layout <- igraph::layout_with_fr(g)
          }
          
          message(sprintf("EGA detected %d communities from raw data (separation = %.1f)", n_communities, separation))
          
          if (return_communities) {
            return(list(layout = layout, communities = communities, n_communities = n_communities))
          }
          return(layout)
        }
      }
    }
    
    # Fallback to igraph community detection with multiple algorithms
    if (requireNamespace("igraph", quietly = TRUE)) {
      g <- igraph::graph_from_adjacency_matrix(W_abs, mode = "undirected", weighted = TRUE)
      
      # Try multiple community detection algorithms
      comm_obj <- NULL
      algorithm_used <- "none"
      
      # 1. Try Louvain (usually best for modularity)
      comm_obj <- tryCatch({
        algorithm_used <- "Louvain"
        igraph::cluster_louvain(g)
      }, error = function(e) NULL)
      
      # 2. Try Walktrap if Louvain fails
      if (is.null(comm_obj)) {
        comm_obj <- tryCatch({
          algorithm_used <- "Walktrap"
          igraph::cluster_walktrap(g)
        }, error = function(e) NULL)
      }
      
      # 3. Try Fast-greedy if both fail
      if (is.null(comm_obj)) {
        comm_obj <- tryCatch({
          algorithm_used <- "Fast-greedy"
          igraph::cluster_fast_greedy(g)
        }, error = function(e) NULL)
      }
      
      if (!is.null(comm_obj)) {
        communities <- igraph::membership(comm_obj)
        n_communities <- length(unique(communities))
        
        message(sprintf("%s algorithm detected %d communities from network (separation = %.1f)", 
                        algorithm_used, n_communities, separation))
        
        # Create community-aware layout
        if (n_communities > 1) {
          # Modify weights to emphasize community structure
          W_community <- W_abs
          within_factor <- separation
          between_factor <- 1 / separation
          
          for (i in 1:nrow(W_community)) {
            for (j in 1:ncol(W_community)) {
              if (communities[i] == communities[j]) {
                # Increase weight for within-community edges
                W_community[i, j] <- W_community[i, j] * within_factor
              } else {
                # Decrease weight for between-community edges
                W_community[i, j] <- W_community[i, j] * between_factor
              }
            }
          }
          
          g_comm <- igraph::graph_from_adjacency_matrix(W_community, mode = "undirected", weighted = TRUE)
          layout <- igraph::layout_with_fr(g_comm, weights = igraph::E(g_comm)$weight)
          
        } else {
          # Just use regular FR layout if only 1 community
          layout <- igraph::layout_with_fr(g)
        }
        
        if (return_communities) {
          return(list(layout = layout, communities = communities, n_communities = n_communities))
        }
        return(layout)
      }
    }
    
    # Final fallback to spring
    message("Community detection failed, using spring layout")
    if (return_communities) {
      return(list(layout = "spring", communities = NULL, n_communities = 1))
    }
    return("spring")
  }, error = function(e) {
    warning("EGA layout failed, using spring layout: ", e$message)
    if (return_communities) {
      return(list(layout = "spring", communities = NULL, n_communities = 1))
    }
    return("spring")
  })
}

detect_bridge_communities <- function(W, algo = "walktrap") {
  g <- tryCatch(
    igraph::graph_from_adjacency_matrix(abs(W), mode = "undirected", weighted = TRUE),
    error = function(e) NULL
  )
  if (is.null(g)) return(NULL)

  if (algo == "ega") {
    df_num <- as.data.frame(W)
    ega <- tryCatch(EGAnet::EGA(df_num, plot.EGA = FALSE), error = function(e) NULL)
    if (!is.null(ega) && !is.null(ega$wc)) {
      return(as.character(ega$wc))
    }
    # Fallback to walktrap if EGA fails
    algo <- "walktrap"
  }

  comm_obj <- switch(algo,
    "walktrap" = tryCatch(igraph::cluster_walktrap(g), error = function(e) NULL),
    "louvain"  = tryCatch(igraph::cluster_louvain(g),  error = function(e) NULL),
    NULL
  )
  # Fallback to walktrap if algo failed
  if (is.null(comm_obj))
    comm_obj <- tryCatch(igraph::cluster_walktrap(g), error = function(e) NULL)
  if (is.null(comm_obj)) return(NULL)

  as.character(igraph::membership(comm_obj))
}

compute_pca_layout <- function(W, data = NULL) {
  # Compute PCA-based layout
  tryCatch({
    # If we have original data, use that for PCA
    if (!is.null(data)) {
      data_numeric <- data[, sapply(data, is.numeric), drop = FALSE]
      if (ncol(data_numeric) >= 2 && nrow(data_numeric) > ncol(data_numeric)) {
        pca <- prcomp(data_numeric, scale. = TRUE, center = TRUE)
        # Use first two components
        layout <- pca$rotation[, 1:2, drop = FALSE]
        colnames(layout) <- c("x", "y")
        return(layout)
      }
    }
    
    # Otherwise use network structure
    # Perform eigen decomposition of the network
    W_sym <- (W + t(W)) / 2  # Make symmetric
    diag(W_sym) <- 0
    
    # Get eigenvalues/vectors
    eigen_result <- eigen(W_sym)
    
    # Use first two eigenvectors as coordinates
    layout <- eigen_result$vectors[, 1:2, drop = FALSE]
    
    # Scale to reasonable range
    layout <- scale(layout) * 2
    colnames(layout) <- c("x", "y")
    
    return(layout)
  }, error = function(e) {
    warning("PCA layout failed, using spring layout: ", e$message)
    return("spring")
  })
}

# ============================================================
# NSON: Nested Specificity-Oriented Network – helper functions
# ============================================================

compute_predictability <- function(data, W, type = c("continuous", "binary")) {
  type <- match.arg(type)
  data <- as.data.frame(data)
  var_names <- colnames(W)
  p <- length(var_names)
  pred <- numeric(p)
  for (j in seq_len(p)) {
    y_nm <- var_names[j]
    neighbors <- var_names[-j][abs(W[j, -j]) > 1e-10]
    if (length(neighbors) == 0) { pred[j] <- 0; next }
    df_fit <- data[, c(y_nm, neighbors), drop = FALSE]
    df_fit <- df_fit[complete.cases(df_fit), , drop = FALSE]
    if (nrow(df_fit) < length(neighbors) + 2) { pred[j] <- NA_real_; next }
    is_bin <- type == "binary" || all(na.omit(df_fit[[y_nm]]) %in% c(0, 1))
    fml <- as.formula(paste("`", y_nm, "`", " ~ ", paste(paste0("`", neighbors, "`"), collapse = "+"), sep=""))
    if (is_bin) {
      fit <- tryCatch(glm(fml, data = df_fit, family = binomial()), error = function(e) NULL)
      if (is.null(fit)) { pred[j] <- NA_real_; next }
      pred[j] <- round(1 - fit$deviance / fit$null.deviance, 4)
    } else {
      fit <- tryCatch(lm(fml, data = df_fit), error = function(e) NULL)
      if (is.null(fit)) { pred[j] <- NA_real_; next }
      pred[j] <- round(summary(fit)$r.squared, 4)
    }
  }
  data.frame(variable = var_names, predictability = pred, stringsAsFactors = FALSE)
}

compute_gnss <- function(data, type = c("binary", "continuous"),
                         z_threshold = 0,
                         missing = c("pairwise", "listwise")) {
  type    <- match.arg(type)
  missing <- match.arg(missing)
  data    <- as.data.frame(data)

  if (type == "continuous") {
    data_z  <- as.data.frame(scale(data))
    act_mat <- as.data.frame((data_z > z_threshold) * 1.0)
    colnames(act_mat) <- colnames(data)
  } else {
    act_mat <- data
  }

  p         <- ncol(act_mat)
  var_names <- colnames(act_mat)
  gnss_vals <- numeric(p)
  act_rates <- numeric(p)

  for (j in seq_len(p)) {
    y_col <- act_mat[[j]]
    if (missing == "listwise") {
      complete_rows <- complete.cases(act_mat)
      act_rates[j]  <- mean(y_col[complete_rows] == 1, na.rm = TRUE)
      mask          <- complete_rows & !is.na(y_col) & (y_col == 1)
    } else {
      act_rates[j]  <- mean(y_col == 1, na.rm = TRUE)
      mask          <- !is.na(y_col) & (y_col == 1)
    }
    if (sum(mask) == 0) { gnss_vals[j] <- NA; next }
    subset_rows      <- act_mat[mask, -j, drop = FALSE]
    cond_activations <- colMeans(subset_rows == 1, na.rm = TRUE)
    gnss_vals[j]     <- mean(cond_activations, na.rm = TRUE)

    # Console audit: show all P(X=1 | Y=1) that were averaged into GNSS(Y)
    cond_lines <- paste(
      sprintf("    P(%s=1 | %s=1) = %.4f", names(cond_activations), var_names[j], cond_activations),
      collapse = "\n")
    message(sprintf(
      "\nGNSS audit — %s  [n_active=%d, activation_rate=%.4f]\n%s\n  --> GNSS(%s) = %.4f",
      var_names[j], sum(mask), act_rates[j], cond_lines, var_names[j], gnss_vals[j]))
  }

  data.frame(
    variable        = var_names,
    activation_rate = round(act_rates, 4),
    gnss            = round(gnss_vals, 4),
    GNSS_rank       = rank(-gnss_vals, na.last = "keep", ties.method = "average"),
    stringsAsFactors = FALSE
  )
}

compute_kway_gnss <- function(data, type = c("binary","continuous"),
                               z_threshold = 0, max_k = 4L,
                               missing = c("pairwise","listwise"),
                               max_combos = 1000L) {
  type    <- match.arg(type)
  missing <- match.arg(missing)
  data    <- as.data.frame(data)
  data    <- data[, vapply(data, function(x) is.numeric(x) || is.logical(x), logical(1)),
                  drop = FALSE]
  if (ncol(data) < 2) return(NULL)

  if (type == "continuous") {
    dz      <- as.data.frame(scale(data))
    act_mat <- as.matrix((dz > z_threshold) * 1.0)
    colnames(act_mat) <- colnames(data)
  } else {
    act_mat <- as.matrix(data); storage.mode(act_mat) <- "double"
  }

  p         <- ncol(act_mat)
  var_names <- colnames(act_mat)
  max_k     <- min(as.integer(max_k), p - 1L)
  if (max_k < 2L) return(NULL)

  k_seq <- seq(2L, max_k)
  rows  <- vector("list", p * length(k_seq))
  ri    <- 0L

  for (j in seq_len(p)) {
    y   <- act_mat[, j]
    if (missing == "listwise") {
      ok   <- complete.cases(act_mat)
      mask <- ok & !is.na(y) & (y == 1)
    } else {
      mask <- !is.na(y) & (y == 1)
    }
    n_act    <- sum(mask)
    peers    <- seq_len(p)[-j]
    sub_full <- act_mat[mask, peers, drop = FALSE]

    for (k in k_seq) {
      ri   <- ri + 1L
      km1  <- k - 1L

      if (n_act == 0L || length(peers) < km1) {
        rows[[ri]] <- data.frame(variable=var_names[j], k=k,
                                 and_gnss=NA_real_, or_gnss=NA_real_,
                                 n_active=n_act, n_combos=0L,
                                 stringsAsFactors=FALSE)
        next
      }

      all_cb  <- combn(length(peers), km1)          # km1 × C(p-1,km1)
      n_total <- ncol(all_cb)
      if (n_total > max_combos) {
        samp   <- sample(n_total, max_combos, replace=FALSE)
        all_cb <- all_cb[, samp, drop=FALSE]
        n_used <- max_combos
      } else {
        n_used <- n_total
      }

      and_p <- vapply(seq_len(ncol(all_cb)), function(ci) {
        sub <- sub_full[, all_cb[, ci], drop=FALSE]
        mean(rowSums(sub == 1) == km1, na.rm=TRUE)
      }, numeric(1))

      or_p <- vapply(seq_len(ncol(all_cb)), function(ci) {
        sub <- sub_full[, all_cb[, ci], drop=FALSE]
        mean(rowSums(sub == 1) > 0L, na.rm=TRUE)
      }, numeric(1))

      rows[[ri]] <- data.frame(
        variable = var_names[j], k = k,
        and_gnss = round(mean(and_p, na.rm=TRUE), 4),
        or_gnss  = round(mean(or_p,  na.rm=TRUE), 4),
        n_active = n_act, n_combos = n_used,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows[seq_len(ri)])
}

compute_cond_prob_matrix <- function(data, type = c("binary","continuous"),
                                      z_threshold = 0,
                                      missing = c("pairwise","listwise")) {
  type    <- match.arg(type)
  missing <- match.arg(missing)
  data    <- as.data.frame(data)

  if (type == "continuous") {
    dz      <- as.data.frame(scale(data))
    act_mat <- as.matrix((dz > z_threshold) * 1.0)
    colnames(act_mat) <- colnames(data)
  } else {
    act_mat <- as.matrix(data); storage.mode(act_mat) <- "double"
  }

  p         <- ncol(act_mat)
  var_names <- colnames(act_mat)
  # cp[i,j] = P(X_j=1 | X_i=1): row = conditioning variable, col = target
  cp_mat    <- matrix(NA_real_, p, p, dimnames = list(var_names, var_names))

  for (i in seq_len(p)) {
    y_col <- act_mat[, i]
    mask  <- if (missing == "listwise") {
      complete.cases(act_mat) & !is.na(y_col) & (y_col == 1)
    } else {
      !is.na(y_col) & (y_col == 1)
    }
    if (sum(mask) == 0) next
    sub         <- act_mat[mask, , drop = FALSE]
    cp_mat[i, ] <- colMeans(sub == 1, na.rm = TRUE)
    cp_mat[i, i] <- NA_real_   # self-conditional is trivially 1; exclude
  }
  cp_mat
}

orient_network_by_gnss <- function(weight_matrix, gnss_table, tolerance = 0) {
  W         <- as.matrix(weight_matrix)
  p         <- ncol(W)
  var_names <- colnames(W)
  gnss_vals <- setNames(gnss_table$gnss, gnss_table$variable)

  rows <- vector("list", 0L)
  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      w <- W[i, j]
      if (is.na(w) || abs(w) < 1e-10) next
      vi <- var_names[i]; vj <- var_names[j]
      gi <- gnss_vals[vi]; gj <- gnss_vals[vj]
      if (is.na(gi) || is.na(gj)) {
        direction <- "undirected (NA GNSS)"
      } else if (abs(gi - gj) <= tolerance) {
        direction <- "undirected"
      } else if (gi < gj) {
        direction <- paste0(vi, " -> ", vj)
      } else {
        direction <- paste0(vj, " -> ", vi)
      }
      rows[[length(rows) + 1]] <- data.frame(
        from          = vi,
        to            = vj,
        weight        = round(w, 4),
        sign          = ifelse(w > 0, "positive", "negative"),
        gnss_from     = if (!is.na(gi)) round(gi, 4) else NA_real_,
        gnss_to       = if (!is.na(gj)) round(gj, 4) else NA_real_,
        direction     = direction,
        abs_gnss_diff = if (!is.na(gi) && !is.na(gj)) round(abs(gi - gj), 4) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) return(data.frame(
    from = character(), to = character(), weight = numeric(), sign = character(),
    gnss_from = numeric(), gnss_to = numeric(), direction = character(),
    abs_gnss_diff = numeric(), stringsAsFactors = FALSE
  ))
  do.call(rbind, rows)
}

# ============================================================
# ENHANCED UI with Modern Styling
# ============================================================
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      /* Modern Color Palette */
      :root {
        --primary-color: #2C3E50;
        --secondary-color: #3498DB;
        --accent-color: #E74C3C;
        --success-color: #27AE60;
        --warning-color: #F39C12;
        --light-bg: #ECF0F1;
        --card-bg: #FFFFFF;
        --border-color: #BDC3C7;
        --text-dark: #2C3E50;
        --text-light: #7F8C8D;
      }
      
      /* Global Styling */
      body {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        background: white;
        color: var(--text-dark);
        padding: 20px;
      }
      
      .container-fluid {
        max-width: 1400px;
        margin: 0 auto;
      }
      
      /* Main Title */
      .title-panel {
        background: linear-gradient(135deg, #2C3E50 0%, #3498DB 100%);
        color: white;
        padding: 30px;
        border-radius: 15px;
        box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        margin-bottom: 30px;
        text-align: center;
      }
      
      .title-panel h2 {
        margin: 0;
        font-weight: 300;
        font-size: 2.5em;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
      }
      
      .title-panel .subtitle {
        margin-top: 10px;
        font-size: 1.1em;
        opacity: 0.9;
      }
      
      /* Tab Navigation */
      .nav-tabs {
        border-bottom: 3px solid var(--secondary-color);
        background: white;
        padding: 10px 10px 0 10px;
        border-radius: 10px 10px 0 0;
      }
      
      .nav-tabs > li > a {
        color: var(--text-dark);
        font-weight: 500;
        border: none;
        padding: 12px 25px;
        margin-right: 5px;
        border-radius: 8px 8px 0 0;
        transition: all 0.3s ease;
      }
      
      .nav-tabs > li > a:hover {
        background: var(--light-bg);
        border: none;
      }
      
      .nav-tabs > li.active > a {
        background: var(--secondary-color) !important;
        color: white !important;
        border: none;
        box-shadow: 0 -3px 10px rgba(52, 152, 219, 0.3);
      }
      
      .tab-content {
        background: white;
        padding: 30px;
        border-radius: 0 0 10px 10px;
        box-shadow: 0 5px 20px rgba(0,0,0,0.1);
      }
      
      /* Card/Panel Styling */
      .well, .panel {
        background: var(--card-bg);
        border: 1px solid var(--border-color);
        border-radius: 10px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.08);
        padding: 20px;
        margin-bottom: 20px;
      }
      
      .info-box {
        background: linear-gradient(135deg, #2C3E50 0%, #3498DB 100%);
        color: white;
        padding: 20px;
        border-radius: 10px;
        margin-bottom: 20px;
        box-shadow: 0 4px 15px rgba(52, 152, 219, 0.4);
      }
      
      .info-box h4 {
        margin-top: 0;
        font-weight: 600;
        border-bottom: 2px solid rgba(255,255,255,0.3);
        padding-bottom: 10px;
      }
      
      /* Split Variable Box */
      .split-box {
        border: 2px dashed var(--warning-color);
        padding: 20px;
        margin: 15px 0;
        border-radius: 10px;
        background: linear-gradient(135deg, #FFF3E0 0%, #FFE0B2 100%);
        box-shadow: 0 2px 10px rgba(243, 156, 18, 0.2);
      }
      
      .split-box h4 {
        color: var(--warning-color);
        margin-top: 0;
        font-weight: 600;
      }
      
      /* Variable Selection Boxes */
      .var-selection-box {
        background: var(--light-bg);
        border-radius: 10px;
        padding: 15px;
        min-height: 400px;
      }
      
      .var-selection-box h4 {
        background: var(--primary-color);
        color: white;
        padding: 10px 15px;
        margin: -15px -15px 15px -15px;
        border-radius: 10px 10px 0 0;
        font-weight: 500;
      }
      
      .var-selection-box select {
        width: 100%;
        border: 2px solid var(--border-color);
        border-radius: 8px;
        padding: 10px;
        font-size: 14px;
        background: white;
      }
      
      .var-selection-box select:focus {
        border-color: var(--secondary-color);
        outline: none;
        box-shadow: 0 0 0 3px rgba(52, 152, 219, 0.1);
      }
      
      /* Buttons */
      .btn {
        border-radius: 8px;
        padding: 10px 20px;
        font-weight: 500;
        transition: all 0.3s ease;
        border: none;
        box-shadow: 0 2px 5px rgba(0,0,0,0.1);
      }
      
      .btn:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 10px rgba(0,0,0,0.2);
      }
      
      .btn-primary {
        background: linear-gradient(135deg, var(--secondary-color) 0%, #2980B9 100%);
        color: white;
      }
      
      .btn-primary:hover {
        background: linear-gradient(135deg, #2980B9 0%, #2471A3 100%);
      }
      
      .btn-warning {
        background: linear-gradient(135deg, var(--warning-color) 0%, #E67E22 100%);
        color: white;
      }
      
      .btn-warning:hover {
        background: linear-gradient(135deg, #E67E22 0%, #D35400 100%);
      }
      
      .btn-success {
        background: linear-gradient(135deg, var(--success-color) 0%, #229954 100%);
        color: white;
      }
      
      .btn-default {
        background: white;
        color: var(--text-dark);
        border: 2px solid var(--border-color);
      }
      
      .btn-default:hover {
        background: var(--light-bg);
        border-color: var(--primary-color);
      }
      
      .btn-sm {
        padding: 6px 12px;
        font-size: 13px;
      }
      
      .btn-block {
        width: 100%;
        margin-top: 10px;
        white-space: normal !important;
        word-break: break-word;
      }
      
      /* Form Controls */
      .form-control, .selectize-input {
        border: 2px solid var(--border-color);
        border-radius: 8px;
        padding: 10px;
        transition: all 0.3s ease;
      }
      
      .form-control:focus, .selectize-input.focus {
        border-color: var(--secondary-color);
        box-shadow: 0 0 0 3px rgba(52, 152, 219, 0.1);
        outline: none;
      }
      
      /* File Input */
      .btn-file {
        background: linear-gradient(135deg, var(--secondary-color) 0%, #2980B9 100%);
        color: white;
      }
      
      /* Checkboxes and Radio Buttons */
      .checkbox label, .radio label {
        font-weight: 500;
        color: var(--text-dark);
      }
      
      /* Output Boxes */
      pre, .shiny-text-output {
        background: #F8F9FA;
        border: 1px solid var(--border-color);
        border-radius: 8px;
        padding: 15px;
        font-family: 'Consolas', 'Monaco', monospace;
        font-size: 13px;
        color: var(--text-dark);
        max-height: 400px;
        overflow-y: auto;
      }
      
      /* DataTable Styling */
      .dataTables_wrapper {
        padding: 15px;
        background: white;
        border-radius: 10px;
        overflow-x: auto;
        max-width: 100%;
      }
      
      table.dataTable {
        border-collapse: collapse !important;
        width: 100% !important;
        max-width: 100%;
      }
      
      table.dataTable thead th {
        background: var(--primary-color);
        color: white;
        font-weight: 500;
        padding: 12px;
        white-space: nowrap;
      }
      
      table.dataTable tbody td {
        padding: 10px;
        border-bottom: 1px solid var(--border-color);
      }
      
      table.dataTable tbody tr:hover {
        background: var(--light-bg);
      }
      
      /* Table container for regular tables */
      .shiny-html-output table {
        max-width: 100%;
        overflow-x: auto;
        display: block;
      }
      
      .table-responsive {
        overflow-x: auto;
        max-width: 100%;
        border: 1px solid var(--border-color);
        border-radius: 8px;
      }
      
      /* Sidebar Styling */
      .well {
        background: linear-gradient(to bottom, #ffffff 0%, #f8f9fa 100%);
      }
      
      /* Status Boxes */
      .status-good {
        background: #D5F4E6;
        border-left: 4px solid var(--success-color);
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
      }
      
      .status-warning {
        background: #FFF3E0;
        border-left: 4px solid var(--warning-color);
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
      }
      
      .status-error {
        background: #FFEBEE;
        border-left: 4px solid var(--accent-color);
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
      }
      
      /* Help Text */
      .help-block {
        color: var(--text-light);
        font-size: 12px;
        font-style: italic;
        margin-top: 5px;
      }
      
      /* Compact Form Groups */
      .panel-compact .form-group {
        margin-bottom: 8px;
      }
      
      /* Section Headers */
      hr {
        border: none;
        height: 2px;
        background: linear-gradient(to right, transparent, var(--border-color), transparent);
        margin: 25px 0;
      }
      
      h4 {
        color: var(--primary-color);
        font-weight: 600;
        margin-top: 20px;
        margin-bottom: 15px;
      }
      
      /* Action Button Group */
      .btn-group-custom {
        display: flex;
        gap: 10px;
        margin: 10px 0;
      }
      
      /* Plot Output */
      .shiny-plot-output {
        border: 1px solid var(--border-color);
        border-radius: 10px;
        padding: 10px;
        background: white;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
      }
      
      /* Responsive adjustments */
      @media (max-width: 768px) {
        .title-panel h2 {
          font-size: 1.8em;
        }
        
        .tab-content {
          padding: 15px;
        }
      }
      
      /* Loading indicator */
      .shiny-busy-indicator {
        position: fixed;
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);
        background: white;
        padding: 30px;
        border-radius: 15px;
        box-shadow: 0 10px 40px rgba(0,0,0,0.3);
        z-index: 10000;
      }
      
      /* Scrollbar styling */
      ::-webkit-scrollbar {
        width: 10px;
        height: 10px;
      }
      
      ::-webkit-scrollbar-track {
        background: var(--light-bg);
        border-radius: 5px;
      }
      
      ::-webkit-scrollbar-thumb {
        background: var(--secondary-color);
        border-radius: 5px;
      }
      
      ::-webkit-scrollbar-thumb:hover {
        background: #2980B9;
      }
      
      /* Color picker styling */
      .colourpicker-input {
        display: flex;
        align-items: center;
        gap: 10px;
      }
      
      .colourpicker-input input[type='text'] {
        font-family: 'Consolas', 'Monaco', monospace;
        font-size: 13px;
      }
    "))
  ),
  
  # Title
  div(class = "title-panel",
      h2("Network Analysis Suite"),
      div(class = "subtitle", "Advanced bootnet estimation with split-group comparison")
  ),
  
  tabsetPanel(id = "tabs",
              
              # ---------------- Data Input ----------------
              tabPanel("Data Input",
                       icon = icon("upload"),
                       br(),
                       fluidRow(
                         column(4,
                                div(class = "info-box",
                                    h4("Upload Your Data"),
                                    fileInput("file", NULL, 
                                              accept = c(".csv",".tsv",".txt",".xlsx",".xls"),
                                              buttonLabel = "Browse...",
                                              placeholder = "No file selected"),
                                    uiOutput("sheet_picker"),
                                    actionButton("clear_data", "Clear / Reset", 
                                                 class = "btn-warning btn-block",
                                                 icon = icon("refresh"))
                                )
                         ),
                         column(4,
                                div(class = "var-selection-box",
                                    h4("Import Options"),
                                    checkboxInput("show_advanced", "Show advanced options", FALSE),
                                    conditionalPanel("input.show_advanced",
                                                     radioButtons("txt_delim", "Delimiter", inline = TRUE,
                                                                  choices = c("Auto"="auto","Comma"=",","Semicolon"=";","Tab"="\t","Space"=" "), 
                                                                  selected="auto"),
                                                     textInput("na_tokens", "NA strings", value = "NA,NaN,.,"),
                                                     checkboxInput("has_header", "First row = headers", TRUE),
                                                     radioButtons("decimal_mark", "Decimal mark", inline = TRUE,
                                                                  choices = c("Dot (.)" = ".", "Comma (,)" = ","), selected = "."),
                                                     textInput("encoding", "Encoding", value = "UTF-8")
                                    ),
                                    checkboxInput("show_preview", "Show preview", TRUE)
                                )
                         ),
                         column(4,
                                div(class = "info-box",
                                    h4("Import Status"),
                                    verbatimTextOutput("import_status", placeholder = TRUE)
                                )
                         )
                       ),
                       hr(),
                       conditionalPanel("input.show_preview",
                                        h4("Data Preview (first 1000 rows)"), 
                                        DT::DTOutput("previewDT"),
                                        br(), 
                                        h4("Summary Statistics"), 
                                        verbatimTextOutput("preview_summary")
                       )
              ),
              
              # ---------------- Variables ----------------
              tabPanel("Variables",
                       icon = icon("sliders"),
                       br(),
                       div(class = "info-box",
                           h4("ℹ️ Variable Type Guide"),
                           tags$ul(
                             tags$li(tags$b("Continuous (Gaussian):"), " Real-valued variables (e.g., age, weight, scores). Used by all estimators."),
                             tags$li(tags$b("Categorical:"), " Unordered categories (e.g., gender, diagnosis). For MGM only."),
                             tags$li(tags$b("Binary/Count (Poisson):"), " Non-negative integers (e.g., number of events, frequency). For MGM only.")
                           ),
                           helpText("⚠️ EBICglasso and pcor require ALL variables to be Continuous. MGM can mix all three types.")
                       ),
                       br(),
                       fluidRow(
                         column(2, class = "panel-compact",
                                div(class = "var-selection-box",
                                    h4("Available"),
                                    selectInput("available_vars", NULL, choices = character(0),
                                                multiple = TRUE, size = 14, width = "100%", selectize = FALSE),
                                    actionButton("auto_detect", "🔍 Auto-detect",
                                                 class = "btn-success btn-block btn-sm"),
                                    helpText("💡 Cmd/Ctrl to multi-select")
                                )
                         ),
                         column(2, class = "panel-compact",
                                div(class = "var-selection-box",
                                    h4("Continuous (all)"),
                                    selectInput("continuous_vars", NULL, choices = character(0),
                                                multiple = TRUE, size = 10, width = "100%", selectize = FALSE),
                                    div(class = "btn-group-custom",
                                        actionButton("add_continuous", "Add", class = "btn btn-primary btn-sm"),
                                        actionButton("remove_continuous", "Remove", class = "btn btn-default btn-sm")
                                    ),
                                    helpText("Gaussian (g)")
                                )
                         ),
                         column(2, class = "panel-compact",
                                div(class = "var-selection-box",
                                    h4("Categorical (MGM)"),
                                    selectInput("discrete_vars", NULL, choices = character(0),
                                                multiple = TRUE, size = 10, width = "100%", selectize = FALSE),
                                    div(class = "btn-group-custom",
                                        actionButton("add_discrete", "Add", class = "btn btn-primary btn-sm"),
                                        actionButton("remove_discrete", "Remove", class = "btn btn-default btn-sm")
                                    ),
                                    helpText("Categorical (c)")
                                )
                         ),
                         column(2, class = "panel-compact",
                                div(class = "var-selection-box",
                                    h4("Count (MGM)"),
                                    selectInput("count_vars", NULL, choices = character(0),
                                                multiple = TRUE, size = 10, width = "100%", selectize = FALSE),
                                    div(class = "btn-group-custom",
                                        actionButton("add_count", "Add", class = "btn btn-primary btn-sm"),
                                        actionButton("remove_count", "Remove", class = "btn btn-default btn-sm")
                                    ),
                                    helpText("Poisson (p)")
                                )
                         ),
                         column(4, class = "panel-compact",
                                div(class = "split-box",
                                    h4("✂️ Split Variable"),
                                    selectInput("split_var", "Select split variable:",
                                                choices = c("None" = ""), width = "100%"),
                                    radioButtons("split_sample_mode", "Sample size:",
                                                choices = c("Use all data" = "all",
                                                          "Equal sample sizes" = "equal"),
                                                selected = "all",
                                                inline = TRUE),
                                    conditionalPanel(
                                      "input.split_sample_mode == 'equal'",
                                      checkboxInput("split_random_seed", "Use fixed seed (reproducible)", value = TRUE),
                                      conditionalPanel(
                                        "input.split_random_seed == true",
                                        numericInput("split_seed_value", "Seed value:", value = 123, min = 1, max = 10000, step = 1)
                                      ),
                                      helpText(HTML("<strong>ℹ️ Equal sampling:</strong><br/>
                                                    Randomly samples from larger group<br/>
                                                    to match smallest group size.<br/>
                                                    <strong>Fixed seed:</strong> Same selection each run (reproducible)<br/>
                                                    <strong>No seed:</strong> Different selection each run (truly random)"))
                                    ),
                                    div(class = "btn-group-custom",
                                        actionButton("set_split", "Set Split", class = "btn btn-warning btn-sm"),
                                        actionButton("clear_split", "Clear", class = "btn btn-default btn-sm")
                                    ),
                                    helpText("🔀 Compare networks across groups")
                                )
                         )
                       ),
                       hr(),
                       fluidRow(
                         column(12,
                                div(class = "well",
                                    h4("🔧 Missing Data Handling"),
                                    radioButtons("missing_policy", "Missing data policy", inline = TRUE,
                                                choices = c("Listwise deletion" = "listwise",
                                                           "Pairwise deletion" = "pairwise",
                                                           "Keep (as-is)" = "keep"),
                                                selected = "listwise"),
                                    helpText(HTML("<strong>Listwise:</strong> Remove rows with any missing values<br/>
                                                  <strong>Pairwise:</strong> Use all available pairs for correlations (works with EBICglasso, pcor, cor)<br/>
                                                  <strong>Keep:</strong> Pass missing values to estimator (may cause errors)"))
                                )
                         )
                       ),
                       hr(),
                       fluidRow(
                         column(9,
                                div(class = "status-good",
                                    icon("check-circle"),
                                    verbatimTextOutput("vars_status")
                                )
                         ),
                         column(3,
                                div(class = "status-warning",
                                    icon("info-circle"),
                                    verbatimTextOutput("split_status")
                                )
                         )
                       )
              ),
              
              # -------------- Estimation & Visualization --------------
              tabPanel("Estimation & Visualization",
                       icon = icon("network-wired"),
                       br(),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      style = "position: sticky; top: 0; max-height: 100vh; overflow-y: auto;",
                                      div(class = "well",
                                          h4("🔧 Estimation"),
                                          shinyWidgets::pickerInput(
                                            "estimator", "Estimator",
                                            choices = c("Correlation Network"="cor",
                                                        "Partial Correlation"="pcor",
                                                        "EBICglasso (GGM)"="EBICglasso",
                                                        "Huge (High-dimensional)"="huge",
                                                        "IsingFit (Binary data)"="IsingFit",
                                                        "Mixed Graphical Model"="mgm",
                                                        "Relative Importance"="relimp"),
                                            selected = "EBICglasso"
                                          ),
                                          
                                          conditionalPanel("input.estimator == 'cor' || input.estimator == 'pcor'",
                                                           selectInput("corMethod", "Correlation method",
                                                                       choices = c("Auto (cor_auto)" = "cor_auto",
                                                                                 "Pearson (cor)" = "cor",
                                                                                 "Covariance (cov)" = "cov",
                                                                                 "Nonparanormal (npn)" = "npn",
                                                                                 "Spearman" = "spearman"),
                                                                       selected = "cor_auto"),
                                                           selectInput("threshold_method", "Threshold",
                                                                      choices = c("None (show all)" = "none",
                                                                                "Manual cutoff" = "manual",
                                                                                "Significant (p < 0.05)" = "sig",
                                                                                "Bonferroni" = "bonferroni",
                                                                                "Holm" = "holm",
                                                                                "Hochberg" = "hochberg",
                                                                                "Hommel" = "hommel",
                                                                                "BH (False Discovery Rate)" = "BH"),
                                                                      selected = "none"),
                                                           conditionalPanel("input.threshold_method == 'manual'",
                                                                          numericInput("threshold_value", "Cutoff value", value = 0.1, min = 0, max = 1, step = 0.05)
                                                           )
                                          ),

                                          conditionalPanel("input.estimator == 'EBICglasso'",
                                                           selectInput("corMethod_ebic", "Correlation method",
                                                                       choices = c("Auto (cor_auto)" = "cor_auto",
                                                                                 "Pearson (cor)" = "cor",
                                                                                 "Covariance (cov)" = "cov",
                                                                                 "Nonparanormal (npn)" = "npn",
                                                                                 "Spearman" = "spearman"),
                                                                       selected = "cor"),
                                                           numericInput("tuning", "EBIC gamma", value = 0.5, step = 0.05),
                                                           checkboxInput("refit", "Refit without LASSO", FALSE),
                                                           selectInput("ebic_samplesize", "Sample size for EBIC (with pairwise deletion)",
                                                                      choices = c("Maximum (total N, like JASP)" = "maximum",
                                                                                "Minimum (complete cases only)" = "minimum",
                                                                                "Pairwise average" = "pairwise_average"),
                                                                      selected = "maximum"),
                                                           helpText("Affects EBIC penalty when using pairwise deletion. 'Maximum' matches JASP's default behavior.")
                                          ),

                                          conditionalPanel("input.estimator == 'huge'",
                                                           selectInput("huge_method", "Method",
                                                                      choices = c("Graphical Lasso (glasso)" = "glasso",
                                                                                "Meinshausen-Bühlmann (mb)" = "mb",
                                                                                "Column-wise regression (ct)" = "ct"),
                                                                      selected = "glasso"),
                                                           helpText("glasso = graphical lasso; mb = neighborhood selection; ct = column-wise regression"),
                                                           selectInput("huge_criterion", "Selection criterion",
                                                                      choices = c("EBIC (Extended BIC)" = "ebic",
                                                                                "RIC (Rotation Information Criterion)" = "ric",
                                                                                "STARS (Stability Selection)" = "stars"),
                                                                      selected = "stars"),
                                                           conditionalPanel("input.huge_criterion == 'ebic'",
                                                                          numericInput("huge_ebic_gamma", "EBIC gamma", value = 0.5, min = 0, max = 1, step = 0.05),
                                                                          helpText("Higher gamma = sparser network")
                                                           ),
                                                           conditionalPanel("input.huge_criterion == 'ric'",
                                                                          numericInput("huge_ric_num", "Number of rotations", value = 10, min = 1, max = 100, step = 1),
                                                                          helpText("Number of random rotations for RIC")
                                                           ),
                                                           conditionalPanel("input.huge_criterion == 'stars'",
                                                                          numericInput("huge_stars_thresh", "StARS threshold", value = 0.1, min = 0.01, max = 0.5, step = 0.01),
                                                                          helpText("Variability threshold (default 0.1, conservative 0.05)"),
                                                                          numericInput("huge_stars_subsample", "Subsample ratio", value = 10, min = 2, max = 100, step = 1),
                                                                          helpText("Number of subsamples for stability selection")
                                                           )
                                          ),

                                          conditionalPanel("input.estimator == 'IsingFit'",
                                                           numericInput("ising_gamma", "EBIC gamma", value = 0.25, min = 0, max = 1, step = 0.05),
                                                           helpText("Hyperparameter for EBIC. Higher = sparser network"),
                                                           checkboxInput("ising_and", "AND-rule", value = TRUE),
                                                           helpText("If TRUE, edge is only included if both regression models agree"),
                                                           selectInput("ising_split", "Data binarization",
                                                                      choices = c("None (data already binary)" = "none",
                                                                                 "Median split" = "median",
                                                                                 "Mean split" = "mean"),
                                                                      selected = "none"),
                                                           helpText("How to convert non-binary data to binary. Choose 'None' if your data is already 0/1."),
                                                           checkboxInput("ising_plot", "Plot progress", value = FALSE),
                                                           helpText("Show estimation progress (slower)"),
                                                           helpText("ℹ️ IsingFit estimates networks for binary (0/1) data using logistic regression with L1 regularization")
                                          ),

                                          conditionalPanel("input.estimator == 'mgm'",
                                                           numericInput("order", "Interaction order", value = 2, min = 1),
                                                           helpText("Order 1 = pairwise only, Order 2 = pairwise + 2-way interactions, etc."),
                                                           shinyWidgets::pickerInput("mgm_rule", "Rule",
                                                                                     choices = c("AND","OR"), selected = "AND"),
                                                           helpText("AND = edge if all node-specific models agree; OR = edge if any agrees"),
                                                           selectInput("mgm_criterion", "Selection criterion",
                                                                      choices = c("EBIC (Extended BIC)" = "EBIC",
                                                                                "CV (Cross-Validation)" = "CV",
                                                                                "RIC (Rotation Info Criterion)" = "RIC",
                                                                                "STARS (Stability Selection)" = "STARS"),
                                                                      selected = "EBIC"),
                                                           helpText("Note: RIC/STARS use EBIC with stability-based selection"),
                                                           conditionalPanel("input.mgm_criterion == 'EBIC'",
                                                                          numericInput("tuning_mgm", "EBIC gamma", value = 0.25, min = 0, max = 1, step = 0.05),
                                                                          helpText("Higher gamma = sparser network (more conservative)")
                                                           ),
                                                           conditionalPanel("input.mgm_criterion == 'CV'",
                                                                          numericInput("mgm_cv_folds", "CV folds", value = 10, min = 2, max = 20, step = 1),
                                                                          helpText("Number of folds for k-fold cross-validation")
                                                           )
                                          ),
                                          
                                          conditionalPanel("input.estimator == 'relimp'",
                                                           checkboxInput("relimp_norm", "Normalized", TRUE),
                                                           shinyWidgets::pickerInput("structureDefault", "Structure baseline",
                                                                                     choices = c("none","EBICglasso","pcor","mgm","cor","TMFG","ggmModSelect","LoGo",
                                                                                                 "IsingFit","IsingSampler","adalasso","huge","ncvRegularize","nodeRegresIC"),
                                                                                     selected = "none"),
                                                           numericInput("relimpEstMin", "Min threshold", value = 0, step = 0.001),
                                                           checkboxInput("relimp_smart_preprocess", "Smart preprocessing", TRUE),
                                                           helpText("💡 Auto-converts dichotomous variables (0/1, Yes/No) to numeric and removes non-numeric columns. Uncheck to keep only explicitly numeric variables."),
                                                           checkboxInput("relimp_dummy", "Dummy coding", FALSE),
                                                           checkboxInput("relimp_drop_nonfinite", "Drop non-finite", TRUE),
                                                           checkboxInput("relimp_drop_zerovar", "Drop zero-variance", TRUE),
                                                           helpText("ℹ️ Variable names with special characters (e.g., 'HADS-A') are automatically converted to safe names (e.g., 'HADS.A') for R formulas.")
                                          ),
                                          
                                          hr(),
                                          checkboxInput("center_scale", "Center & scale", TRUE),
                                          
                                          h4("Appearance"),
                                          colourpicker::colourInput("pieColorChoice", "Pie color (predictability)",
                                                                    value = "#ADD8E6",
                                                                    showColour = "both",
                                                                    palette = "square",
                                                                    returnName = FALSE),
                                          colourpicker::colourInput("nodeColor", "Node color",
                                                                    value = "#FFFFFF",
                                                                    showColour = "both",
                                                                    palette = "square",
                                                                    returnName = FALSE),
                                          colourpicker::colourInput("nodeBorderColor", "Node border color",
                                                                    value = "#8F8F8F",
                                                                    showColour = "both",
                                                                    palette = "square",
                                                                    returnName = FALSE),
                                          checkboxInput("predictability_pies", "Show predictability", TRUE),
                                          conditionalPanel("input.predictability_pies == true",
                                                           sliderInput("pieBorder", "Pie chart fill",
                                                                      min = 0, max = 1, value = 0.3, step = 0.05),
                                                           helpText("0 = full node (no pie), 0.3 = ring (default), 1 = full pie chart")
                                          ),
                                          
                                          h4("📐 Layout"),
                                          shinyWidgets::pickerInput("layout", "Type",
                                                                    choices = c("Spring (Fruchterman-Reingold)"="spring",
                                                                                "Circle"="circle",
                                                                                "Groups"="groups",
                                                                                "EGA (Community Detection)"="ega",
                                                                                "PCA (Principal Components)"="pca"), 
                                                                    selected = "spring"),
                                          conditionalPanel("input.layout == 'ega'",
                                                           checkboxInput("ega_color_communities", "Color nodes by community", TRUE),
                                                           sliderInput("ega_separation", "Community separation strength", 
                                                                       min = 1, max = 10, value = 3, step = 0.5),
                                                           helpText("💡 Higher values = more space between communities. Lower = tighter clustering. Default = 3."),
                                                           helpText("⚠️ Coloring only works when >1 community is detected. Check console messages after estimation.")
                                          ),
                                          helpText("💡 Spring: general purpose. EGA: spaces communities apart (console shows # detected). PCA: variance-based positioning."),
                                          shinyWidgets::pickerInput("theme", "Theme",
                                                                    choices = c("colorblind","classic","gray","Borkulo"), 
                                                                    selected = "classic"),
                                          checkboxInput("use_numbered_nodes", "Use numbered nodes with legend", FALSE),
                                          helpText("💡 Shows numbers on nodes (1,2,3...) with full names in legend table below plot"),
                                          conditionalPanel("!input.use_numbered_nodes",
                                                           checkboxInput("show_labels", "Show node labels", TRUE)
                                          ),
                                          conditionalPanel("input.show_labels || input.use_numbered_nodes",
                                                           checkboxInput("label_scale_equal", "Equal label scaling", TRUE),
                                                           numericInput("label_cex", "Label size", value = 4.5, step = 0.1),
                                                           checkboxInput("label_bold", "Bold labels", TRUE)
                                          ),
                                          checkboxInput("show_node_means", "Show node mean levels table", FALSE),
                                          conditionalPanel("input.show_node_means && !output.has_split",
                                                           checkboxInput("show_means_on_plot", "Display means next to node labels", FALSE),
                                                           helpText("⚠️ Only available for single networks (not split groups)")
                                          ),
                                          helpText("💡 Displays mean value for each variable below the plot"),
                                          checkboxInput("show_edge_table", "Show edge table", FALSE),
                                          conditionalPanel("input.show_edge_table",
                                                           helpText("💡 Displays all edges with weights. Works with all estimators including RIN and split analyses.")
                                          ),
                                          checkboxInput("edgeLabels", "Edge labels", FALSE),
                                          numericInput("edgeLabelCex", "Edge label size", value = 0.9, step = 0.1),
                                          numericInput("vsize", "Node size", value = 6, step = 1),
                                          checkboxInput("nodeSizeByMean", "Size by mean", TRUE),
                                          conditionalPanel("input.nodeSizeByMean",
                                                           checkboxInput("useReverseScaling", "Use reverse scaling (for change scores)", FALSE),
                                                           helpText("💡 Check when negative values indicate improvement (e.g., symptom reduction). Largest reduction → largest node."),
                                                           checkboxInput("useGlobalScaling", "Scale nodes globally across groups", FALSE),
                                                           helpText("💡 Check to make node sizes comparable between groups. Unchecked = each network scaled independently.")
                                          ),
                                          sliderInput("nodeSizeRange", "Size range", min = 2, max = 20, value = c(4,10)),
                                          numericInput("esize", "Edge thickness", value = 5, step = 1),
                                          numericInput("asize", "Arrow size", value = 3, step = 1),
                                          numericInput("edgeScaleMin", "Min edge strength", value = 0, step = 0.01),
                                          numericInput("edgeScaleMax", "Max edge strength", value = 1, step = 0.01),
                                          numericInput("minEdgeVis", "Min |edge| visible", value = 0, step = 0.01),

                                          hr(),
                                          h4("Bootstrap Edge Significance"),
                                          checkboxInput("boot_sig_enable",
                                                        "Retain significant edges only (bootstrap CI)",
                                                        FALSE),
                                          conditionalPanel("input.boot_sig_enable",
                                            helpText("Runs a nonparametric bootstrap to test each edge. Edges whose bootstrap distribution does not exclude zero are set to zero."),
                                            numericInput("boot_sig_nboots", "Bootstrap samples",
                                                         value = 1000, min = 100, max = 5000, step = 100),
                                            numericInput("boot_sig_alpha", "Alpha level",
                                                         value = 0.05, min = 0.001, max = 0.2, step = 0.005),
                                            selectInput("boot_sig_correction",
                                                        "Correction for multiple comparisons",
                                                        choices = c("None"                    = "none",
                                                                    "Bonferroni"              = "bonferroni",
                                                                    "Holm"                    = "holm",
                                                                    "Hochberg"                = "hochberg",
                                                                    "BH (False Discovery Rate)" = "BH"),
                                                        selected = "none"),
                                            actionButton("run_boot_sig",
                                                         "Run Bootstrap Significance Test",
                                                         class = "btn-warning btn-block",
                                                         icon  = icon("sync")),
                                            helpText("Note: may take 1–10 minutes depending on network size and number of samples.")
                                          ),

                                          hr(),
                                          actionButton("run", "Estimate Network",
                                                       class = "btn-primary btn-block btn-lg",
                                                       icon = icon("play")),
                                          br(),
                                          downloadButton("downloadRDS", "Download Results",
                                                         class = "btn-success btn-block"),
                                          br(),
                                          downloadButton("downloadRscript", "Export R Script",
                                                         class = "btn-info btn-block",
                                                         icon  = icon("code"))
                                      )
                         ),
                         
                         mainPanel(width = 9,
                                   tabsetPanel(
                                     tabPanel("Network",
                                              br(),
                                              conditionalPanel(
                                                condition = "output.has_split",
                                                div(class = "status-warning",
                                                    h4("Network Comparison by Split Variable"),
                                                    verbatimTextOutput("split_info")
                                                )
                                              ),
                                              plotOutput("networkPlot", height = "720px"),
                                              br(),
                                              conditionalPanel(
                                                condition = "input.use_numbered_nodes",
                                                div(class = "well",
                                                    h4("Node Legend"),
                                                    helpText("Node numbers correspond to the following variables:"),
                                                    DT::dataTableOutput("node_legend_table")
                                                )
                                              ),
                                              br(),
                                              conditionalPanel(
                                                condition = "input.show_node_means",
                                                div(class = "well",
                                                    h4("Node Mean Levels"),
                                                    helpText("Shows raw (unstandardized) mean values for each variable."),
                                                    helpText("For split-group analysis: displays per-group means with between-subjects ANOVA (F and p values). Significant p-values are color-coded."),
                                                    DT::dataTableOutput("node_means_table")
                                                )
                                              ),
                                              br(),
                                              conditionalPanel(
                                                condition = "input.show_edge_table",
                                                div(class = "well",
                                                    h4("Edge Table"),
                                                    helpText("All edges (connections) in the network with their weights."),
                                                    conditionalPanel(
                                                      "output.has_split",
                                                      helpText("For split analyses: displays edges for each group separately.")
                                                    ),
                                                    DT::dataTableOutput("edge_table")
                                                )
                                              ),
                                              br(),
                                              conditionalPanel(
                                                condition = "output.has_split && input.predictability_pies",
                                                div(class = "well",
                                                    h4("Node Predictability by Group (R²)"),
                                                    helpText("Proportion of variance explained for each node in each group:"),
                                                    DT::dataTableOutput("split_predictability_table")
                                                )
                                              ),
                                              br(),
                                              div(class = "well",
                                                  verbatimTextOutput("netSummary")
                                              )
                                     ),
                                     tabPanel("Centrality",
                                              br(),
                                              div(class = "well",
                                                  h4("Centrality Measures"),
                                                  fluidRow(
                                                    column(4,
                                                           shinyWidgets::pickerInput("cent_include", "Include", 
                                                                                     multiple = TRUE,
                                                                                     choices = CENT_ALL, 
                                                                                     selected = c("Strength", "Closeness", "Betweenness"))
                                                    ),
                                                    column(4,
                                                           shinyWidgets::pickerInput("cent_orderBy", "Order by",
                                                                                     choices = CENT_ALL, selected = "Strength")
                                                    ),
                                                    column(4,
                                                           checkboxInput("cent_decreasing", "Descending", TRUE)
                                                    )
                                                  )
                                              ),
                                              plotOutput("centPlot", height = "640px"),
                                              verbatimTextOutput("centWarn")
                                     ),
                                     tabPanel("Bridge Symptoms",
                                              br(),
                                              div(class = "well",
                                                h4("Communities"),
                                                fluidRow(
                                                  column(4,
                                                    radioButtons("bridge_comm_method", "Community source",
                                                                 choices = c("Auto-detect" = "auto", "Manual" = "manual"),
                                                                 selected = "auto", inline = TRUE)
                                                  ),
                                                  column(4,
                                                    shinyWidgets::pickerInput("bridge_auto_algo", "Algorithm",
                                                                              choices = c("Walktrap" = "walktrap",
                                                                                          "Louvain"  = "louvain",
                                                                                          "EGA"      = "ega"),
                                                                              selected = "walktrap")
                                                  ),
                                                  column(4,
                                                    conditionalPanel("input.bridge_comm_method == 'manual'",
                                                      textInput("bridge_comm_manual", "Community vector",
                                                                placeholder = "e.g. 1,1,2,2,3")
                                                    )
                                                  )
                                                ),
                                                actionButton("run_bridge", "Run Bridge Analysis",
                                                             class = "btn-primary", icon = icon("project-diagram"))
                                              ),
                                              verbatimTextOutput("bridge_comm_info"),
                                              plotOutput("bridge_plot", height = "480px"),
                                              br(),
                                              DT::DTOutput("bridge_table"),
                                              hr(),
                                              div(class = "well",
                                                h4("Bridge Stability (Case-Dropping Bootstrap)"),
                                                fluidRow(
                                                  column(4,
                                                    numericInput("bridge_boots", "Bootstraps", value = 1000, min = 100, step = 100)
                                                  ),
                                                  column(4,
                                                    numericInput("bridge_ncores", "Cores", value = 4, min = 1, max = 16)
                                                  ),
                                                  column(4,
                                                    br(),
                                                    actionButton("run_bridge_stability", "Run Stability",
                                                                 class = "btn-warning", icon = icon("spinner"))
                                                  )
                                                )
                                              ),
                                              plotOutput("bridge_stability_plot", height = "400px"),
                                              verbatimTextOutput("bridge_cs_output")
                                     ),
                                     tabPanel("Network Comparison",
                                              br(),
                                              conditionalPanel(
                                                condition = "!output.has_split",
                                                div(class = "well",
                                                    h4("ℹ️ Network Comparison Test"),
                                                    helpText("Network Comparison Test (NCT) is only available when a split variable is selected."),
                                                    helpText("Please select a split variable in the 'Data Input' tab to compare networks across groups.")
                                                )
                                              ),
                                              conditionalPanel(
                                                condition = "output.has_split",
                                                div(class = "well",
                                                    h4("⚙️ Network Comparison Test Settings"),
                                                    fluidRow(
                                                      column(3,
                                                             numericInput("nct_iterations", "Number of permutations",
                                                                         value = 1000, min = 100, max = 10000, step = 100),
                                                             helpText("Default: 1000. More iterations = more accurate p-values but slower")
                                                      ),
                                                      column(3,
                                                             selectInput("nct_test_type", "Comparison type",
                                                                        choices = c("Independent (unpaired)" = "independent",
                                                                                  "Paired (same subjects)" = "paired"),
                                                                        selected = "independent"),
                                                             conditionalPanel(
                                                               "input.nct_test_type == 'paired'",
                                                               uiOutput("nct_subject_id_selector"),
                                                               helpText(HTML("<strong>ℹ️ Subject ID column:</strong><br/>
                                                                             Select to auto-match subjects.<br/>
                                                                             Leave blank if data pre-sorted."))
                                                             )
                                                      ),
                                                      column(3,
                                                             checkboxInput("nct_test_edges", "Test individual edges", value = FALSE),
                                                             helpText("Compare each edge between networks (slower)")
                                                      ),
                                                      column(3,
                                                             checkboxInput("nct_test_centrality", "Test centrality measures", value = FALSE),
                                                             helpText("Uses centrality measures from Centrality tab (Strength, Closeness, Betweenness, ExpectedInfluence)")
                                                      )
                                                    ),
                                                    fluidRow(
                                                      column(12,
                                                             div(style = "margin-top: 10px; padding: 10px; background-color: #e3f2fd; border-left: 4px solid #2196F3;",
                                                                 strong("NCT Compatible Estimators:"), br(),
                                                                 "✓ EBICglasso (Gaussian Graphical Models)", br(),
                                                                 "✓ IsingFit (Binary/Ising networks)", br(),
                                                                 "✓ Correlation & Partial Correlation", br(),
                                                                 br(),
                                                                 strong("⚠️ Not compatible with:"), br(),
                                                                 "✗ Relative Importance (RIN)", br(),
                                                                 "✗ Mixed Graphical Models (MGM)", br(),
                                                                 "✗ Huge estimator", br(),
                                                                 br(),
                                                                 strong("Parameters: "),
                                                                 "NCT will use your gamma/tuning parameters from network estimation"
                                                             )
                                                      )
                                                    ),
                                                    actionButton("run_nct", "Run Network Comparison Test",
                                                                class = "btn-primary btn-lg", icon = icon("chart-line"))
                                                ),
                                                br(),
                                                conditionalPanel(
                                                  condition = "output.nct_available",
                                                  div(class = "well",
                                                      h4("📊 Network Comparison Results"),
                                                      verbatimTextOutput("nct_summary"),
                                                      hr(),
                                                      h5("Global Strength Test"),
                                                      tableOutput("nct_global_table"),
                                                      hr(),
                                                      conditionalPanel(
                                                        condition = "input.nct_test_edges",
                                                        h5("Edge Differences"),
                                                        DT::dataTableOutput("nct_edges_table"),
                                                        hr()
                                                      ),
                                                      conditionalPanel(
                                                        condition = "input.nct_test_centrality",
                                                        h5("Centrality Differences"),
                                                        DT::dataTableOutput("nct_centrality_table")
                                                      )
                                                  )
                                                )
                                              )
                                     ),
                                     tabPanel("Data",
                                              br(),
                                              div(class = "well",
                                                  h4("Active Dataset"),
                                                  div(class = "table-responsive",
                                                      tableOutput("headData")
                                                  )
                                              ),
                                              verbatimTextOutput("dataInfo")
                                     )
                                   )
                         )
                       )
              ),

              # ---------------- NSON ----------------
              tabPanel("NSON",
                       icon = icon("project-diagram"),
                       br(),
                       div(class = "status-warning",
                           h4(HTML("&#9888; Conceptual Warning")),
                           p(HTML("<strong>NSON arrows indicate nested specificity, not causality.</strong> An arrow X &rarr; Y means GNSS(Y) &gt; GNSS(X): Y occupies a more nested/specific position in the network. This does <em>not</em> imply X causes Y, precedes Y in time, or that intervening on X will change Y.")),
                           tags$details(
                             tags$summary(
                               style = "cursor: pointer; font-weight: bold; color: #E67E22; margin-top: 4px; user-select: none;",
                               HTML("&#43; How GNSS is computed &amp; interpretation (click to expand)")
                             ),
                             br(),
                             h5("Binary data"),
                             p(HTML("For each variable Y, GNSS(Y) is the average conditional activation of all other variables X given that Y is active:")),
                             tags$blockquote(tags$code(HTML(
                               "GNSS(Y) = 1/(p&minus;1) &times; &sum;<sub>X&ne;Y</sub> P(X=1 | Y=1)"
                             ))),
                             p(HTML("where P(X=1&nbsp;|&nbsp;Y=1) is the proportion of observations where X=1 among those where Y=1.")),
                             p(HTML("A <strong>high GNSS</strong> means that when Y is active, many other variables also tend to be active — Y occupies a more nested/specific position and appears at the <em>receiving end</em> of arrows.")),
                             p(HTML("A <strong>low GNSS</strong> means Y is more broadly/generally distributed — it tends to appear at the <em>origin</em> of arrows.")),
                             hr(),
                             h5("Continuous data"),
                             p("Continuous variables are first z-score standardized:"),
                             tags$blockquote(tags$code(HTML(
                               "Z<sub>ij</sub> = (x<sub>ij</sub> &minus; &mu;<sub>j</sub>) / &sigma;<sub>j</sub>"
                             ))),
                             p(HTML("where x<sub>ij</sub> is participant i's score on variable j, &mu;<sub>j</sub> the mean, and &sigma;<sub>j</sub> the standard deviation.")),
                             p(HTML("Standardized values are then binarized with threshold &tau; to create an activation matrix A:")),
                             tags$blockquote(
                               tags$code(HTML("A<sub>ij</sub> = 1 &nbsp; if Z<sub>ij</sub> &gt; &tau;")), tags$br(),
                               tags$code(HTML("A<sub>ij</sub> = 0 &nbsp; if Z<sub>ij</sub> &le; &tau;"))
                             ),
                             p(HTML("Default &tau;&nbsp;=&nbsp;0: values above the sample mean are considered active. GNSS is then computed on A exactly as in the binary case. The threshold &tau; can be adjusted in the settings panel.")),
                             hr(),
                             h5("Interpretation of arrows"),
                             tags$blockquote(tags$code(HTML(
                               "X &rarr; Y &nbsp;&nbsp; means &nbsp;&nbsp; GNSS(Y) &gt; GNSS(X)"
                             ))),
                             p(HTML("Y occupies a more nested/specific position than X within the estimated conditional dependency structure. In practical terms, individuals with elevated Y tend to also show elevated X more often than the reverse.")),
                             hr(),
                             h5("Activation rate"),
                             tags$ul(
                               tags$li(HTML("<strong>Binary:</strong> ActivationRate(Y) = P(Y=1) — the prevalence of active observations.")),
                               tags$li(HTML("<strong>Continuous:</strong> ActivationRate(Y) = P(Z<sub>Y</sub> &gt; &tau;) — proportion of observations exceeding the z-threshold. At &tau;&nbsp;=&nbsp;0, approximately the proportion scoring above the sample mean."))
                             )
                           )
                       ),
                       br(),
                       sidebarLayout(
                         sidebarPanel(width = 3,
                                      style = "position: sticky; top: 0; max-height: 100vh; overflow-y: auto;",
                                      div(class = "well",
                                          h4("NSON Settings"),
                                          numericInput("nson_z_threshold",
                                                       "z-threshold (continuous data)",
                                                       value = 0, min = -3, max = 3, step = 0.5),
                                          helpText("Threshold for converting continuous z-scores to activation. Common values: 0, 0.5, 1."),
                                          numericInput("nson_tolerance",
                                                       "GNSS tolerance (ε)",
                                                       value = 0, min = 0, max = 0.2, step = 0.01),
                                          helpText("Edges with |GNSS difference| ≤ ε are left undirected."),
                                          shinyWidgets::pickerInput(
                                            "nson_orient_by", "Orient arrows by",
                                            choices = c(
                                              "GNSS (nested specificity)" = "gnss",
                                              "Activation rate"           = "activation_rate"),
                                            selected = "gnss"),
                                          checkboxInput("nson_show_gnss_labels",
                                                        "Show GNSS values as node labels", FALSE),
                                          checkboxInput("nson_export_activation",
                                                        "Enable activation matrix export", FALSE),
                                          hr(),
                                          h4("Layout & Theme"),
                                          shinyWidgets::pickerInput(
                                            "nson_layout", "Layout type",
                                            choices = c(
                                              "Spring (Fruchterman-Reingold)" = "spring",
                                              "Circle"                        = "circle",
                                              "EGA (Community Detection)"     = "ega",
                                              "PCA (Principal Components)"    = "pca",
                                              "Sync with main view"           = "main"),
                                            selected = "spring"),
                                          conditionalPanel(
                                            "input.nson_layout == 'ega'",
                                            sliderInput("nson_ega_separation",
                                                        "Community separation",
                                                        min = 1, max = 10,
                                                        value = 3, step = 0.5),
                                            helpText("Higher = more space between communities.")
                                          ),
                                          shinyWidgets::pickerInput(
                                            "nson_theme", "Theme",
                                            choices = c("classic", "colorblind",
                                                        "gray", "Borkulo", "Smurf"),
                                            selected = "classic"),
                                          hr(),
                                          h4("Nodes"),
                                          colourpicker::colourInput(
                                            "nson_node_color", "Node color",
                                            value = "#FFFFFF",
                                            showColour = "both", palette = "square",
                                            returnName = FALSE),
                                          colourpicker::colourInput(
                                            "nson_border_color", "Node border color",
                                            value = "#8F8F8F",
                                            showColour = "both", palette = "square",
                                            returnName = FALSE),
                                          numericInput("nson_vsize", "Node size",
                                                       value = 6, min = 1, max = 30, step = 0.5),
                                          helpText("Same scale as main network tab."),
                                          checkboxInput("nson_size_by_mean",
                                                        "Size nodes by mean", FALSE),
                                          conditionalPanel(
                                            "input.nson_size_by_mean",
                                            sliderInput("nson_size_range", "Size range",
                                                        min = 2, max = 20, value = c(4, 10)),
                                            checkboxInput("nson_reverse_scaling",
                                                          "Reverse scaling (change scores)", FALSE)
                                          ),
                                          hr(),
                                          h4("Edges & Arrows"),
                                          numericInput("nson_esize", "Edge thickness",
                                                       value = 3, min = 0.1, max = 20, step = 0.5),
                                          numericInput("nson_asize", "Arrow size",
                                                       value = 0.5, min = 0.05, max = 3, step = 0.1),
                                          helpText("Arrow size applies to the NSON directed plot only."),
                                          numericInput("nson_min_edge", "Min |edge| visible",
                                                       value = 0, min = 0, step = 0.01),
                                          numericInput("nson_edge_scale_min",
                                                       "Edge rescale min", value = 0,
                                                       min = 0, step = 0.01),
                                          numericInput("nson_edge_scale_max",
                                                       "Edge rescale max", value = 1,
                                                       min = 0, step = 0.01),
                                          hr(),
                                          h4("Labels"),
                                          checkboxInput("nson_show_labels",
                                                        "Show node labels", TRUE),
                                          conditionalPanel(
                                            "input.nson_show_labels || input.nson_show_gnss_labels",
                                            numericInput("nson_label_cex", "Label size",
                                                         value = 4.5, min = 0.5, max = 10,
                                                         step = 0.1),
                                            checkboxInput("nson_label_bold",
                                                          "Bold labels", TRUE)
                                          ),
                                          hr(),
                                          h4("Plot Zoom"),
                                          sliderInput("nson_zoom", NULL,
                                                      min = 50, max = 250, value = 100,
                                                      step = 10, post = "%"),
                                          hr(),
                                          actionButton("run_nson", "Compute NSON",
                                                       class = "btn-primary btn-block btn-lg",
                                                       icon = icon("play")),
                                          br(),
                                          downloadButton("download_nson_edges",
                                                         "Export Edge Table",
                                                         class = "btn-success btn-block"),
                                          br(),
                                          conditionalPanel(
                                            "input.nson_export_activation",
                                            downloadButton("download_nson_activation",
                                                           "Export Activation Matrix",
                                                           class = "btn-info btn-block")
                                          )
                                      )
                         ),
                         mainPanel(width = 9,
                                   tabsetPanel(
                                     tabPanel("GNSS Table",
                                              br(),
                                              div(class = "well",
                                                  h4("Global Nested Specificity Scores"),
                                                  helpText(HTML(
                                                    "GNSS (binary) or zGNSS (continuous): average conditional activation of all other nodes when this node is active/elevated.<br/>
                                                    <strong>Higher score = more nested/specific position in the network.</strong>"
                                                  )),
                                                  DT::DTOutput("nson_gnss_table")
                                              )
                                     ),
                                     tabPanel("Classic Network",
                                              br(),
                                              helpText("Standard undirected network — same edges as the NSON plot but without arrowheads."),
                                              uiOutput("nson_classic_plot_ui")
                                     ),
                                     tabPanel("NSON Plot",
                                              br(),
                                              shinyWidgets::radioGroupButtons(
                                                "nson_plot_mode", "Arrow source",
                                                choices = c("GNSS (nested specificity)" = "gnss",
                                                            "Conditional Probability"   = "cp"),
                                                selected = "gnss", status = "primary",
                                                size = "sm"
                                              ),
                                              conditionalPanel("input.nson_plot_mode === 'cp'",
                                                br(),
                                                div(class = "well",
                                                  fluidRow(
                                                    column(3, shinyWidgets::pickerInput(
                                                      "cp_mode", "Graph mode",
                                                      choices = c("Full (all pairs)"            = "full",
                                                                  "Sparse (network edges only)" = "sparse"),
                                                      selected = "sparse")),
                                                    column(3, sliderInput("cp_min_prob",
                                                      "Min probability threshold",
                                                      min = 0, max = 1, value = 0.1, step = 0.01)),
                                                    column(2, sliderInput("cp_curvature",
                                                      "Edge curvature",
                                                      min = 0, max = 1, value = 0.35, step = 0.05)),
                                                    column(2,
                                                      checkboxInput("cp_edge_labels",
                                                        "Show probability on arrows", FALSE),
                                                      checkboxInput("cp_reverse_arrows",
                                                        HTML("Reverse arrows<br/><small style='color:#888'>(X&rarr;Y shows P(X|Y))</small>"), FALSE)),
                                                    column(2, br(),
                                                      actionButton("run_cp_graph", "Compute",
                                                        class = "btn-primary btn-block",
                                                        icon  = icon("share-alt")))
                                                  )
                                                )
                                              ),
                                              uiOutput("nson_directed_plot_ui"),
                                              br(),
                                              conditionalPanel("input.nson_plot_mode !== 'cp'",
                                                div(class = "info-box",
                                                  h5("Interpretation"),
                                                  tags$ul(
                                                    tags$li("Edges represent conditional dependencies estimated by the selected network model."),
                                                    tags$li("Edge width represents the magnitude of the estimated edge weight."),
                                                    tags$li("Edge color represents the sign of the edge (green = positive, red = negative)."),
                                                    tags$li("Arrowheads represent nested specificity orientation based on GNSS/zGNSS."),
                                                    tags$li("Arrows do not represent causal effects, temporal order, or intervention effects.")
                                                  )
                                                )
                                              ),
                                              conditionalPanel(
                                                "input.nson_plot_mode === 'cp' && !input.cp_reverse_arrows",
                                                div(class = "info-box",
                                                  h5("Interpretation — standard direction"),
                                                  tags$ul(
                                                    tags$li(HTML("<strong>Arrow X &rarr; Y = P(Y=1 | X=1)</strong>: given that X is active, how likely is Y to be active? The arrow points <em>from the condition to the outcome</em>.")),
                                                    tags$li(HTML("Both X&rarr;Y and Y&rarr;X are drawn; their widths encode P(Y|X) and P(X|Y) respectively, so the two arcs in a pair can differ.")),
                                                    tags$li(HTML("<strong>Do not read arrows as causal.</strong> A thick X&rarr;Y only means X and Y frequently co-occur, with X active — it says nothing about what would happen if X were intervened on."))
                                                  )
                                                )
                                              ),
                                              conditionalPanel(
                                                "input.nson_plot_mode === 'cp' && input.cp_reverse_arrows",
                                                div(class = "info-box",
                                                  h5("Interpretation — reversed direction"),
                                                  tags$ul(
                                                    tags$li(HTML("<strong>Arrow X &rarr; Y = P(X=1 | Y=1)</strong>: given that Y is active, how likely is X to be active? The arrow points <em>from the outcome back to its condition</em>.")),
                                                    tags$li(HTML("Both X&rarr;Y and Y&rarr;X are drawn; their widths encode P(X|Y) and P(Y|X) respectively.")),
                                                    tags$li(HTML("<strong>Warning:</strong> reversed arrows are easy to misread as causal. A thick X&rarr;Y only means Y is a reliable context for X — it does <em>not</em> imply Y causes or precedes X."))
                                                  )
                                                )
                                              )
                                     ),
                                     tabPanel("Edge Table",
                                              br(),
                                              div(class = "well",
                                                  h4("Oriented Edge Table"),
                                                  helpText("All network edges with GNSS-based orientation. For oriented edges, 'from' has lower GNSS (less specific) and 'to' has higher GNSS (more specific)."),
                                                  DT::DTOutput("nson_edge_table")
                                              )
                                     ),
                                     tabPanel("Conditional Probability Graph",
                                              icon = icon("share-alt"),
                                              br(),
                                              div(class = "status-warning",
                                                h4(HTML("&#9432; Conditional Probability Graph")),
                                                p(HTML("The graph is rendered in the <strong>NSON Plot</strong> tab — switch the <em>Arrow source</em> toggle to <em>Conditional Probability</em> to view it.")),
                                                hr(),
                                                h5("Arrow direction"),
                                                tags$ul(
                                                  tags$li(HTML("<strong>Standard (default):</strong> X &rarr; Y = <strong>P(Y=1 | X=1)</strong> — the arrow points from the condition to the outcome. A thick X&rarr;Y means Y is likely active when X is active.")),
                                                  tags$li(HTML("<strong>Reversed (checkbox):</strong> X &rarr; Y = <strong>P(X=1 | Y=1)</strong> — the arrow points from the outcome back to the condition. A thick X&rarr;Y means X is likely active when Y is active."))
                                                ),
                                                p(HTML("<strong>In both modes</strong>, X&rarr;Y and Y&rarr;X are drawn as separate arcs with independent widths, reflecting the asymmetry of conditional probability. <strong>Full</strong> mode shows all variable pairs; <strong>Sparse</strong> restricts to pairs connected in the estimated network.")),
                                                div(class = "alert alert-warning", style = "margin-top:8px;",
                                                  HTML("<strong>&#9888; Not causal.</strong> Arrow width encodes co-occurrence frequency under a specific conditioning event — it does not imply that activating one node will change another, nor does it encode temporal order or intervention effects.")
                                                )
                                              ),
                                              br(),
                                              downloadButton("download_cp_matrix",
                                                             "Export CP Matrix",
                                                             class = "btn-success")
                                     ),
                                     tabPanel("Syndromic Depth & Centrality (k-way)",
                                              icon = icon("layer-group"),
                                              br(),
                                              div(class = "status-warning",
                                                h4(HTML("&#9888; Syndromic Depth &amp; Centrality (k-way)")),
                                                p(HTML("Expands pairwise Nested Specificity (k=2) into higher-order combinatorial spaces (k=3, k=4, &hellip;). Each chart is sorted by descending k=2 AND score to preserve a uniform node ranking.")),
                                                tags$details(
                                                  tags$summary(
                                                    style = "cursor:pointer;font-weight:bold;color:#E67E22;margin-top:4px;user-select:none;",
                                                    HTML("&#43; Conceptual Warning &amp; Methodological Interpretation (click to expand)")
                                                  ),
                                                  br(),
                                                  tags$ol(
                                                    tags$li(HTML("<strong>AND Rule (Syndrome Severity):</strong> measures a symptom&rsquo;s power to anchor tight overlapping clusters. Acts as an empirical proxy for a Networked Rasch Model — symptoms with tall columns retain high joint probability even under extreme co-activation demands, marking items of high psychometric &ldquo;difficulty&rdquo; and severe clinical stratum.")),
                                                    br(),
                                                    tags$li(HTML("<strong>OR Rule (Systemic Vulnerability):</strong> maps divergent causal pathways. High scores indicate a gateway node that virtually guarantees systemic comorbidity, even if the exact presentation varies between patients."))
                                                  ),
                                                  br(),
                                                  div(style = "background:#FDF2E9;border-left:4px solid #E67E22;padding:10px;border-radius:4px;",
                                                    HTML("<strong>&#9888; WARNING ON INTERPRETATION:</strong> As k increases under the AND rule, probabilities naturally compress toward zero due to the &ldquo;Curse of Dimensionality&rdquo; (the rare-event trap). If your sample size is small or the Ising network skeleton is highly sparse, deep k-way combinations may yield zero-variance rows or collapse into NA values. Use higher k settings carefully relative to your sample power.")
                                                  )
                                                )
                                              ),
                                              br(),
                                              fluidRow(
                                                column(3,
                                                  div(class = "well",
                                                    h4("k-way Settings"),
                                                    uiOutput("kway_max_k_ui"),
                                                    uiOutput("kway_group_selector_ui"),
                                                    shinyWidgets::pickerInput(
                                                      "kway_chart_choice", "Display chart",
                                                      choices = c(
                                                        "AND — Syndromic Decay (line)"              = "and_grouped",
                                                        "AND — Total Syndromic Footprint (stacked)" = "and_stacked",
                                                        "OR — Vulnerability Spread (line)"          = "or_grouped",
                                                        "OR — Total Comorbidity Footprint (stacked)" = "or_stacked"
                                                      ),
                                                      selected = "and_grouped"
                                                    ),
                                                    shinyWidgets::pickerInput(
                                                      "kway_sort_by", "Sort nodes by",
                                                      choices = c(
                                                        "k=2 AND score (desc)"     = "k2_and",
                                                        "k=2 OR score (desc)"      = "k2_or",
                                                        "Max-k AND score (desc)"   = "kmax_and",
                                                        "Max-k OR score (desc)"    = "kmax_or",
                                                        "Alphabetical"             = "alpha"
                                                      ),
                                                      selected = "k2_and"
                                                    ),
                                                    sliderInput("kway_text_size", "Label size",
                                                                min = 5, max = 16, value = 9, step = 1),
                                                    hr(),
                                                    actionButton("run_kway", "Compute k-way GNSS",
                                                                 class = "btn-primary btn-block btn-lg",
                                                                 icon  = icon("layer-group")),
                                                    br(),
                                                    downloadButton("download_kway", "Export Table",
                                                                   class = "btn-success btn-block")
                                                  )
                                                ),
                                                column(9,
                                                  plotOutput("kway_main_plot", height = "550px")
                                                )
                                              )
                                     )
                                   )
                         )
                       )
              )
  )
)

# ============================================================
# SERVER (unchanged logic, only uses enhanced UI)
# ============================================================
server <- function(input, output, session){

  # Shared layout: updated by the main network plot, consumed by NSON when "main" layout selected
  main_net_layout <- reactiveVal(NULL)

  # All server logic remains identical to original
  # [Previous server code goes here - truncated for brevity]
  # The server logic is identical to your original code
  
  output$sheet_picker <- renderUI({
    f <- input$file; if (is.null(f)) return(NULL)
    ext <- tolower(tools::file_ext(f$name))
    if (ext %in% c("xlsx","xls")) {
      sheets <- tryCatch(readxl::excel_sheets(f$datapath), error = function(e) NULL)
      if (is.null(sheets)) return(helpText("Could not read Excel sheets."))
      selectInput("sheet", "Excel sheet", choices = sheets, selected = sheets[1])
    } else NULL
  })
  
  observeEvent(input$clear_data, {
    updateTabsetPanel(session, "tabs", selected = "Data Input")
    updateCheckboxInput(session, "show_advanced", value = FALSE)
    updateRadioButtons(session, "missing_policy", selected = "listwise")
    updateCheckboxInput(session, "show_preview", value = TRUE)
    raw_cache(NULL); dat_cache(NULL); import_msgs(c("Reset done."))
    values$data <- NULL; values$available <- character(0)
    values$continuous <- character(0); values$discrete <- character(0); values$count <- character(0)
    values$split_variable <- NULL
    updateSelectInput(session, "available_vars",  choices = character(0))
    updateSelectInput(session, "continuous_vars", choices = character(0))
    updateSelectInput(session, "discrete_vars",   choices = character(0))
    updateSelectInput(session, "count_vars",      choices = character(0))
    updateSelectInput(session, "split_var", choices = c("None" = ""))
  })
  
  na_tokens_vec <- reactive({
    toks <- unlist(strsplit(input$na_tokens %||% "NA,NaN,.", ",")); toks <- trim_ws(toks); toks[nzchar(toks)]
  })
  raw_cache <- reactiveVal(NULL)
  import_msgs <- reactiveVal(character(0))
  
  raw_dat <- reactive({
    f <- input$file; msgs <- character(0)
    if (is.null(f)) { import_msgs(c("No file uploaded yet.")); return(NULL) }
    ext <- tolower(tools::file_ext(f$name))
    enc <- input$encoding %||% "UTF-8"
    dec <- if (identical(input$decimal_mark, ",")) "," else "."
    locale_obj <- readr::locale(decimal_mark = dec, encoding = enc)
    na_vec <- na_tokens_vec()

    df <- NULL
    if (ext %in% c("xlsx","xls")) {
      sh <- input$sheet %||% 1
      df <- tryCatch(readxl::read_excel(f$datapath, sheet = sh, na = na_vec),
                     error = function(e){
                       msg <- paste0("Excel read error: ", e$message)
                       showNotification(msg, type = "error", duration = 15)
                       msgs <<- c(msgs, msg); NULL
                     })
    } else if (ext %in% c("csv","tsv","txt")) {
      delim <- input$txt_delim %||% "auto"
      if (identical(delim, "auto")) {
        df <- tryCatch(vroom::vroom(f$datapath, na = na_vec, col_names = isTRUE(input$has_header),
                                    locale = vroom::vroom_locale(decimal_mark = dec, encoding = enc),
                                    altrep = FALSE),
                       error = function(e){
                         msgs <<- c(msgs, paste0("vroom auto-delim failed: ", e$message, " — trying readr..."))
                         tryCatch(readr::read_delim(f$datapath, delim = ",", na = na_vec,
                                                    col_names = isTRUE(input$has_header), locale = locale_obj),
                                  error = function(e2){ msgs <<- c(msgs, paste0("readr failed: ", e2$message)); NULL })
                       })
      } else {
        df <- tryCatch(readr::read_delim(f$datapath, delim = delim, na = na_vec,
                                         col_names = isTRUE(input$has_header), locale = locale_obj,
                                         guess_max = 10000, progress = FALSE),
                       error = function(e){ msgs <<- c(msgs, paste0("readr failed: ", e$message)); NULL })
      }
    } else { msgs <- c(msgs, paste0("Unsupported extension: .", ext)) }
    
    if (is.null(df)) {
      final_msg <- if (length(msgs) > 0) paste(msgs, collapse="\n") else "File read failed."
      showNotification(final_msg, type = "error", duration = 20)
      import_msgs(final_msg); return(NULL)
    }
    df <- as.data.frame(df, check.names = TRUE, stringsAsFactors = FALSE)
    msgs <- c(msgs, paste0("✓ Loaded: ", nrow(df), " rows × ", ncol(df), " columns"))
    if (ncol(df) < 2) msgs <- c(msgs, "⚠ Warning: fewer than 2 columns")

    # Auto-fix columns where values are numeric + trailing letter (e.g. "2b", "1a").
    # This handles instruments like the Beck Depression Inventory where items 16/18
    # have a/b sub-versions encoded as "2b", "1a", etc.
    # Columns that are purely non-numeric letters (e.g. "a"/"b" version indicators)
    # are left unchanged and will be treated as categorical.
    fixed_cols <- character(0)
    for (col in colnames(df)) {
      x <- df[[col]]
      if (!is.character(x)) next
      non_na <- x[!is.na(x) & nchar(trimws(x)) > 0]
      if (length(non_na) == 0) next
      # Match values that are numeric optionally followed by a single letter (e.g. "2b", "1a", "3")
      is_num_letter <- grepl("^-?[0-9]+\\.?[0-9]*[a-zA-Z]?$", trimws(non_na))
      has_letter    <- grepl("^-?[0-9]+\\.?[0-9]*[a-zA-Z]$",  trimws(non_na))
      if (all(is_num_letter) && any(has_letter)) {
        df[[col]] <- suppressWarnings(as.numeric(sub("[a-zA-Z]+$", "", trimws(x))))
        fixed_cols <- c(fixed_cols, col)
      }
    }
    if (length(fixed_cols) > 0) {
      msgs <- c(msgs, paste0("ℹ Auto-stripped letter suffixes in ",
                             length(fixed_cols), " column(s) (e.g. '2b'→2): ",
                             paste(fixed_cols, collapse = ", ")))
    }

    raw_cache(df); import_msgs(msgs); df
  })
  
  dat_cache <- reactiveVal(NULL)
  dat_unbalanced <- reactive({
    df <- raw_dat(); if (is.null(df)) return(NULL)
    if (identical(input$missing_policy, "listwise")) {
      cc <- stats::complete.cases(df)
      if (!any(cc)) {
        import_msgs(c(isolate(import_msgs()),
                      "❌ Listwise deletion removed all rows (0 complete cases).",
                      "   → Switch to 'pairwise' in Import Options > Show advanced options."))
        dat_cache(NULL); return(NULL)
      }
      out <- df[cc, , drop = FALSE]; dat_cache(out); out
    } else { dat_cache(df); df }
  })

  # Apply equal sampling if split variable is set and mode is "equal"
  dat <- reactive({
    df <- dat_unbalanced()
    if (is.null(df)) return(NULL)

    # Check if split variable is set and equal sampling is requested
    split <- split_var()
    sample_mode <- input$split_sample_mode %||% "all"

    if (!is.null(split) && split != "" && sample_mode == "equal" && split %in% colnames(df)) {
      message("\n=== Equal Sample Size Balancing ===")

      # Get group sizes (sort to ensure consistent ordering)
      split_values <- sort(unique(na.omit(df[[split]])))
      group_sizes <- table(df[[split]])

      message(paste("Split variable:", split))
      message(paste("Groups:", paste(names(group_sizes), collapse = ", ")))
      message(paste("Original sizes:", paste(group_sizes, collapse = ", ")))

      # Find minimum group size
      min_size <- min(group_sizes)
      message(paste("Target size (minimum):", min_size))

      # Get seed settings
      use_seed <- isTRUE(input$split_random_seed)
      seed_value <- as.integer(input$split_seed_value)

      # Default values if inputs not ready
      if (is.na(use_seed)) use_seed <- TRUE
      if (is.na(seed_value) || is.null(seed_value)) seed_value <- 123L

      if (use_seed) {
        message(paste("Using fixed seed:", seed_value, "(reproducible sampling)"))
      } else {
        message("Using random seed (truly random sampling - different each run)")
      }

      # Sample from each group
      balanced_indices <- c()
      for (group_val in split_values) {
        group_idx <- which(df[[split]] == group_val)

        if (length(group_idx) > min_size) {
          # Set seed if reproducible mode is on
          if (use_seed) {
            set.seed(seed_value)
          }

          # Randomly sample min_size observations
          sampled_idx <- sample(group_idx, size = min_size, replace = FALSE)
          balanced_indices <- c(balanced_indices, sampled_idx)
          message(paste("  Group", group_val, ": sampled", min_size, "from", length(group_idx)))
        } else {
          # Keep all observations
          balanced_indices <- c(balanced_indices, group_idx)
          message(paste("  Group", group_val, ": kept all", length(group_idx)))
        }
      }

      # Return balanced dataset
      df_balanced <- df[balanced_indices, , drop = FALSE]
      message(paste("Total observations after balancing:", nrow(df_balanced)))
      message("=== End Balancing ===\n")

      return(df_balanced)
    }

    # No balancing needed
    return(df)
  })
  
  output$import_status <- renderText({ paste(import_msgs(), collapse = "\n") })
  output$previewDT <- DT::renderDT({
    df <- raw_dat(); shiny::validate(shiny::need(!is.null(df), "No file loaded."))
    DT::datatable(head(df, 1000), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$preview_summary <- renderPrint({
    df <- raw_dat(); if (is.null(df)) return(cat("No data."))
    cat("Rows:", nrow(df), " Cols:", ncol(df), "\n\n")
    types <- sapply(df, function(x) paste(class(x), collapse=","))
    nas <- sapply(df, function(x) sum(is.na(x)))
    info <- data.frame(Variable = colnames(df), Type = unlist(types), NA_count = unlist(nas))
    print(info, row.names = FALSE)
  })
  
  values <- reactiveValues(data = NULL, available = character(0), continuous = character(0),
                           discrete = character(0), count = character(0), split_variable = NULL, node_legend = NULL,
                           node_means = NULL, split_predictability = NULL)
  
  observeEvent(dat(), {
    df <- dat(); if (is.null(df)) return()
    values$data <- df
    vars <- colnames(df)
    values$continuous <- intersect(values$continuous, vars)
    values$discrete   <- intersect(values$discrete, vars)
    values$count      <- intersect(values$count, vars)
    if (!is.null(values$split_variable) && !(values$split_variable %in% vars)) {
      values$split_variable <- NULL
    }
    taken <- union(union(values$continuous, values$discrete), values$count)
    if (!is.null(values$split_variable)) taken <- union(taken, values$split_variable)
    values$available  <- setdiff(vars, taken)
    updateSelectInput(session, "available_vars",  choices = values$available,  selected = NULL)
    updateSelectInput(session, "continuous_vars", choices = values$continuous, selected = NULL)
    updateSelectInput(session, "discrete_vars",   choices = values$discrete,   selected = NULL)
    updateSelectInput(session, "count_vars",      choices = values$count,      selected = NULL)
    all_vars <- colnames(df)
    updateSelectInput(session, "split_var",
                      choices = c("None" = "", all_vars),
                      selected = values$split_variable %||% "")
  }, ignoreInit = FALSE)
  
  observeEvent(input$auto_detect, {
    req(values$data)
    continuous <- c(); discrete <- c(); count <- c()
    for(col in names(values$data)) {
      if (!is.null(values$split_variable) && col == values$split_variable) next
      if (is.numeric(values$data[[col]])) {
        vals <- na.omit(values$data[[col]])
        n_unique <- length(unique(vals))

        # Check if values are all non-negative integers (potential count variable)
        is_integer_like <- all(vals == floor(vals))
        is_nonneg <- all(vals >= 0)

        if (is_integer_like && is_nonneg && n_unique <= 20) {
          # Likely a count variable (Poisson)
          count <- c(count, col)
        } else if (n_unique > 10) {
          # Continuous (Gaussian)
          continuous <- c(continuous, col)
        } else {
          # Few unique values - categorical
          discrete <- c(discrete, col)
        }
      } else {
        # Non-numeric - categorical
        discrete <- c(discrete, col)
      }
    }
    values$continuous <- unique(continuous)
    values$discrete   <- unique(discrete)
    values$count      <- unique(count)
    taken <- union(union(values$continuous, values$discrete), values$count)
    if (!is.null(values$split_variable)) taken <- union(taken, values$split_variable)
    values$available  <- setdiff(names(values$data), taken)
    updateSelectInput(session, "available_vars",  choices = values$available,  selected = NULL)
    updateSelectInput(session, "continuous_vars", choices = values$continuous, selected = NULL)
    updateSelectInput(session, "discrete_vars",   choices = values$discrete,   selected = NULL)
    updateSelectInput(session, "count_vars",      choices = values$count,      selected = NULL)
    showNotification(paste("✓ Auto-classified:", length(values$continuous), "continuous,",
                           length(values$discrete), "categorical,", length(values$count), "count"),
                     type = "message", duration = 3)
  })
  
  observeEvent(input$set_split, {
    split_val <- input$split_var
    if (is.null(split_val)) return()
    current_split <- values$split_variable
    if (is.null(current_split)) current_split <- ""
    if (split_val != "" && split_val != current_split) {
      old_split <- values$split_variable
      values$split_variable <- split_val
      values$available <- setdiff(values$available, split_val)
      values$continuous <- setdiff(values$continuous, split_val)
      values$discrete <- setdiff(values$discrete, split_val)
      values$count <- setdiff(values$count, split_val)
      if (!is.null(old_split) && old_split != "") {
        values$available <- c(values$available, old_split)
      }
      updateSelectInput(session, "available_vars", choices = values$available, selected = NULL)
      updateSelectInput(session, "continuous_vars", choices = values$continuous, selected = NULL)
      updateSelectInput(session, "discrete_vars", choices = values$discrete, selected = NULL)
      updateSelectInput(session, "count_vars", choices = values$count, selected = NULL)
      showNotification(paste("✓ Split variable:", split_val), type = "message", duration = 3)
    }
  })
  
  observeEvent(input$clear_split, {
    if (!is.null(values$split_variable) && values$split_variable != "") {
      old_split <- values$split_variable
      values$split_variable <- NULL
      values$available <- c(values$available, old_split)
      updateSelectInput(session, "available_vars", choices = values$available, selected = old_split)
      updateSelectInput(session, "split_var", selected = "")
      showNotification("✓ Split cleared", type = "message", duration = 3)
    }
  })
  
  observeEvent(input$add_continuous, {
    s <- input$available_vars; if(!is.null(s) && length(s) > 0) {
      values$continuous <- unique(c(values$continuous, s))
      values$available  <- setdiff(values$available, s)
      updateSelectInput(session, "continuous_vars", choices = values$continuous, selected = s)
      updateSelectInput(session, "available_vars",  choices = values$available,  selected = NULL)
    }
  })
  observeEvent(input$remove_continuous, {
    s <- input$continuous_vars; if(!is.null(s) && length(s) > 0) {
      values$available  <- unique(c(values$available, s))
      values$continuous <- setdiff(values$continuous, s)
      updateSelectInput(session, "continuous_vars", choices = values$continuous, selected = NULL)
      updateSelectInput(session, "available_vars",  choices = values$available,  selected = s)
    }
  })
  observeEvent(input$add_discrete, {
    s <- input$available_vars; if(!is.null(s) && length(s) > 0) {
      values$discrete  <- unique(c(values$discrete, s))
      values$available <- setdiff(values$available, s)
      updateSelectInput(session, "discrete_vars",  choices = values$discrete,  selected = s)
      updateSelectInput(session, "available_vars", choices = values$available, selected = NULL)
    }
  })
  observeEvent(input$remove_discrete, {
    s <- input$discrete_vars; if(!is.null(s) && length(s) > 0) {
      values$available <- unique(c(values$available, s))
      values$discrete  <- setdiff(values$discrete, s)
      updateSelectInput(session, "discrete_vars",  choices = values$discrete,  selected = NULL)
      updateSelectInput(session, "available_vars", choices = values$available, selected = s)
    }
  })
  observeEvent(input$add_count, {
    s <- input$available_vars; if(!is.null(s) && length(s) > 0) {
      values$count     <- unique(c(values$count, s))
      values$available <- setdiff(values$available, s)
      updateSelectInput(session, "count_vars",     choices = values$count,     selected = s)
      updateSelectInput(session, "available_vars", choices = values$available, selected = NULL)
    }
  })
  observeEvent(input$remove_count, {
    s <- input$count_vars; if(!is.null(s) && length(s) > 0) {
      values$available <- unique(c(values$available, s))
      values$count     <- setdiff(values$count, s)
      updateSelectInput(session, "count_vars",     choices = values$count,     selected = NULL)
      updateSelectInput(session, "available_vars", choices = values$available, selected = s)
    }
  })

  output$vars_status <- renderText({
    paste0("Available: ", length(values$available),
           " | Continuous: ", length(values$continuous),
           " | Categorical: ", length(values$discrete),
           " | Count: ", length(values$count))
  })
  
  output$split_status <- renderText({
    if (!is.null(values$split_variable) && values$split_variable != "") {
      df <- dat()
      if (!is.null(df) && values$split_variable %in% colnames(df)) {
        n_levels <- length(unique(na.omit(df[[values$split_variable]])))
        paste0("Split: ", values$split_variable, " (", n_levels, " groups)")
      } else {
        "Split: None"
      }
    } else {
      "Split: None"
    }
  })
  
  output$has_split <- reactive({
    !is.null(values$split_variable) && values$split_variable != ""
  })
  outputOptions(output, "has_split", suspendWhenHidden = FALSE)
  
  cont_vars <- reactive({ values$continuous })
  disc_vars <- reactive({ values$discrete })
  count_vars <- reactive({ values$count })
  split_var <- reactive({
    sv <- values$split_variable
    if (is.null(sv) || sv == "") return(NULL)
    return(sv)
  })
  
  output$headData <- renderTable({ head(dat() %||% data.frame()) })
  output$dataInfo <- renderPrint({
    df <- dat(); if (is.null(df)) return(cat("No data available."))
    str(df)
  })
  
  prepped_data <- reactive({
    df <- dat(); if (is.null(df)) return(NULL)
    keep <- unique(c(cont_vars(), disc_vars(), count_vars()))
    if (length(keep) > 1) df <- df[, intersect(keep, colnames(df)), drop = FALSE]
    if (ncol(df) < 2) return(NULL)

    # Convert categorical variables to factors
    for (nm in intersect(disc_vars(), colnames(df))) if (!is.factor(df[[nm]])) df[[nm]] <- factor(df[[nm]])

    # Center and scale continuous variables only (not count variables)
    if (isTRUE(input$center_scale)) {
      cont_idx <- colnames(df) %in% cont_vars() & sapply(df, is_numeric_like)
      if (any(cont_idx)) df[, cont_idx] <- as.data.frame(scale(df[, cont_idx, drop = FALSE]))
    }
    df
  })
  
  means_for_nodes <- reactive({
    df <- dat()
    if (is.null(df)) return(NULL)
    split <- split_var()
    if (!is.null(split) && split != "" && split %in% colnames(df)) {
      split_values <- unique(na.omit(df[[split]]))
      means_list <- list()
      for (group in split_values) {
        group_idx <- which(df[[split]] == group)
        group_df <- df[group_idx, , drop = FALSE]
        m <- vapply(group_df, function(col) {
          if (is.numeric(col)) {
            mean(col, na.rm = TRUE)
          } else if (is.logical(col)) {
            mean(col, na.rm = TRUE)
          } else if (is.factor(col) && nlevels(col) == 2) {
            v <- as.numeric(col) - 1
            mean(v, na.rm = TRUE)
          } else {
            NA_real_
          }
        }, numeric(1))
        names(m) <- colnames(group_df)
        means_list[[as.character(group)]] <- m
      }
      return(means_list)
    } else {
      m <- vapply(df, function(col) {
        if (is.numeric(col)) {
          mean(col, na.rm = TRUE)
        } else if (is.logical(col)) {
          mean(col, na.rm = TRUE)
        } else if (is.factor(col) && nlevels(col) == 2) {
          v <- as.numeric(col) - 1
          mean(v, na.rm = TRUE)
        } else {
          NA_real_
        }
      }, numeric(1))
      names(m) <- colnames(df)
      return(m)
    }
  })
  
  # ---- Bootstrap edge significance state ----
  boot_sig_result <- reactiveVal(NULL)

  # Reset bootstrap result whenever the network is re-estimated
  observeEvent(input$run, { boot_sig_result(NULL) }, ignoreInit = TRUE)

  observeEvent(input$run_boot_sig, {
    req(input$boot_sig_enable)
    res <- est_result()
    if (is.null(res)) {
      showNotification("Estimate a network first before running the bootstrap test.", type = "warning")
      return()
    }

    net <- if (res$split) {
      showNotification(
        "Split-group network detected. Bootstrap will run on the first group network.",
        type = "warning", duration = 6
      )
      res$groups[[1]]
    } else {
      res$network
    }

    if (is.null(net)) {
      showNotification("No estimated network found.", type = "error")
      return()
    }

    nBoots <- as.integer(input$boot_sig_nboots %||% 1000)

    withProgress(message = paste0("Bootstrap in progress (", nBoots, " samples)..."), value = 0.1, {
      result <- tryCatch({
        bootnet::bootnet(net,
                         nBoots     = nBoots,
                         type       = "nonparametric",
                         statistics = "edge",
                         nCores     = 1,
                         verbose    = FALSE)
      }, error = function(e) {
        showNotification(paste("Bootstrap failed:", conditionMessage(e)),
                         type = "error", duration = 10)
        NULL
      })
      setProgress(1)
    })

    boot_sig_result(result)

    if (!is.null(result)) {
      showNotification(
        paste0("Bootstrap complete (", nBoots, " samples). Network updated."),
        type = "message", duration = 5
      )
    }
  })

  est_result <- eventReactive(input$run, {
    df <- prepped_data()
    shiny::validate(shiny::need(!is.null(df), "Data not ready: ensure at least two variables."))
    split <- split_var()
    if (!is.null(split) && split != "") {
      full_df <- dat()
      if (is.null(full_df) || !(split %in% colnames(full_df))) {
        showNotification("Split variable not found", type = "error")
        return(NULL)
      }
      split_values <- sort(unique(na.omit(full_df[[split]])))
      if (length(split_values) > 10) {
        showNotification("Too many groups (max 10)", type = "error")
        return(NULL)
      }
      if (length(split_values) < 2) {
        showNotification("Need at least 2 groups", type = "error")
        return(NULL)
      }
      results <- list()
      failed_groups <- character(0)
      
      for (group in split_values) {
        group_idx <- which(full_df[[split]] == group)
        group_df <- df[group_idx, , drop = FALSE]
        
        if (nrow(group_df) < 10) {
          showNotification(paste("Group", group, "has too few observations (n =", nrow(group_df), ")"), 
                           type = "warning", duration = 4)
          failed_groups <- c(failed_groups, as.character(group))
          next
        }
        
        # Estimate network for this group with error handling
        group_net <- tryCatch({
          estimate_single_network(group_df, input)
        }, error = function(e) {
          showNotification(
            paste0("Group '", group, "' failed: ", e$message), 
            type = "warning", 
            duration = 6
          )
          NULL
        })
        
        if (!is.null(group_net)) {
          results[[as.character(group)]] <- group_net
        } else {
          failed_groups <- c(failed_groups, as.character(group))
        }
      }
      
      if (length(results) == 0) {
        showNotification("No valid groups with sufficient data", type = "error")
        return(NULL)
      }
      
      if (length(failed_groups) > 0) {
        showNotification(
          paste0("⚠️ ", length(failed_groups), " group(s) failed: ", 
                 paste(failed_groups, collapse = ", "), 
                 "\nSuccessfully estimated: ", length(results), " group(s)"), 
          type = "warning",
          duration = 8
        )
      }
      return(list(split = TRUE, groups = results, split_var = split))
    } else {
      return(list(split = FALSE, network = estimate_single_network(df, input)))
    }
  }, ignoreInit = TRUE)

  # Helper function to apply thresholding to correlation/pcor networks
  apply_threshold <- function(net, df, method, manual_value = NULL) {
    # Get the weight matrix
    W <- tryCatch(getWmat(net), error = function(e) NULL)
    if (is.null(W)) return(net)

    n <- nrow(df)
    p <- ncol(W)

    if (method == "manual") {
      # Simple manual cutoff
      threshold_val <- manual_value %||% 0.1
      W[abs(W) < threshold_val] <- 0
    } else {
      # Calculate p-values for correlations
      # Use Fisher's z-transformation for correlation significance
      r_to_p <- function(r, n) {
        if (is.na(r) || r == 0) return(1)
        t_stat <- r * sqrt(n - 2) / sqrt(1 - r^2)
        2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)
      }

      # Create p-value matrix
      pval_mat <- matrix(1, nrow = p, ncol = p)
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          pval_mat[i, j] <- pval_mat[j, i] <- r_to_p(W[i, j], n)
        }
      }

      # Apply multiple testing correction
      pvals_vec <- pval_mat[upper.tri(pval_mat)]

      if (method == "sig") {
        # No correction, just p < 0.05
        sig_threshold <- 0.05
        pvals_adj <- pvals_vec
      } else if (method %in% c("bonferroni", "holm", "hochberg", "hommel", "BH")) {
        # R's p.adjust function
        pvals_adj <- p.adjust(pvals_vec, method = tolower(method))
        sig_threshold <- 0.05
      } else {
        pvals_adj <- pvals_vec
        sig_threshold <- 0.05
      }

      # Create significance mask
      sig_mask <- matrix(TRUE, nrow = p, ncol = p)
      sig_mask[upper.tri(sig_mask)] <- pvals_adj < sig_threshold
      sig_mask[lower.tri(sig_mask)] <- t(sig_mask)[lower.tri(sig_mask)]
      diag(sig_mask) <- TRUE

      # Zero out non-significant edges
      W[!sig_mask] <- 0

      message(paste("Threshold method:", method))
      message(paste("Edges before threshold:", sum(abs(pval_mat[upper.tri(pval_mat)]) < 1)))
      message(paste("Edges after threshold:", sum(abs(W[upper.tri(W)]) > 0)))
    }

    # Update the network object with thresholded weights
    net$graph <- W
    return(net)
  }

  # Helper: bootstrap CI-based edge significance thresholding
  # Uses the empirical bootstrap distribution for each edge to compute a two-sided
  # p-value (2 * min(P(boot <= 0), P(boot >= 0))), then applies p.adjust.
  apply_boot_threshold <- function(net, boot_result, alpha = 0.05, correction = "none") {
    W <- tryCatch(bootnet::getWmat(net), error = function(e) NULL)
    if (is.null(W) || is.null(boot_result)) return(net)
    p <- ncol(W)
    if (p < 2) return(net)

    boots <- boot_result$boots
    if (is.null(boots) || length(boots) == 0) return(net)

    pairs    <- which(upper.tri(W), arr.ind = TRUE)
    n_edges  <- nrow(pairs)
    if (n_edges == 0) return(net)

    pvals <- numeric(n_edges)
    for (k in seq_len(n_edges)) {
      i <- pairs[k, 1]; j <- pairs[k, 2]
      boot_ests <- sapply(boots, function(b) {
        g <- b$graph
        if (!is.null(g) && i <= nrow(g) && j <= ncol(g)) g[i, j] else NA_real_
      })
      boot_ests <- boot_ests[!is.na(boot_ests)]
      if (length(boot_ests) == 0) { pvals[k] <- 1; next }
      # Two-sided p-value from empirical CDF at 0
      F0 <- mean(boot_ests <= 0)
      pvals[k] <- 2 * min(F0, 1 - F0)
    }

    pvals_adj <- if (correction != "none") p.adjust(pvals, method = correction) else pvals

    n_before <- sum(abs(W[upper.tri(W)]) > 0)
    for (k in seq_len(n_edges)) {
      if (pvals_adj[k] >= alpha) {
        i <- pairs[k, 1]; j <- pairs[k, 2]
        W[i, j] <- 0; W[j, i] <- 0
      }
    }
    n_after <- sum(abs(W[upper.tri(W)]) > 0)
    message(sprintf("Bootstrap threshold (%s, alpha=%.3f): %d -> %d edges",
                    correction, alpha, n_before, n_after))

    net$graph <- W
    net
  }

  # ---- Fallback debug helpers ----

  # Resolve a parameter value; record whether a fallback to default was applied
  resolve <- function(val, default_val) {
    is_fb <- is.null(val) || (is.character(val) && identical(val, ""))
    list(value = if (is_fb) default_val else val, fallback = is_fb)
  }

  # Print all resolved parameters to the R console, flagging any fallbacks
  debug_params <- function(params_list) {
    message("\n========== Effective Analysis Parameters ==========")
    any_fb <- FALSE
    for (nm in names(params_list)) {
      entry <- params_list[[nm]]
      if (isTRUE(entry$fallback)) {
        message(sprintf("  [FALLBACK] %-28s NULL/empty -> %s", nm, deparse(entry$value)))
        any_fb <- TRUE
      } else {
        message(sprintf("            %-28s %s", nm, deparse(entry$value)))
      }
    }
    if (!any_fb) message("  (no fallbacks - all parameters explicitly set)")
    message("====================================================\n")
  }

  # ---- R script code generator ----

  generate_r_code <- function(input, res, boot_res, df_prepped, df_raw) {
    def          <- input$estimator %||% "EBICglasso"
    ts           <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    split_v      <- input$split_var %||% ""
    use_split    <- !is.null(split_v) && split_v != "" && !is.null(res$split) && res$split
    cv           <- isTRUE(input$center_scale)
    n_rows       <- if (!is.null(df_prepped)) nrow(df_prepped) else "?"

    # Variable lists
    cv_vars  <- cont_vars()
    dv_vars  <- disc_vars()
    cov_vars <- count_vars()
    all_vars <- unique(c(cv_vars, dv_vars, cov_vars))

    fmt_vec <- function(v) if (length(v) == 0) "character(0)" else
      paste0('c(', paste0('"', v, '"', collapse = ", "), ')')

    lines <- character(0)
    L <- function(...) { lines <<- c(lines, sprintf(...)) }
    Lraw <- function(x) { lines <<- c(lines, x) }
    Lblank <- function() { lines <<- c(lines, "") }

    # Header
    L("# ============================================================")
    L("# BootRIN Export  |  Estimator: %s  |  %s", def, ts)
    L("# Variables: %d  |  N = %s  |  Split: %s",
      length(all_vars), n_rows, if (use_split) split_v else "none")
    L("# ============================================================")
    Lblank()

    # 1. Packages
    L("# ---- 1. Packages ----")
    L("library(bootnet)")
    L("library(qgraph)")
    if (def == "relimp")   L("library(relaimpo)")
    if (def == "mgm")      L("library(mgm)")
    if (def == "IsingFit") L("library(IsingFit)")
    if (def == "huge")     L("library(huge)")
    Lblank()

    # 2. Data
    L("# ---- 2. Load data ----")
    L("# Replace the path below with the location of your data file")
    if (!is.null(df_raw)) {
      # Try to get a filename from the session info if available
      L("# df <- read.csv(\"your_data.csv\")")
    } else {
      L("# df <- read.csv(\"your_data.csv\")")
    }
    Lblank()
    L("# Variables used in the analysis")
    if (length(cv_vars)  > 0) L("vars_continuous <- %s", fmt_vec(cv_vars))
    if (length(dv_vars)  > 0) L("vars_discrete   <- %s", fmt_vec(dv_vars))
    if (length(cov_vars) > 0) L("vars_count      <- %s", fmt_vec(cov_vars))
    L("all_vars        <- %s", fmt_vec(all_vars))
    Lblank()
    L("df_analysis <- df[, all_vars, drop = FALSE]")
    Lblank()

    # 3. Split
    if (use_split) {
      L("# ---- 3. Split groups by '%s' ----", split_v)
      grp_vals <- names(res$groups)
      for (i in seq_along(grp_vals)) {
        g  <- grp_vals[i]
        L("group_%s <- df_analysis[df[[\"%s\"]] == %s, ]",
          gsub("[^A-Za-z0-9_]", "_", g), split_v,
          if (is.numeric(suppressWarnings(as.numeric(g)))) g else paste0('"', g, '"'))
      }
      Lblank()
    }

    # 4. Preprocessing
    L("# ---- 4. Preprocessing ----")
    L("# center_scale: %s", cv)
    if (cv && length(cv_vars) > 0) {
      L("df_analysis[, vars_continuous] <- scale(df_analysis[, vars_continuous])")
      if (use_split) {
        grp_vals <- names(res$groups)
        for (g in grp_vals) {
          gn <- gsub("[^A-Za-z0-9_]", "_", g)
          L("group_%s[, vars_continuous] <- scale(group_%s[, vars_continuous])", gn, gn)
        }
      }
    }
    if (def == "relimp") {
      Lblank()
      L("# relimp: column name sanitization + keep numeric columns only")
      if (use_split) {
        grp_vals <- names(res$groups)
        for (g in grp_vals) {
          gn <- gsub("[^A-Za-z0-9_]", "_", g)
          L("colnames(group_%s) <- make.names(colnames(group_%s), unique = TRUE)", gn, gn)
          L("group_%s <- group_%s[, sapply(group_%s, is.numeric), drop = FALSE]", gn, gn, gn)
          if (isTRUE(input$relimp_drop_zerovar)) {
            L("group_%s <- group_%s[, vapply(group_%s, function(x) var(x, na.rm=TRUE) > 1e-10, logical(1))]",
              gn, gn, gn)
          }
          if (isTRUE(input$relimp_drop_nonfinite)) {
            L("group_%s <- group_%s[apply(group_%s, 1, function(r) all(is.finite(r))), ]",
              gn, gn, gn)
          }
        }
      } else {
        L("colnames(df_analysis) <- make.names(colnames(df_analysis), unique = TRUE)")
        L("df_analysis <- df_analysis[, sapply(df_analysis, is.numeric), drop = FALSE]")
        if (isTRUE(input$relimp_drop_zerovar)) {
          L("df_analysis <- df_analysis[, vapply(df_analysis, function(x) var(x, na.rm=TRUE) > 1e-10, logical(1))]")
        }
        if (isTRUE(input$relimp_drop_nonfinite)) {
          L("df_analysis <- df_analysis[apply(df_analysis, 1, function(r) all(is.finite(r))), ]")
        }
      }
    }
    Lblank()

    # 5. Estimation
    L("# ---- 5. Estimate network ----")

    est_call <- function(data_name) {
      missing_m <- input$missing_policy %||% "listwise"
      if (def == "EBICglasso") {
        cor_m  <- input$corMethod_ebic  %||% "cor"
        gamma  <- input$tuning          %||% 0.5
        refit  <- isTRUE(input$refit)
        ss     <- input$ebic_samplesize %||% "maximum"
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default     = \"EBICglasso\",")
        L("  corMethod   = \"%s\",",   cor_m)
        L("  tuning      = %s,",       gamma)
        L("  refit       = %s,",       refit)
        L("  missing     = \"%s\",",   missing_m)
        L("  sampleSize  = \"%s\",",   ss)
        L("  weighted    = TRUE,")
        L("  signed      = TRUE)")
      } else if (def == "pcor") {
        cor_m <- input$corMethod %||% "cor_auto"
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default   = \"pcor\",")
        L("  corMethod = \"%s\",", cor_m)
        L("  missing   = \"%s\")", missing_m)
      } else if (def == "cor") {
        cor_m <- input$corMethod %||% "cor_auto"
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default   = \"cor\",")
        L("  corMethod = \"%s\",", cor_m)
        L("  missing   = \"%s\")", missing_m)
        thr_m <- input$threshold_method %||% "none"
        if (thr_m != "none") {
          Lblank()
          L("# Significance thresholding: method = '%s'", thr_m)
          L("n    <- nrow(%s)", data_name)
          L("W    <- %s$graph", data_name)
          L("p    <- ncol(W)")
          if (thr_m == "manual") {
            thr_v <- input$threshold_value %||% 0.1
            L("W[abs(W) < %s] <- 0", thr_v)
          } else {
            L("r_to_p <- function(r, n) {")
            L("  if (is.na(r) || r == 0) return(1)")
            L("  t_stat <- r * sqrt(n - 2) / sqrt(1 - r^2)")
            L("  2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)")
            L("}")
            L("pval_mat <- matrix(1, p, p)")
            L("for (i in 1:(p-1)) for (j in (i+1):p)")
            L("  pval_mat[i,j] <- pval_mat[j,i] <- r_to_p(W[i,j], n)")
            L("pvals <- pval_mat[upper.tri(pval_mat)]")
            if (thr_m == "sig") {
              L("pvals_adj <- pvals")
            } else {
              L("pvals_adj <- p.adjust(pvals, method = \"%s\")", thr_m)
            }
            L("sig_mask <- matrix(TRUE, p, p)")
            L("sig_mask[upper.tri(sig_mask)] <- pvals_adj < 0.05")
            L("sig_mask[lower.tri(sig_mask)] <- t(sig_mask)[lower.tri(sig_mask)]")
            L("W[!sig_mask] <- 0")
          }
          L("%s$graph <- W", data_name)
        }
      } else if (def == "huge") {
        hm  <- input$huge_method    %||% "glasso"
        hc  <- input$huge_criterion %||% "stars"
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default   = \"huge\",")
        L("  huge      = \"%s\",", hm)
        L("  criterion = \"%s\"%s", hc,
          if (hc == "ebic")
            sprintf(",\n  ebic.gamma = %s)", input$huge_ebic_gamma %||% 0.5)
          else if (hc == "stars")
            sprintf(",\n  stars.thresh = %s)", input$huge_stars_thresh %||% 0.1)
          else
            ")")
      } else if (def == "IsingFit") {
        gamma <- input$ising_gamma %||% 0.25
        rule  <- if (isTRUE(input$ising_and)) "AND" else "OR"
        split_m <- input$ising_split %||% "none"
        sp_val  <- switch(split_m, median = '"median"', mean = '"mean"', "FALSE")
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default = \"IsingFit\",")
        L("  gamma   = %s,",  gamma)
        L("  rule    = \"%s\",", rule)
        L("  split   = %s)", sp_val)
      } else if (def == "mgm") {
        ord  <- input$order        %||% 2
        rule <- input$mgm_rule     %||% "AND"
        crit <- input$mgm_criterion %||% "EBIC"
        tun  <- input$tuning_mgm   %||% 0.25
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default   = \"mgm\",")
        L("  order     = %s,", ord)
        L("  rule      = \"%s\",", rule)
        L("  criterion = \"%s\",", crit)
        if (crit == "EBIC") L("  tuning    = %s,", tun)
        L("  missing   = \"%s\")", missing_m)
      } else if (def == "relimp") {
        norm   <- isTRUE(input$relimp_norm)
        sd_val <- input$structureDefault %||% "none"
        tun    <- input$tuning           %||% 0.5
        L("%s <- estimateNetwork(%s,", data_name, data_name)
        L("  default          = \"relimp\",")
        L("  normalized       = %s,", norm)
        L("  structureDefault = \"%s\",", sd_val)
        L("  tuning           = %s)", tun)
        # relimp min threshold post-proc
        min_thr <- as.numeric(input$relimpEstMin %||% 0)
        if (min_thr > 0) {
          Lblank()
          L("# Minimum edge threshold (relimpEstMin = %s)", min_thr)
          L("W <- as.matrix(%s$graph); storage.mode(W) <- \"double\"; W[is.na(W)] <- 0", data_name)
          L("W[abs(W) < %s] <- 0", min_thr)
          L("%s$graph <- W", data_name)
        }
      }
    }

    if (use_split) {
      grp_vals <- names(res$groups)
      for (g in grp_vals) {
        gn <- gsub("[^A-Za-z0-9_]", "_", g)
        L("# Group: %s", g)
        L("net_%s <- group_%s", gn, gn)
        est_call(paste0("net_", gn))
        Lblank()
      }
    } else {
      L("net <- df_analysis")
      est_call("net")
    }
    Lblank()

    # 6. Bootstrap significance post-processing
    if (!is.null(boot_res) && isTRUE(input$boot_sig_enable)) {
      boot_alpha <- input$boot_sig_alpha      %||% 0.05
      boot_corr  <- input$boot_sig_correction %||% "none"
      boot_n     <- input$boot_sig_nboots     %||% 1000
      L("# ---- 6. Bootstrap edge significance ----")
      L("# Reproduce the bootstrap threshold applied in the app")
      if (use_split) {
        for (g in names(res$groups)) {
          gn <- gsub("[^A-Za-z0-9_]", "_", g)
          L("boot_%s <- bootnet::bootnet(net_%s, nBoots = %s,", gn, gn, boot_n)
          L("                            type = \"nonparametric\", statistics = \"edge\", nCores = 1)")
        }
      } else {
        L("boot_result <- bootnet::bootnet(net, nBoots = %s,", boot_n)
        L("                               type = \"nonparametric\", statistics = \"edge\", nCores = 1)")
      }
      Lblank()
      L("# Apply bootstrap threshold: correction = '%s', alpha = %s", boot_corr, boot_alpha)
      L("apply_boot_threshold <- function(net, boot_result, alpha = %s, correction = \"%s\") {",
        boot_alpha, boot_corr)
      L("  W     <- net$graph")
      L("  boots <- boot_result$boots")
      L("  pairs <- which(upper.tri(W), arr.ind = TRUE)")
      L("  pvals <- numeric(nrow(pairs))")
      L("  for (k in seq_len(nrow(pairs))) {")
      L("    i <- pairs[k,1]; j <- pairs[k,2]")
      L("    be <- sapply(boots, function(b) if(!is.null(b$graph)) b$graph[i,j] else NA_real_)")
      L("    be <- be[!is.na(be)]")
      L("    if (!length(be)) { pvals[k] <- 1; next }")
      L("    F0 <- mean(be <= 0)")
      L("    pvals[k] <- 2 * min(F0, 1 - F0)")
      L("  }")
      L("  pa <- if (correction != \"none\") p.adjust(pvals, method = correction) else pvals")
      L("  for (k in seq_len(nrow(pairs))) {")
      L("    if (pa[k] >= alpha) { i <- pairs[k,1]; j <- pairs[k,2]; W[i,j] <- W[j,i] <- 0 }")
      L("  }")
      L("  net$graph <- W; net")
      L("}")
      Lblank()
      if (use_split) {
        for (g in names(res$groups)) {
          gn <- gsub("[^A-Za-z0-9_]", "_", g)
          L("net_%s <- apply_boot_threshold(net_%s, boot_%s)", gn, gn, gn)
        }
      } else {
        L("net <- apply_boot_threshold(net, boot_result)")
      }
      Lblank()
    }

    # 7. Inspect
    L("# ---- 7. Inspect results ----")
    is_directed <- def %in% c("relimp", "mgm", "IsingFit")
    if (use_split) {
      for (g in names(res$groups)) {
        gn <- gsub("[^A-Za-z0-9_]", "_", g)
        L("print(net_%s$graph)", gn)
        L("qgraph::qgraph(net_%s$graph, layout = \"spring\", directed = %s,",
          gn, is_directed)
        L("               title = \"Group: %s\")", g)
      }
    } else {
      L("print(net$graph)")
      L("qgraph::qgraph(net$graph, layout = \"spring\", directed = %s)", is_directed)
    }

    paste(lines, collapse = "\n")
  }

  estimate_single_network <- function(df, input) {
    def <- input$estimator

    # ---- Fallback debug: log all effective parameters to console ----
    .fb <- list(
      estimator    = resolve(input$estimator,      "EBICglasso"),
      missing      = resolve(input$missing_policy, "listwise"),
      center_scale = resolve(input$center_scale,   FALSE)
    )
    if (identical(def, "EBICglasso")) {
      .fb$corMethod   <- resolve(input$corMethod_ebic,  "cor")
      .fb$tuning      <- resolve(input$tuning,           0.5)
      .fb$refit       <- resolve(input$refit,            FALSE)
      .fb$sampleSize  <- resolve(input$ebic_samplesize,  "maximum")
    } else if (identical(def, "pcor")) {
      .fb$corMethod   <- resolve(input$corMethod, "cor_auto")
    } else if (identical(def, "cor")) {
      .fb$corMethod        <- resolve(input$corMethod,        "cor_auto")
      .fb$threshold_method <- resolve(input$threshold_method, "none")
      .fb$threshold_value  <- resolve(input$threshold_value,  0.1)
    } else if (identical(def, "huge")) {
      .fb$huge_method    <- resolve(input$huge_method,       "glasso")
      .fb$huge_criterion <- resolve(input$huge_criterion,    "stars")
      .fb$stars_thresh   <- resolve(input$huge_stars_thresh, 0.1)
    } else if (identical(def, "IsingFit")) {
      .fb$gamma  <- resolve(input$ising_gamma, 0.25)
      .fb$rule   <- resolve(input$ising_and,   TRUE)
      .fb$split  <- resolve(input$ising_split, "none")
    } else if (identical(def, "mgm")) {
      .fb$order     <- resolve(input$order,         2)
      .fb$rule      <- resolve(input$mgm_rule,      "AND")
      .fb$criterion <- resolve(input$mgm_criterion, "EBIC")
      .fb$tuning    <- resolve(input$tuning_mgm,    0.25)
    } else if (identical(def, "relimp")) {
      .fb$normalized       <- resolve(input$relimp_norm,             TRUE)
      .fb$structureDefault <- resolve(input$structureDefault,        "none")
      .fb$tuning           <- resolve(input$tuning,                  0.5)
      .fb$smart_preprocess <- resolve(input$relimp_smart_preprocess, TRUE)
      .fb$relimpEstMin     <- resolve(input$relimpEstMin,            0)
    }
    debug_params(.fb)

    # Validation: Check if non-continuous variables are included for EBICglasso/pcor/cor/huge
    if (def %in% c("EBICglasso","pcor","cor","huge")){
      non_continuous <- setdiff(colnames(df), cont_vars())
      if (length(non_continuous) > 0) {
        showNotification(
          paste0("⚠️ ", def, " requires ALL variables to be Continuous (Gaussian).\n",
                 "Non-continuous variables found: ", paste(non_continuous, collapse=", "), "\n",
                 "Please move these to 'Available Variables' or remove them from the analysis."),
          type = "error",
          duration = 10
        )
        return(NULL)
      }

      df2 <- df[, sapply(df, is_numeric_like), drop = FALSE]
      if (ncol(df2) < 2) {
        showNotification("Need at least 2 continuous variables", type = "error")
        return(NULL)
      }

      # Get correlation method based on estimator
      cor_method <- if (def == "EBICglasso") {
        input$corMethod_ebic %||% "cor"
      } else {
        input$corMethod %||% "cor_auto"
      }

      # Get missing data handling method from Variables panel
      missing_method <- input$missing_policy %||% "listwise"

      # Debug info about sample sizes
      n_total <- nrow(df2)
      n_complete <- sum(complete.cases(df2))
      n_with_na <- n_total - n_complete

      message(paste("\n=== ESTIMATION DEBUG ==="))
      message(paste("Estimator:", def))
      message(paste("Correlation method:", cor_method))
      message(paste("Missing policy:", missing_method))
      message(paste("Total rows:", n_total))
      message(paste("Complete cases:", n_complete))
      message(paste("Rows with NA:", n_with_na))
      message(paste("Gamma (tuning):", input$tuning))

      # Calculate pairwise sample sizes if using pairwise
      if (missing_method == "pairwise" && n_with_na > 0) {
        # Calculate pairwise complete observations for each pair of variables
        pairwise_ns <- matrix(0, ncol(df2), ncol(df2))
        for (i in 1:ncol(df2)) {
          for (j in i:ncol(df2)) {
            pairwise_ns[i,j] <- pairwise_ns[j,i] <- sum(complete.cases(df2[, c(i,j)]))
          }
        }
        diag(pairwise_ns) <- n_total

        message(paste("Pairwise N - Min:", min(pairwise_ns[lower.tri(pairwise_ns)])))
        message(paste("Pairwise N - Max:", max(pairwise_ns[lower.tri(pairwise_ns)])))
        message(paste("Pairwise N - Mean:", round(mean(pairwise_ns[lower.tri(pairwise_ns)]), 1)))
        message(paste("JASP 'Maximum' equivalent: Use N =", n_total))
      }
      message("========================\n")

      # Additional debug: show data summary and correlations
      message("\n=== DATA SUMMARY BEFORE ESTIMATION ===")
      message(paste("Variables:", paste(colnames(df2), collapse = ", ")))
      message(paste("Dimensions:", nrow(df2), "rows x", ncol(df2), "cols"))

      # Show summary statistics for each variable
      for (var in colnames(df2)) {
        var_data <- df2[[var]]
        message(sprintf("  %s: mean=%.3f, sd=%.3f, range=[%.3f, %.3f], NA=%d",
                       var,
                       mean(var_data, na.rm = TRUE),
                       sd(var_data, na.rm = TRUE),
                       min(var_data, na.rm = TRUE),
                       max(var_data, na.rm = TRUE),
                       sum(is.na(var_data))))
      }

      # Compute and show correlation matrix
      if (missing_method == "pairwise") {
        test_cor <- cor(df2, use = "pairwise.complete.obs")
      } else {
        test_cor <- cor(df2, use = "complete.obs")
      }
      message("\n=== PEARSON CORRELATION MATRIX ===")
      print(round(test_cor, 3))
      message(paste("\nCorrelation range:", round(min(test_cor[lower.tri(test_cor)]), 3),
                   "to", round(max(test_cor[lower.tri(test_cor)]), 3)))
      message(paste("Mean |r|:", round(mean(abs(test_cor[lower.tri(test_cor)])), 3)))
      message("=====================================\n")

      # Store original variable names for restoration after bootnet
      # (bootnet may abbreviate names in the network matrix)
      original_var_names <- colnames(df2)

      if (def == "EBICglasso"){
        # Get sample size method for pairwise deletion
        samplesize_method <- input$ebic_samplesize %||% "pairwise_average"

        # Determine bootnet's sampleSize parameter (must be character string)
        # bootnet accepts: "maximum", "minimum", "pairwise_average", "pairwise_maximum", "pairwise_minimum"
        if (missing_method == "pairwise" && n_with_na > 0) {
          bootnet_samplesize <- switch(samplesize_method,
            "maximum" = "maximum",                      # Use total N
            "minimum" = "minimum",                      # Use complete cases N
            "pairwise_average" = "pairwise_average",    # Use average pairwise N (bootnet default)
            "pairwise_average"  # default fallback
          )
          message(paste("Using sampleSize method:", bootnet_samplesize, "(UI selection:", samplesize_method, ")"))
          message(paste("This will use N =", switch(bootnet_samplesize,
                                                      "maximum" = n_total,
                                                      "minimum" = n_complete,
                                                      paste0("~", round(mean(pairwise_ns[lower.tri(pairwise_ns)]), 1))),
                        "for EBIC calculation"))
        } else {
          bootnet_samplesize <- "maximum"  # Default for non-pairwise
        }

        # Call bootnet with JASP-compatible parameters
        # JASP defaults: weighted=TRUE, signed=TRUE
        # NOTE: JASP does NOT pass lambda.min.ratio - uses bootnet's default (0.01)
        net <- bootnet::estimateNetwork(df2, default = "EBICglasso", corMethod = cor_method,
                                 tuning = input$tuning, refit = input$refit,
                                 missing = missing_method,
                                 sampleSize = bootnet_samplesize,
                                 weighted = TRUE,
                                 signed = TRUE)

        # Restore original variable names (bootnet may abbreviate them)
        if (!is.null(net$graph) && !is.null(original_var_names)) {
          if (ncol(net$graph) == length(original_var_names)) {
            colnames(net$graph) <- original_var_names
            rownames(net$graph) <- original_var_names
          }
        }
        net

      } else if (def == "pcor") {
        net <- bootnet::estimateNetwork(df2, default = "pcor", corMethod = cor_method,
                                missing = missing_method)

        # Restore original variable names (bootnet may abbreviate them)
        if (!is.null(net$graph) && !is.null(original_var_names)) {
          if (ncol(net$graph) == length(original_var_names)) {
            colnames(net$graph) <- original_var_names
            rownames(net$graph) <- original_var_names
          }
        }
        net

      } else if (def == "cor") {
        # Correlation network with optional thresholding
        net <- bootnet::estimateNetwork(df2, default = "cor", corMethod = cor_method,
                                        missing = missing_method)

        # Apply threshold if specified
        threshold_method <- input$threshold_method %||% "none"
        if (threshold_method != "none") {
          net <- apply_threshold(net, df2, threshold_method, input$threshold_value)
        }

        # Restore original variable names (bootnet may abbreviate them)
        if (!is.null(net$graph) && !is.null(original_var_names)) {
          if (ncol(net$graph) == length(original_var_names)) {
            colnames(net$graph) <- original_var_names
            rownames(net$graph) <- original_var_names
          }
        }
        net

      } else if (def == "huge") {
        # Huge: High-dimensional undirected graph estimation
        message("\n=== Huge Estimation ===")
        huge_method <- input$huge_method %||% "glasso"
        huge_criterion <- input$huge_criterion %||% "stars"

        message(paste("Method:", huge_method))
        message(paste("Criterion:", huge_criterion))

        net <- bootnet::estimateNetwork(
          df2,
          default = "huge",
          huge = huge_method,
          criterion = huge_criterion,
          ebic.gamma = if (huge_criterion == "ebic") input$huge_ebic_gamma %||% 0.5 else NULL,
          stars.thresh = if (huge_criterion == "stars") input$huge_stars_thresh %||% 0.1 else NULL,
          stars.subsample.ratio = if (huge_criterion == "stars") {
            n <- nrow(df2)
            ratio_input <- input$huge_stars_subsample %||% 10
            # Convert number of subsamples to ratio
            # Default is 10*sqrt(n)/n when n>144, 0.8 when n<=144
            if (n > 144) {
              min(10 * sqrt(n) / n, 0.9)
            } else {
              0.8
            }
          } else NULL,
          rep.num = if (huge_criterion == "ric") input$huge_ric_num %||% 10 else NULL
        )

        # Restore original variable names (bootnet may abbreviate them)
        if (!is.null(net$graph) && !is.null(original_var_names)) {
          if (ncol(net$graph) == length(original_var_names)) {
            colnames(net$graph) <- original_var_names
            rownames(net$graph) <- original_var_names
          }
        }
        net
      }
    } else if (def == "IsingFit") {
      message("\n=== IsingFit Estimation ===")
      ising_gamma <- input$ising_gamma %||% 0.25
      ising_and <- input$ising_and %||% TRUE

      message(paste("EBIC gamma:", ising_gamma))
      message(paste("AND-rule:", ising_and))

      # IsingFit requires binary (0/1) data - NOT scaled/centered
      # The df parameter is prepped_data which might be scaled - we need raw data

      # Get variable names from the passed df
      var_names <- colnames(df)

      # Get raw unscaled data
      full_raw <- dat()
      if (is.null(full_raw)) {
        showNotification("Data not available", type = "error")
        return(NULL)
      }

      # Use rownames to get the exact same rows from raw data
      # This ensures split analyses get the correct group's data
      row_indices <- as.integer(rownames(df))
      df_binary <- full_raw[row_indices, var_names, drop = FALSE]

      message(paste("Using rows:", min(row_indices), "to", max(row_indices), "( n =", length(row_indices), ")"))

      # Convert to numeric binary (0/1) - ensure it's truly binary
      for (col_name in colnames(df_binary)) {
        col <- df_binary[[col_name]]
        if (is.factor(col)) {
          # Convert factor to numeric (levels become 0, 1, 2, ...)
          df_binary[[col_name]] <- as.numeric(col) - 1
        } else if (is.logical(col)) {
          # Convert logical to 0/1
          df_binary[[col_name]] <- as.numeric(col)
        } else if (is.numeric(col)) {
          # Ensure it's numeric
          df_binary[[col_name]] <- as.numeric(col)
        }
      }

      message(paste("Variables:", ncol(df_binary)))
      message(paste("Observations:", nrow(df_binary)))

      # Check data range to verify it's binary
      data_ranges <- sapply(df_binary, function(x) paste(min(x, na.rm=TRUE), "-", max(x, na.rm=TRUE)))
      message("Data ranges (should be 0-1 for binary):")
      message(paste(capture.output(print(head(data_ranges, 10))), collapse="\n"))

      # Get binarization option
      split_method <- input$ising_split %||% "none"

      # Convert split method to bootnet parameter
      split_param <- if (split_method == "none") {
        FALSE  # Data is already binary
      } else if (split_method == "median") {
        "median"
      } else if (split_method == "mean") {
        "mean"
      } else {
        FALSE
      }

      message(paste("Binarization method:", split_method))

      # bootnet's IsingFit default accepts: tuning (gamma), rule (AND/OR), and split
      bootnet::estimateNetwork(
        df_binary,
        default = "IsingFit",
        tuning = ising_gamma,
        rule = if(ising_and) "AND" else "OR",
        split = split_param
      )
    } else if (def == "mgm"){
      vars <- colnames(df)

      # Determine variable types: g (Gaussian), c (categorical), p (Poisson/count)
      typ <- ifelse(vars %in% cont_vars(), "g",
                    ifelse(vars %in% count_vars(), "p",
                           ifelse(vars %in% disc_vars(), "c",
                                  ifelse(sapply(df, is_numeric_like),"g","c"))))

      # Prepare data for MGM: convert categorical variables to integer codes
      df_mgm <- df
      for (i in seq_along(vars)) {
        if (typ[i] == "c") {
          # Convert categorical to integer codes (1, 2, 3, ...)
          df_mgm[[vars[i]]] <- as.integer(factor(df[[vars[i]]]))
        }
      }

      # Determine levels for each variable
      lvl <- vapply(seq_along(vars), function(i){
        if (typ[i] == "g" || typ[i] == "p") {
          1L  # Gaussian and Poisson have level 1
        } else {
          length(levels(factor(df[[i]])))  # Categorical has number of categories
        }
      }, integer(1))

      # Show summary of variable types for MGM
      message("\n=== MGM Variable Types ===")
      type_summary <- data.frame(
        Variable = vars,
        Type = typ,
        Type_Name = ifelse(typ == "g", "Gaussian",
                          ifelse(typ == "p", "Poisson",
                                 ifelse(typ == "c", "Categorical", "Unknown"))),
        Level = lvl
      )
      message(paste(capture.output(print(type_summary, row.names = FALSE)), collapse = "\n"))
      message("=== END MGM Variable Types ===\n")

      # Prepare criterion parameter based on selection
      criterion <- input$mgm_criterion %||% "EBIC"

      message("\n=== MGM Estimation Parameters ===")
      message(paste("Order:", input$order))
      message(paste("Rule:", input$mgm_rule))
      message(paste("Criterion:", criterion))

      if (criterion == "CV") {
        # Cross-validation with specified folds
        n_folds <- input$mgm_cv_folds %||% 10
        message(paste("CV folds:", n_folds))
        message("=== END MGM Parameters ===\n")
        bootnet::estimateNetwork(df_mgm, default = "mgm", type = typ, level = lvl,
                                 order = input$order, rule = input$mgm_rule,
                                 criterion = "CV", nFolds = n_folds)
      } else {
        # EBIC with gamma parameter (default)
        gamma_val <- input$tuning_mgm %||% 0.25
        message(paste("EBIC gamma:", gamma_val))
        message("=== END MGM Parameters ===\n")
        bootnet::estimateNetwork(df_mgm, default = "mgm", type = typ, level = lvl,
                                 order = input$order, rule = input$mgm_rule,
                                 criterion = "EBIC", tuning = gamma_val)
      }
    } else if (def == "relimp"){
      message("\n=== Relative Importance Estimation ===")
      message(paste("Input data dimensions:", nrow(df), "rows x", ncol(df), "columns"))
      message(paste("Row range:", paste(range(as.integer(rownames(df))), collapse = " to ")))

      X <- tryCatch(
        make_relimp_matrix(
          df,
          smart_preprocess = isTRUE(input$relimp_smart_preprocess),
          dummy = isTRUE(input$relimp_dummy),
          drop_nonfinite = isTRUE(input$relimp_drop_nonfinite),
          drop_zerovar = isTRUE(input$relimp_drop_zerovar)
        ),
        error = function(e) {
          showNotification(paste("Relimp preprocessing:", e$message), type = "error", duration = 6)
          NULL
        }
      )
      if (is.null(X)) return(NULL)

      message(paste("After preprocessing:", nrow(X), "rows x", ncol(X), "columns"))
      
      # Store original names before estimation
      original_names <- attr(X, "original_names")
      if (is.null(original_names)) original_names <- colnames(X)
      
      # Estimate network
      net <- tryCatch({
        # For IsingFit structure baseline, pass additional parameters
        if (input$structureDefault == "IsingFit") {
          # Get binarization option from IsingFit settings
          split_method <- input$ising_split %||% "none"
          split_param <- if (split_method == "none") {
            FALSE
          } else if (split_method == "median") {
            "median"
          } else if (split_method == "mean") {
            "mean"
          } else {
            FALSE
          }

          message("Using IsingFit structure baseline")
          message(paste("Binarization method:", split_method))

          bootnet::estimateNetwork(
            X, default = "relimp",
            normalized = isTRUE(input$relimp_norm),
            structureDefault = "IsingFit",
            tuning = input$tuning,
            # IsingFit-specific parameters
            rule = if(isTRUE(input$ising_and)) "AND" else "OR",
            split = split_param
          )
        } else {
          bootnet::estimateNetwork(
            X, default = "relimp",
            normalized = isTRUE(input$relimp_norm),
            structureDefault = input$structureDefault,
            tuning = input$tuning
          )
        }
      }, error = function(e) {
        showNotification(
          paste0("Relimp estimation failed: ", e$message,
                 "\nTip: Try using a different structure baseline or check your variable names."),
          type = "error",
          duration = 8
        )
        NULL
      })
      
      if (is.null(net)) return(NULL)
      
      # Restore original variable names in the network object
      if (!is.null(net$graph) && !is.null(original_names)) {
        if (ncol(net$graph) == length(original_names)) {
          colnames(net$graph) <- original_names
          rownames(net$graph) <- original_names
        }
      }
      
      return(net)
    } else stop("Unknown estimator")
  }
  
  est_result_thr <- reactive({
    res <- est_result()
    if (is.null(res)) return(NULL)

    boot_res      <- boot_sig_result()
    use_boot      <- isTRUE(input$boot_sig_enable) && !is.null(boot_res)
    alpha_val     <- input$boot_sig_alpha      %||% 0.05
    correction_val <- input$boot_sig_correction %||% "none"

    if (res$split) {
      for (g in names(res$groups)) {
        if (identical(input$estimator, "relimp")) {
          res$groups[[g]] <- apply_relimp_threshold(res$groups[[g]],
                                                    thr = as.numeric(input$relimpEstMin %||% 0))
        }
        if (use_boot) {
          res$groups[[g]] <- tryCatch(
            apply_boot_threshold(res$groups[[g]], boot_res, alpha_val, correction_val),
            error = function(e) res$groups[[g]]
          )
        }
      }
      res
    } else {
      if (identical(input$estimator, "relimp")) {
        res$network <- apply_relimp_threshold(res$network,
                                              thr = as.numeric(input$relimpEstMin %||% 0))
      }
      if (use_boot) {
        res$network <- tryCatch(
          apply_boot_threshold(res$network, boot_res, alpha_val, correction_val),
          error = function(e) res$network
        )
      }
      res
    }
  })
  
  predictability_vec <- reactive({
    df <- prepped_data()
    if (is.null(df)) return(NULL)
    split <- split_var()
    full_df <- dat()
    if (!is.null(split) && split != "" && !is.null(full_df) && split %in% colnames(full_df)) {
      split_values <- unique(na.omit(full_df[[split]]))
      pred_list <- list()
      for (group in split_values) {
        group_idx <- which(full_df[[split]] == group)
        group_df <- df[group_idx, , drop = FALSE]
        if (nrow(group_df) < 10) {
          pred_list[[as.character(group)]] <- NULL
          next
        }
        group_df_num <- group_df[, sapply(group_df, is_numeric_like), drop = FALSE]
        if (ncol(group_df_num) < 2) {
          pred_list[[as.character(group)]] <- NULL
          next
        }
        r2 <- predictability_R2(group_df_num)
        out <- rep(0, ncol(group_df))
        if (ncol(group_df_num) > 0 && length(r2) > 0) {
          idx <- match(colnames(group_df_num), colnames(group_df))
          out[idx[!is.na(idx)]] <- r2
        }
        names(out) <- colnames(group_df)
        pred_list[[as.character(group)]] <- out
      }
      return(pred_list)
    } else {
      df_num <- df[, sapply(df, is_numeric_like), drop = FALSE]
      r2 <- predictability_R2(df_num)
      out <- rep(0, ncol(df))
      if (ncol(df_num) > 0 && length(r2) > 0) {
        idx <- match(colnames(df_num), colnames(df))
        out[idx[!is.na(idx)]] <- r2
      }
      names(out) <- colnames(df)
      return(out)
    }
  })
  
  rescale_edge_weights <- function(W, minS, maxS){
    A <- abs(W)
    if (!is.finite(minS)) minS <- 0
    if (!is.finite(maxS) || maxS <= minS) maxS <- max(A, na.rm = TRUE)
    if (!is.finite(maxS) || maxS == 0) return(W*0)
    S <- (A - minS) / (maxS - minS); S[S < 0] <- 0; S[S > 1] <- 1
    S * sign(W)
  }
  safe_get_W <- function(net){
    W <- tryCatch(bootnet::getWmat(net), error = function(e) NULL)
    if (!is.null(W)) return(W)
    if (!is.null(net$graph) && is.matrix(net$graph)) return(net$graph)
    NULL
  }
  
  output$networkPlot <- renderPlot({
    res <- est_result_thr()
    shiny::validate(shiny::need(!is.null(res), "Estimate a network first."))
    df <- prepped_data()
    shiny::validate(shiny::need(!is.null(df), "No usable data."))
    
    # Sort groups alphabetically/numerically for consistent ordering
    # Do this FIRST before any computations so all tables/plots use same order
    if (res$split && !is.null(res$groups) && length(res$groups) > 1) {
      # Sort group names (0, 1, 2... or A, B, C...)
      sorted_groups <- sort(names(res$groups))
      # Reorder res$groups
      res$groups <- res$groups[sorted_groups]
      message(paste("Groups ordered alphabetically/numerically:", paste(sorted_groups, collapse=", ")))
    }
    
    # Compute node means from RAW data (before standardization)
    if (isTRUE(input$show_node_means)) {
      raw_df <- dat()  # Get raw data BEFORE preprocessing
      
      if (!is.null(raw_df) && !is.null(df)) {
        # Get the variables that are actually in the network
        network_vars <- colnames(df)
        
        # Filter raw data to only include network variables
        if (all(network_vars %in% colnames(raw_df))) {
          raw_network_data <- raw_df[, network_vars, drop = FALSE]

          # Convert all columns to numeric (for categorical/factor variables)
          # This matches the MGM preprocessing where categorical vars become integers
          for (col in network_vars) {
            if (!is.numeric(raw_network_data[[col]])) {
              # Convert factors/characters to integer codes
              raw_network_data[[col]] <- as.numeric(as.integer(factor(raw_network_data[[col]])))
            }
          }

          # For split analysis, compute per-group means and ANOVA
          if (res$split && !is.null(res$split_var) && res$split_var %in% colnames(raw_df)) {
            split_var_data <- raw_df[[res$split_var]]
            # Use the already-sorted group order from res$groups
            groups <- names(res$groups)
            
            # Initialize data frame for results
            means_df <- data.frame(Variable = network_vars, stringsAsFactors = FALSE)
            
            # Compute means for each group
            for (g in groups) {
              group_idx <- which(split_var_data == g)
              group_data <- raw_network_data[group_idx, , drop = FALSE]
              group_means <- colMeans(group_data, na.rm = TRUE)
              means_df[[paste0("Group_", g, "_Mean")]] <- round(group_means, 3)
              group_sds <- apply(group_data, 2, sd, na.rm = TRUE)
              means_df[[paste0("Group_", g, "_SD")]] <- round(group_sds, 3)
            }
            
            # Perform ANOVA for each variable
            f_values <- numeric(length(network_vars))
            p_values <- numeric(length(network_vars))
            
            for (i in seq_along(network_vars)) {
              var_name <- network_vars[i]
              if (var_name %in% colnames(raw_df)) {
                # Create data frame for ANOVA
                anova_data <- data.frame(
                  value = raw_df[[var_name]],
                  group = factor(split_var_data)
                )
                anova_data <- anova_data[complete.cases(anova_data), ]
                
                if (nrow(anova_data) > 0 && length(unique(anova_data$group)) > 1) {
                  # Perform one-way ANOVA (between-subjects)
                  anova_result <- tryCatch({
                    aov_fit <- aov(value ~ group, data = anova_data)
                    summary(aov_fit)
                  }, error = function(e) NULL)
                  
                  if (!is.null(anova_result)) {
                    f_values[i] <- anova_result[[1]]$`F value`[1]
                    p_values[i] <- anova_result[[1]]$`Pr(>F)`[1]
                  } else {
                    f_values[i] <- NA
                    p_values[i] <- NA
                  }
                } else {
                  f_values[i] <- NA
                  p_values[i] <- NA
                }
              }
            }
            
            means_df$F_value <- round(f_values, 3)
            means_df$p_value <- round(p_values, 4)
            
            values$node_means <- means_df
            
          } else {
            # No split - compute overall means
            overall_means <- colMeans(raw_network_data, na.rm = TRUE)
            overall_sds <- apply(raw_network_data, 2, sd, na.rm = TRUE)
            
            values$node_means <- data.frame(
              Variable = network_vars,
              Mean = round(overall_means, 3),
              SD = round(overall_sds, 3),
              stringsAsFactors = FALSE
            )
          }
        } else {
          values$node_means <- NULL
        }
      } else {
        values$node_means <- NULL
      }
    } else {
      values$node_means <- NULL
    }
    
    # Compute split predictability if applicable - USE SAME METHOD AS PIES!
    if (res$split && isTRUE(input$predictability_pies)) {
      # Use the SAME predictability_vec() that's used for pies
      all_pred <- predictability_vec()

      if (!is.null(all_pred) && is.list(all_pred)) {
        # Use the already-sorted group order from res$groups
        group_names_sorted <- names(res$groups)

        # Convert list to data frame - use network variables, not all df variables
        # Get network variables from the first group
        all_vars <- NULL
        for (group_name in group_names_sorted) {
          net <- res$groups[[group_name]]
          if (!is.null(net)) {
            W <- safe_get_W(net)
            if (!is.null(W) && !is.null(colnames(W))) {
              all_vars <- colnames(W)
              break
            }
          }
        }
        # Fallback to df columns if network vars not found
        if (is.null(all_vars)) all_vars <- colnames(df)
        pred_df <- data.frame(Variable = all_vars, stringsAsFactors = FALSE)
        
        for (group_name in group_names_sorted) {
          pred_values <- all_pred[[group_name]]
          if (!is.null(pred_values) && length(pred_values) > 0) {
            # Match to all_vars
            matched_values <- numeric(length(all_vars))
            for (i in seq_along(all_vars)) {
              if (all_vars[i] %in% names(pred_values)) {
                matched_values[i] <- pred_values[all_vars[i]]
              } else {
                matched_values[i] <- NA
              }
            }
            pred_df[[paste0("Group_", group_name)]] <- round(matched_values, 3)
          }
        }
        
        values$split_predictability <- pred_df
      } else {
        values$split_predictability <- NULL
      }
    } else {
      values$split_predictability <- NULL
    }
    
    if (res$split) {
      n_groups <- length(res$groups)
      if (n_groups == 0) {
        plot.new()
        text(0.5, 0.5, "No valid groups to plot", cex = 2)
        return()
      }
      
      # Groups already sorted at beginning of render
      ncol <- min(3, n_groups)
      nrow <- ceiling(n_groups / ncol)
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      par(mfrow = c(nrow, ncol), mar = c(2, 2, 3, 2))
      shared_layout <- NULL
      shared_communities <- NULL
      if (input$layout %in% c("spring", "ega", "pca")) {
        for (group_name in names(res$groups)) {
          net <- res$groups[[group_name]]
          if (!is.null(net)) {
            W <- safe_get_W(net)
            if (!is.null(W)) {
              W <- as.matrix(W)
              storage.mode(W) <- "double"
              W[is.na(W)] <- 0

              # Compute appropriate layout
              if (input$layout == "ega") {
                message(paste("Computing EGA layout for group:", group_name))
                # Get data for this group
                full_df <- dat()
                split <- split_var()
                group_data <- NULL
                if (!is.null(split) && split %in% colnames(full_df)) {
                  group_idx <- which(full_df[[split]] == group_name)
                  # Use full_df filtered, not df
                  group_data <- full_df[group_idx, colnames(W), drop = FALSE]
                } else {
                  group_data <- full_df[, colnames(W), drop = FALSE]
                }
                # Return communities for coloring
                ega_result <- compute_ega_layout(W, data = group_data, return_communities = TRUE,
                                                 separation = input$ega_separation %||% 3)
                if (is.list(ega_result)) {
                  shared_layout <- ega_result$layout
                  comm_ega  <- ega_result$communities
                  n_comm    <- ega_result$n_communities %||% length(unique(comm_ega))
                  # k-means fallback when EGA finds ≤1 community
                  if (is.null(comm_ega) || n_comm <= 1) {
                    lm_mat <- if (is.matrix(shared_layout)) shared_layout else NULL
                    if (!is.null(lm_mat)) {
                      k_eg  <- if (nrow(W) >= 6) 3L else 2L
                      km_eg <- tryCatch(kmeans(lm_mat, centers = min(k_eg, nrow(lm_mat)-1), nstart=20),
                                        error = function(e) NULL)
                      if (!is.null(km_eg)) {
                        comm_ega <- km_eg$cluster
                        n_comm   <- length(unique(km_eg$cluster))
                      }
                    }
                  }
                  shared_communities <- list(communities = comm_ega, n_communities = n_comm)
                } else {
                  shared_layout <- ega_result
                }
              } else if (input$layout == "pca") {
                # Get data for first group
                full_df <- dat()
                split <- split_var()
                if (!is.null(split) && split %in% colnames(full_df)) {
                  group_idx <- which(full_df[[split]] == group_name)
                  group_data <- full_df[group_idx, colnames(W), drop = FALSE]
                  shared_layout <- compute_pca_layout(W, data = group_data)
                } else {
                  shared_layout <- compute_pca_layout(W, data = full_df[, colnames(W), drop = FALSE])
                }
                # k-means on PCA layout coordinates for community coloring
                if (is.matrix(shared_layout)) {
                  k_pca <- if (nrow(W) >= 6) 3L else 2L
                  km_pca <- tryCatch(
                    kmeans(shared_layout, centers = min(k_pca, nrow(shared_layout) - 1), nstart = 20),
                    error = function(e) NULL)
                  if (!is.null(km_pca)) {
                    shared_communities <- list(
                      communities  = km_pca$cluster,
                      n_communities = length(unique(km_pca$cluster))
                    )
                  }
                }
              } else {
                # Spring layout
                shared_layout <- qgraph::qgraph(W, DoNotPlot = TRUE, layout = "spring")$layout
              }
              if (is.matrix(shared_layout)) main_net_layout(shared_layout)
              break
            }
          }
        }
      }

      # Compute means only for the NETWORK variables (not all variables from data)
      # This is critical for global scaling to work correctly
      # Get network variables from the first group's network matrix
      network_vars <- NULL
      for (group_name in names(res$groups)) {
        net <- res$groups[[group_name]]
        if (!is.null(net)) {
          W <- safe_get_W(net)
          if (!is.null(W) && !is.null(colnames(W))) {
            network_vars <- colnames(W)
            break
          }
        }
      }
      # Fallback to df columns if network_vars not found
      if (is.null(network_vars)) network_vars <- colnames(df)

      all_means_for_global <- NULL

      if (res$split) {
        full_df <- dat()
        split <- split_var()
        if (!is.null(split) && split %in% colnames(full_df)) {
          # Compute means for each group, but only for network variables
          all_means_for_global <- list()
          for (group_name in names(res$groups)) {
            group_idx <- which(full_df[[split]] == group_name)
            group_means <- colMeans(full_df[group_idx, network_vars, drop = FALSE], na.rm = TRUE)
            all_means_for_global[[as.character(group_name)]] <- group_means
          }
        }
      }
      
      all_means <- means_for_nodes()
      all_pred <- predictability_vec()
      
      # Compute global scaling if requested
      global_scaled_means <- NULL
      if (isTRUE(input$nodeSizeByMean) && isTRUE(input$useGlobalScaling) && is.list(all_means_for_global)) {
        message("=== GLOBAL SCALING DEBUG ===")
        message(paste("Number of groups:", length(all_means_for_global)))
        message(paste("Group names:", paste(names(all_means_for_global), collapse=", ")))
        
        # Step 1: Collect ALL means from ALL groups with their identifiers
        all_varnames <- unique(unlist(lapply(all_means_for_global, names)))
        all_groups <- names(all_means_for_global)
        message(paste("Number of NETWORK variables:", length(all_varnames)))
        message(paste("Variables:", paste(head(all_varnames, 5), collapse=", ")))
        
        # Step 2: Create a data frame with all combinations
        means_df <- data.frame()
        for (group_name in all_groups) {
          group_vals <- all_means_for_global[[group_name]]
          message(paste("Group", group_name, "has", length(group_vals), "values"))
          for (var in all_varnames) {
            if (var %in% names(group_vals) && is.finite(group_vals[var])) {
              means_df <- rbind(means_df, data.frame(
                group = group_name,
                variable = var,
                mean_value = group_vals[var],
                stringsAsFactors = FALSE
              ))
            }
          }
        }
        
        message(paste("Total rows in means_df:", nrow(means_df)))
        if (nrow(means_df) > 0) {
          message("Sample rows:")
          message(paste(capture.output(print(head(means_df, 6))), collapse="\n"))
          
          # Step 3: Scale ALL values together in one operation
          use_reverse <- isTRUE(input$useReverseScaling)
          message(paste("Using reverse scaling:", use_reverse))
          size_range <- input$nodeSizeRange
          if (is.null(size_range) || length(size_range) != 2 || anyNA(size_range)) {
            size_range <- c(4, 10)  # default range
          }
          means_df$scaled_size <- scale_to_range(means_df$mean_value,
                                                 rng = size_range,
                                                 reverse = use_reverse)
          
          message("After scaling:")
          message(paste(capture.output(print(head(means_df, 6))), collapse="\n"))
          
          # Step 4: Create lookup structure for each group
          global_scaled_means <- list()
          for (group_name in all_groups) {
            group_subset <- means_df[means_df$group == group_name, ]
            scaled_vals <- setNames(group_subset$scaled_size, group_subset$variable)
            global_scaled_means[[group_name]] <- scaled_vals
            message(paste("Group", group_name, "lookup has", length(scaled_vals), "values"))
            message(paste("Sample:", paste(names(head(scaled_vals, 3)), "=", round(head(scaled_vals, 3), 2), collapse=", ")))
          }
          
          # Message
          global_min <- min(means_df$mean_value, na.rm = TRUE)
          global_max <- max(means_df$mean_value, na.rm = TRUE)
          if (use_reverse) {
            message(sprintf("Global REVERSE scaling: ALL values scaled together, range [%.3f, %.3f]", global_min, global_max))
          } else {
            message(sprintf("Global STANDARD scaling: ALL values scaled together, range [%.3f, %.3f]", global_min, global_max))
          }
          message("=== END GLOBAL SCALING DEBUG ===")
        }
      }
      
      for (group_name in names(res$groups)) {
        net <- res$groups[[group_name]]
        if (is.null(net)) {
          plot.new()
          text(0.5, 0.5, paste("Group", group_name, "\nNo data"), cex = 1.2)
          next
        }
        group_means <- if (is.list(all_means)) all_means[[group_name]] else all_means
        group_pred <- if (is.list(all_pred)) all_pred[[group_name]] else all_pred
        
        # Get pre-computed scaled sizes if using global scaling
        group_scaled_sizes <- NULL
        if (!is.null(global_scaled_means) && group_name %in% names(global_scaled_means)) {
          group_scaled_sizes <- global_scaled_means[[group_name]]
          message(paste("Retrieved scaled sizes for group", group_name, ":", length(group_scaled_sizes), "values"))
        } else {
          message(paste("No scaled sizes for group", group_name))
          message(paste("  global_scaled_means is NULL:", is.null(global_scaled_means)))
          if (!is.null(global_scaled_means)) {
            message(paste("  group_name in names:", group_name %in% names(global_scaled_means)))
          }
        }
        
        plot_single_network(net, df, input,
                            title = paste("Group:", group_name),
                            layout_override = shared_layout,
                            means_override = group_means,
                            pred_override = group_pred,
                            scaled_sizes_override = group_scaled_sizes,
                            communities_override = shared_communities)
      }
    } else {
      plot_single_network(res$network, df, input)
    }
  })
  
  plot_single_network <- function(net, df, input, title = NULL,
                                  layout_override = NULL,
                                  means_override = NULL,
                                  pred_override = NULL,
                                  scaled_sizes_override = NULL,
                                  communities_override = NULL) {
    # Debug at start of function
    message("\n=== PLOT_SINGLE_NETWORK START ===")
    message(paste("Title:", title %||% "(no title)"))
    message(paste("input$predictability_pies:", isTRUE(input$predictability_pies)))
    message(paste("pred_override is.null:", is.null(pred_override)))
    message("=== END PLOT START ===\n")

    # Validate inputs
    if (is.null(net) || is.null(df)) {
      plot.new()
      text(0.5, 0.5, "Network or data is NULL", cex = 1.5)
      return(invisible(NULL))
    }

    W_try <- safe_get_W(net)
    use_direct_qgraph <- !is.null(W_try) && is.numeric(W_try)

    # Use network matrix column names (not df column names) to avoid mismatch
    # For correlation networks, df may have more variables than the network matrix
    varnames <- if (!is.null(W_try) && !is.null(colnames(W_try))) {
      colnames(W_try)
    } else {
      colnames(df)
    }

    if (is.null(varnames) || length(varnames) == 0) {
      plot.new()
      text(0.5, 0.5, "No variable names found", cex = 1.5)
      return(invisible(NULL))
    }
    device_bg <- "white"
    par(bg = device_bg)
    
    pies_all <- if (isTRUE(input$predictability_pies)) {
      pred_vals <- if (!is.null(pred_override)) pred_override else predictability_vec()

      # Diagnostic output
      message("\n=== PREDICTABILITY DEBUG ===")
      message(paste("predictability_pies enabled:", isTRUE(input$predictability_pies)))
      message(paste("pred_vals is.null:", is.null(pred_vals)))
      if (!is.null(pred_vals)) {
        message(paste("pred_vals class:", paste(class(pred_vals), collapse=", ")))
        message(paste("pred_vals length:", length(pred_vals)))
        if (!is.list(pred_vals) && length(pred_vals) > 0) {
          n_nonzero <- sum(pred_vals > 0, na.rm = TRUE)
          n_na <- sum(is.na(pred_vals))
          message(paste("Non-zero R² values:", n_nonzero))
          message(paste("NA values:", n_na))
          if (n_nonzero > 0) {
            message("Sample non-zero values:")
            sample_vals <- head(pred_vals[pred_vals > 0 & !is.na(pred_vals)], 3)
            for (i in seq_along(sample_vals)) {
              message(sprintf("  %s: R² = %.3f", names(sample_vals)[i], sample_vals[i]))
            }
          }
        }
      }
      message("=== END PREDICTABILITY DEBUG ===\n")

      pred_vals
    } else NULL
    
    # Handle labels - use numbers if numbered nodes enabled
    if (isTRUE(input$use_numbered_nodes)) {
      labels_arg <- as.character(1:length(varnames))
      # Store the mapping for the legend table
      values$node_legend <- data.frame(
        Number = 1:length(varnames),
        Variable = varnames,
        stringsAsFactors = FALSE
      )
    } else {
      labels_arg <- if (isTRUE(input$show_labels)) varnames else FALSE
      values$node_legend <- NULL
    }
    
    # Add means to labels if requested (only for single networks, not splits)
    if (isTRUE(input$show_means_on_plot) && isTRUE(input$show_node_means) && 
        !identical(labels_arg, FALSE) && is.null(title)) {
      # Get raw means from dat()
      raw_df <- dat()
      if (!is.null(raw_df) && all(varnames %in% colnames(raw_df))) {
        raw_means <- colMeans(raw_df[, varnames, drop = FALSE], na.rm = TRUE)
        if (isTRUE(input$use_numbered_nodes)) {
          # Add means to numbers: "1 (M=5.2)"
          labels_arg <- paste0(labels_arg, "\n(M=", round(raw_means, 1), ")")
        } else {
          # Add means to variable names: "HADS-A (M=5.2)"
          labels_arg <- paste0(labels_arg, "\n(M=", round(raw_means, 1), ")")
        }
      }
    }
    
    if (isTRUE(input$nodeSizeByMean)) {
      means <- if (!is.null(means_override)) means_override else means_for_nodes()
      vsize_vec <- rep(input$vsize, length(varnames))
      
      if (!is.null(means) && !is.list(means) && length(means) > 0) {
        m <- means[varnames]
        
        # Check if we have pre-computed scaled sizes (from global scaling)
        if (!is.null(scaled_sizes_override)) {
          message(paste("=== Using global scaled sizes for", length(varnames), "variables ==="))
          message(paste("varnames:", paste(head(varnames, 5), collapse=", ")))
          message(paste("scaled_sizes_override has", length(scaled_sizes_override), "values"))
          message(paste("scaled_sizes_override names:", paste(head(names(scaled_sizes_override), 5), collapse=", ")))
          
          # Use pre-computed scaled sizes directly (already scaled globally)
          for (i in seq_along(varnames)) {
            var <- varnames[i]
            if (var %in% names(scaled_sizes_override) && is.finite(scaled_sizes_override[var])) {
              vsize_vec[i] <- scaled_sizes_override[var]
              message(paste("  ", var, "→ size", round(vsize_vec[i], 2)))
            } else {
              message(paste("  ", var, "→ NOT FOUND or not finite"))
            }
          }
          message(paste("Final vsize_vec:", paste(round(head(vsize_vec, 5), 2), collapse=", ")))
          message("=== End scaled sizes debug ===")
        } else {
          # Local scaling within this network only
          # Ensure m is finite before indexing
          m[!is.finite(m)] <- NA
          idx_ok <- which(!is.na(m) & is.finite(m))
          
          if (length(idx_ok) > 0) {
            # Apply reverse scaling if checkbox is checked (for change scores)
            use_reverse <- isTRUE(input$useReverseScaling)
            size_range <- input$nodeSizeRange
            if (is.null(size_range) || length(size_range) != 2 || anyNA(size_range)) {
              size_range <- c(4, 10)  # default range
            }
            vsize_vec[idx_ok] <- scale_to_range(m[idx_ok],
                                                rng = size_range,
                                                reverse = use_reverse)
            # Only show message if not using global scaling
            if (use_reverse) {
              message("Node sizing: Using REVERSE min-max scaling (larger reduction → larger node)")
            } else {
              message("Node sizing: Using STANDARD min-max scaling (larger value → larger node)")
            }
          }
        }
      }
      # Final safety check
      vsize_vec[!is.finite(vsize_vec)] <- input$vsize
    } else {
      vsize_vec <- input$vsize
    }

    pieColorArg <- if (!is.null(pies_all) && !is.list(pies_all)) {
      make_pie_color_arg(input$pieColorChoice %||% "#ADD8E6", length(varnames), bg = device_bg)
    } else {
      sanitize_colors(input$pieColorChoice %||% "#ADD8E6", bg = device_bg)[1]
    }
    
    directed <- identical(input$estimator, "relimp")
    
    # Store community info for EGA
    ega_communities <- NULL
    n_communities <- 1

    # Check if communities were passed from split analysis
    if (!is.null(communities_override)) {
      ega_communities <- communities_override$communities
      n_communities <- communities_override$n_communities
      message(paste("Using shared communities from split analysis:", n_communities, "communities"))
    }

    # Determine layout to use
    layout_to_use <- if (!is.null(layout_override)) {
      # Validate layout_override matches number of nodes
      if (is.matrix(layout_override) && nrow(layout_override) == length(varnames)) {
        layout_override
      } else {
        message(paste("WARNING: layout_override has", nrow(layout_override),
                     "rows but network has", length(varnames), "nodes. Using default layout."))
        input$layout %||% "spring"
      }
    } else {
      # Compute special layouts if needed
      if (!is.null(input$layout) && identical(input$layout, "ega") && isTRUE(use_direct_qgraph)) {
        W <- as.matrix(W_try)
        storage.mode(W) <- "double"
        W[is.na(W)] <- 0

        # Get layout and communities with raw data - use only variables in network
        df_network <- df[, varnames, drop = FALSE]
        ega_result <- compute_ega_layout(W, data = df_network, return_communities = TRUE,
                                         separation = input$ega_separation %||% 3)
        if (is.list(ega_result)) {
          layout_to_use <- ega_result$layout
          ega_communities <- ega_result$communities
          n_communities <- ega_result$n_communities
        } else {
          layout_to_use <- ega_result
        }
        layout_to_use
      } else if (!is.null(input$layout) && identical(input$layout, "pca") && isTRUE(use_direct_qgraph)) {
        W <- as.matrix(W_try)
        storage.mode(W) <- "double"
        W[is.na(W)] <- 0
        df_network <- df[, varnames, drop = FALSE]
        pca_lm <- compute_pca_layout(W, data = df_network)
        if (is.matrix(pca_lm)) {
          k_pca <- if (length(varnames) >= 6) 3L else 2L
          km <- tryCatch(kmeans(pca_lm, centers = min(k_pca, nrow(pca_lm) - 1), nstart = 20),
                         error = function(e) NULL)
          if (!is.null(km)) {
            ega_communities <- km$cluster
            n_communities   <- length(unique(km$cluster))
          }
        }
        pca_lm
      } else {
        input$layout %||% "spring"
      }
    }
    if (is.matrix(layout_to_use)) main_net_layout(layout_to_use)

    # Determine node colors - use communities ONLY if multiple communities detected
    node_colors <- if (!is.null(ega_communities) &&
                       isTRUE(input$ega_color_communities) &&
                       !is.na(n_communities) &&
                       n_communities > 1) {
      # Create pastel color palette for communities
      # Lower saturation (s=0.4) and higher value (v=0.95) creates softer, more pastel colors
      comm_colors <- rainbow(n_communities, s = 0.4, v = 0.95, alpha = 0.85)
      comm_colors[ega_communities]
    } else {
      input$nodeColor
    }

    # Determine node shapes based on variable types
    # circle = continuous/count, square = categorical (3+ levels), triangle = dichotomous (2 levels)
    node_shapes <- sapply(varnames, function(var) {
      if (var %in% disc_vars()) {
        # Check if dichotomous (2 levels) or polytomous (3+ levels)
        n_levels <- length(unique(na.omit(df[[var]])))
        if (n_levels == 2) {
          "triangle"  # Dichotomous
        } else {
          "square"    # Categorical (3+ levels)
        }
      } else if (var %in% count_vars()) {
        "circle"      # Count variables use circles
      } else {
        "circle"      # Continuous variables use circles (default)
      }
    })

    # Debug output for node shapes
    message("\n=== NODE SHAPES DEBUG ===")
    message(paste("Total variables:", length(varnames)))
    message(paste("Continuous vars:", paste(cont_vars(), collapse=", ")))
    message(paste("Categorical vars:", paste(disc_vars(), collapse=", ")))
    message(paste("Count vars:", paste(count_vars(), collapse=", ")))
    shape_summary <- table(node_shapes)
    message("Shape distribution:")
    for (shape in names(shape_summary)) {
      message(sprintf("  %s: %d nodes", shape, shape_summary[shape]))
    }
    # Show which variables get which shapes
    non_circle <- node_shapes[node_shapes != "circle"]
    if (length(non_circle) > 0) {
      message("Non-circle shapes:")
      for (i in seq_along(non_circle)) {
        message(sprintf("  %s → %s", names(non_circle)[i], non_circle[i]))
      }
    }
    message("=== END NODE SHAPES DEBUG ===\n")

    if (use_direct_qgraph) {
      W <- as.matrix(W_try)
      storage.mode(W) <- "double"
      W[is.na(W)] <- 0
      if (is.null(colnames(W)) || is.null(rownames(W))) {
        colnames(W) <- rownames(W) <- varnames
      }
      if (setequal(colnames(W), varnames)) {
        W <- W[varnames, varnames, drop = FALSE]
      }
      
      minS <- suppressWarnings(as.numeric(input$edgeScaleMin))
      if (!is.finite(minS) || minS < 0) minS <- 0
      maxS <- suppressWarnings(as.numeric(input$edgeScaleMax))
      if (!is.finite(maxS) || maxS <= 0) maxS <- 1
      Wviz <- rescale_edge_weights(W, minS, maxS)

      # Extract pies using varnames (original variable names) not colnames(Wviz)
      # because pcor may abbreviate names in the network matrix
      pies <- if (!is.null(pies_all) && !is.list(pies_all)) pies_all[varnames] else NULL

      # Debug pie extraction
      message("\n=== PIE EXTRACTION DEBUG ===")
      message(paste("pies_all is.null:", is.null(pies_all)))
      if (!is.null(pies_all)) {
        message(paste("pies_all length:", length(pies_all)))
        message(paste("pies_all names:", paste(head(names(pies_all), 5), collapse=", ")))
        message(paste("Wviz colnames:", paste(head(colnames(Wviz), 5), collapse=", ")))
        message(paste("Names match?:", all(colnames(Wviz) %in% names(pies_all))))
        message(paste("Wviz has names?:", !is.null(colnames(Wviz))))
      }
      message(paste("pies (after extraction) is.null:", is.null(pies)))
      if (!is.null(pies)) {
        message(paste("pies length:", length(pies)))
        message(paste("pies with values > 0 (BEFORE sanitization):", sum(pies > 0, na.rm = TRUE)))
        message(paste("pies with NA (BEFORE sanitization):", sum(is.na(pies))))
        message(paste("Sample pies values (first 3):", paste(head(pies, 3), collapse=", ")))
      }
      message("=== END PIE EXTRACTION ===\n")

      # Sanitize pie values - remove NA values and ensure they're in [0, 1] range
      if (!is.null(pies) && length(pies) > 0) {
        # Replace NA values with 0
        pies[is.na(pies)] <- 0
        # Ensure all values are in [0, 1] range
        pies[pies < 0] <- 0
        pies[pies > 1] <- 1
        # If all pies are 0, set to NULL (no predictability to show)
        if (all(pies == 0)) {
          message("WARNING: All pies are 0 after sanitization, setting to NULL")
          pies <- NULL
        }
      }

      # Update pieColorArg based on network size (must match number of nodes)
      # pieColor in qgraph must be either length 1 or equal to ncol(Wviz)
      if (!is.null(pies) && length(pies) > 0) {
        # Use single color (length 1) - simpler and always works
        pieColorArg <- sanitize_colors(input$pieColorChoice %||% "#ADD8E6", bg = device_bg)[1]
        message(paste("pieColorArg set to single color:", pieColorArg))
      } else {
        pieColorArg <- NULL
        message("pieColorArg set to NULL (no pies)")
      }

      # Sanitize input parameters to avoid NA/NULL issues
      min_edge_vis <- input$minEdgeVis
      if (is.null(min_edge_vis) || !is.finite(min_edge_vis)) min_edge_vis <- 0
      
      esize_val <- input$esize
      if (is.null(esize_val) || !is.finite(esize_val)) esize_val <- 1.5
      
      asize_val <- input$asize
      if (is.null(asize_val) || !is.finite(asize_val)) asize_val <- 5
      
      label_cex_val <- input$label_cex
      if (is.null(label_cex_val) || !is.finite(label_cex_val)) label_cex_val <- 1
      
      edge_label_cex_val <- input$edgeLabelCex
      if (is.null(edge_label_cex_val) || !is.finite(edge_label_cex_val)) edge_label_cex_val <- 0.5

      # Sanitize theme
      theme_val <- input$theme
      if (is.null(theme_val) || !is.character(theme_val)) theme_val <- "classic"

      # Sanitize node colors
      if (is.null(node_colors) || length(node_colors) == 0) {
        node_colors <- input$nodeColor %||% "#FFFFFF"
      }

      # Debug parameter values before qgraph call
      message("=== QGRAPH PARAMETERS DEBUG ===")
      message(paste("layout_to_use class:", paste(class(layout_to_use), collapse=", ")))
      message(paste("layout_to_use is.null:", is.null(layout_to_use)))
      if (is.character(layout_to_use)) message(paste("layout string:", layout_to_use))
      message(paste("theme:", input$theme))
      message(paste("directed:", isTRUE(directed)))
      message(paste("label_scale_equal:", isTRUE(input$label_scale_equal)))
      message(paste("fade:", TRUE))
      message(paste("minimum:", min_edge_vis))
      # Check if there are any edges to display
      n_edges <- sum(abs(Wviz) > 1e-10, na.rm = TRUE) / 2  # Divide by 2 for undirected, or just count for directed
      message(paste("Number of edges:", n_edges))
      message(paste("edge.labels requested:", isTRUE(input$edgeLabels)))
      message(paste("vsize_vec length:", length(vsize_vec)))
      message(paste("vsize_vec range:", paste(range(vsize_vec, na.rm=TRUE), collapse=" to ")))
      message(paste("node_colors length:", length(node_colors)))
      message("=== END DEBUG ===")

      # Only show edge labels if there are edges AND user requested them
      show_edge_labels <- isTRUE(input$edgeLabels) && n_edges > 0

      # Safe qgraph call with error handling
      tryCatch({
        suppressWarnings(qgraph::qgraph(
          Wviz,
          layout = layout_to_use,
          theme = theme_val,
          directed = isTRUE(directed), parallelEdge = TRUE,
          pie = pies, pieColor = pieColorArg, pieBorder = input$pieBorder %||% 0.3,
          color = node_colors,
          shape = node_shapes,
          border.color = input$nodeBorderColor,
          legend = FALSE,
          labels = labels_arg,
          label.cex = label_cex_val * 0.7,
          label.font = if (isTRUE(input$label_bold)) 2 else 1,
          label.scale.equal = isTRUE(input$label_scale_equal),
          vsize = vsize_vec, esize = esize_val, asize = asize_val,
          fade = TRUE, minimum = min_edge_vis,
          edge.labels = show_edge_labels,
          edge.label.cex = edge_label_cex_val,
          title = title,
          title.cex = 1.2
        ))
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Plotting error:", e$message), cex = 0.8)
        message(paste("ERROR in qgraph:", e$message))
        message(paste("Full error:", paste(capture.output(print(e)), collapse="\n")))
      })
    } else {
      pies <- if (!is.null(pies_all) && !is.list(pies_all)) pies_all[varnames] else NULL

      # Sanitize pie values - remove NA values and ensure they're in [0, 1] range
      if (!is.null(pies) && length(pies) > 0) {
        # Replace NA values with 0
        pies[is.na(pies)] <- 0
        # Ensure all values are in [0, 1] range
        pies[pies < 0] <- 0
        pies[pies > 1] <- 1
        # If all pies are 0, set to NULL (no predictability to show)
        if (all(pies == 0)) pies <- NULL
      }

      # Update pieColorArg based on final pies value
      # pieColor in qgraph must be either length 1 or equal to number of nodes
      if (!is.null(pies) && length(pies) > 0) {
        # Use single color (length 1) - simpler and always works
        pieColorArg <- sanitize_colors(input$pieColorChoice %||% "#ADD8E6", bg = device_bg)[1]
      } else {
        pieColorArg <- NULL
      }

      # Sanitize input parameters to avoid NA/NULL issues
      min_edge_vis <- input$minEdgeVis
      if (is.null(min_edge_vis) || !is.finite(min_edge_vis)) min_edge_vis <- 0
      
      esize_val <- input$esize
      if (is.null(esize_val) || !is.finite(esize_val)) esize_val <- 1.5
      
      asize_val <- input$asize
      if (is.null(asize_val) || !is.finite(asize_val)) asize_val <- 5
      
      label_cex_val <- input$label_cex
      if (is.null(label_cex_val) || !is.finite(label_cex_val)) label_cex_val <- 1
      
      edge_label_cex_val <- input$edgeLabelCex
      if (is.null(edge_label_cex_val) || !is.finite(edge_label_cex_val)) edge_label_cex_val <- 0.5

      # Sanitize theme
      theme_val <- input$theme
      if (is.null(theme_val) || !is.character(theme_val)) theme_val <- "classic"

      # Sanitize node colors
      if (is.null(node_colors) || length(node_colors) == 0) {
        node_colors <- input$nodeColor %||% "#FFFFFF"
      }

      # Check if there are edges to label
      W_check <- tryCatch(getWmat(net), error = function(e) NULL)
      n_edges <- if (!is.null(W_check)) sum(abs(W_check) > 1e-10, na.rm = TRUE) / 2 else 0
      show_edge_labels <- isTRUE(input$edgeLabels) && n_edges > 0

      # Safe plot call with error handling
      tryCatch({
        suppressWarnings(plot(
          net,
          layout = layout_to_use,
          theme = theme_val,
          directed = isTRUE(directed), parallelEdge = TRUE,
          pie = pies, pieColor = pieColorArg, pieBorder = input$pieBorder %||% 0.3,
          color = node_colors,
          shape = node_shapes,
          border.color = input$nodeBorderColor,
          legend = FALSE,
          labels = labels_arg,
          label.cex = label_cex_val * 0.7,
          label.font = if (isTRUE(input$label_bold)) 2 else 1,
          label.scale.equal = isTRUE(input$label_scale_equal),
          vsize = vsize_vec, esize = esize_val, asize = asize_val,
          fade = TRUE, minimum = min_edge_vis,
          edge.labels = show_edge_labels,
          edge.label.cex = edge_label_cex_val,
          title = title,
          title.cex = 1.2
        ))
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Plotting error:", e$message), cex = 1)
      })
    }
    
    invisible(NULL)
  }
  
  output$split_info <- renderPrint({
    res <- est_result_thr()
    if (!is.null(res) && res$split) {
      df <- dat()
      split_var <- res$split_var
      if (!is.null(df) && split_var %in% colnames(df)) {
        split_table <- table(df[[split_var]])
        cat("Split by:", split_var, "\n")
        cat("Group sizes:\n")
        print(split_table)
      }
    }
  })
  
  output$netSummary <- renderPrint({
    res <- est_result_thr()
    if (is.null(res)) return(invisible(NULL))
    
    # Display layout info
    cat("Layout:", input$layout, "\n")
    if (input$layout == "ega") {
      if (isTRUE(input$ega_color_communities)) {
        cat("  Community coloring: ENABLED\n")
      } else {
        cat("  Community coloring: disabled\n")
      }
      cat("  Note: Check console/messages for number of communities detected\n")
      cat("        (Coloring only applies when >1 community found)\n")
    }
    cat("\n")
    
    if (res$split) {
      cat("Network comparison across", length(res$groups), "groups\n")
      cat("Split variable:", res$split_var, "\n\n")
      for (group_name in names(res$groups)) {
        net <- res$groups[[group_name]]
        if (!is.null(net)) {
          cat("Group:", group_name, "\n")
          W <- tryCatch(bootnet::getWmat(net), error = function(e) NULL)
          if (!is.null(W)) {
            cat("  Graph dim:", paste(dim(W), collapse=" x "), "\n")
          }
        }
      }
    } else {
      net <- res$network
      W <- tryCatch(bootnet::getWmat(net), error = function(e) NULL)
      if (!is.null(W)) {
        if (!is.null(colnames(W))) {
          cat("Nodes:\n", paste(colnames(W), collapse=", "), "\n")
        }
        cat("\nGraph dim:", paste(dim(W), collapse=" x "), "\n")
      }
      cat("\nNode predictability (R²):\n")
      pies <- predictability_vec()
      df <- prepped_data()
      message("DEBUG: pies is.null:", is.null(pies))
      message("DEBUG: df is.null:", is.null(df))
      if (!is.null(pies)) {
        message("DEBUG: pies class:", class(pies))
        message("DEBUG: pies length:", length(pies))
        if (length(pies) > 0) {
          message("DEBUG: First few pies:", paste(head(names(pies)), "=", round(head(pies), 3), collapse=", "))
        }
      }
      if (!is.null(pies) && !is.null(df)) {
        out <- data.frame(Node = colnames(df), R2 = round(pies[colnames(df)], 3))
        print(out, row.names = FALSE)
      } else {
        cat("Predictability unavailable.\n")
      }
    }
  })
  
  # Node legend table for numbered nodes
  output$node_legend_table <- DT::renderDataTable({
    req(input$use_numbered_nodes)
    req(values$node_legend)
    
    DT::datatable(
      values$node_legend,
      options = list(
        pageLength = 20,
        searching = TRUE,
        ordering = FALSE,
        dom = 'ftp',
        columnDefs = list(
          list(width = '80px', targets = 0),
          list(width = 'auto', targets = 1)
        )
      ),
      rownames = FALSE,
      class = 'cell-border stripe',
      escape = FALSE
    )
  })
  
  # Node means table
  output$node_means_table <- DT::renderDataTable({
    req(input$show_node_means)
    req(values$node_means)
    
    means_df <- values$node_means
    
    # Determine if this is split analysis (has group columns)
    has_groups <- any(grepl("^Group_", colnames(means_df)))
    
    if (has_groups) {
      # Split analysis - show groups + ANOVA
      DT::datatable(
        means_df,
        options = list(
          pageLength = 20,
          searching = TRUE,
          ordering = TRUE,
          dom = 'ftp',
          columnDefs = list(
            list(className = 'dt-center', targets = 1:(ncol(means_df)-1))
          )
        ),
        rownames = FALSE,
        class = 'cell-border stripe',
        escape = FALSE
      ) %>%
        DT::formatRound(columns = setdiff(colnames(means_df), "Variable"), digits = 3) %>%
        DT::formatStyle(
          'p_value',
          backgroundColor = DT::styleInterval(
            cuts = c(0.001, 0.01, 0.05),
            values = c('#90EE90', '#FFFF99', '#FFE4B5', 'white')
          )
        )
    } else {
      # Single network - show overall means
      DT::datatable(
        means_df,
        options = list(
          pageLength = 20,
          searching = TRUE,
          ordering = TRUE,
          order = list(list(1, 'desc')),  # Sort by Mean descending
          dom = 'ftp',
          columnDefs = list(
            list(className = 'dt-center', targets = 1:2)
          )
        ),
        rownames = FALSE,
        class = 'cell-border stripe',
        escape = FALSE
      ) %>%
        DT::formatRound(columns = c('Mean', 'SD'), digits = 3)
    }
  })
  
  # Split predictability table
  output$split_predictability_table <- DT::renderDataTable({
    req(values$split_predictability)

    DT::datatable(
      values$split_predictability,
      options = list(
        pageLength = 20,
        searching = TRUE,
        ordering = TRUE,
        dom = 'ftp',
        columnDefs = list(
          list(className = 'dt-center', targets = 1:(ncol(values$split_predictability)-1))
        )
      ),
      rownames = FALSE,
      class = 'cell-border stripe',
      escape = FALSE
    ) %>%
      DT::formatRound(columns = 2:ncol(values$split_predictability), digits = 3)
  })

  # Edge table
  output$edge_table <- DT::renderDataTable({
    req(input$show_edge_table)
    res <- est_result_thr()
    req(res)

    # Function to extract edges from a weight matrix
    extract_edges <- function(W, group_name = NULL) {
      if (is.null(W)) return(NULL)

      W <- as.matrix(W)
      node_names <- colnames(W)
      if (is.null(node_names)) {
        node_names <- paste0("V", 1:ncol(W))
      }

      # Get upper triangle indices (avoid duplicates)
      edges_list <- list()
      for (i in 1:(nrow(W)-1)) {
        for (j in (i+1):ncol(W)) {
          weight <- W[i, j]
          # Only include non-zero edges
          if (!is.na(weight) && abs(weight) > 1e-10) {
            edge_info <- list(
              From = node_names[i],
              To = node_names[j],
              Weight = weight
            )
            if (!is.null(group_name)) {
              edge_info$Group <- group_name
            }
            edges_list[[length(edges_list) + 1]] <- edge_info
          }
        }
      }

      if (length(edges_list) == 0) return(NULL)
      do.call(rbind, lapply(edges_list, as.data.frame, stringsAsFactors = FALSE))
    }

    # Check if split analysis
    if (res$split && !is.null(res$groups)) {
      # Combine edges from all groups
      all_edges <- list()
      for (group_name in names(res$groups)) {
        net <- res$groups[[group_name]]
        W <- safe_get_W(net)
        if (!is.null(W)) {
          group_edges <- extract_edges(W, group_name = group_name)
          if (!is.null(group_edges)) {
            all_edges[[group_name]] <- group_edges
          }
        }
      }

      if (length(all_edges) == 0) {
        return(DT::datatable(data.frame(Message = "No edges found"), rownames = FALSE))
      }

      edges_df <- do.call(rbind, all_edges)

      DT::datatable(
        edges_df,
        options = list(
          pageLength = 25,
          searching = TRUE,
          ordering = TRUE,
          order = list(list(3, 'desc')),  # Sort by Weight (absolute) descending
          dom = 'ftp'
        ),
        rownames = FALSE,
        class = 'cell-border stripe',
        escape = FALSE
      ) %>%
        DT::formatRound(columns = 'Weight', digits = 4)

    } else {
      # Single network
      net <- res$network
      W <- safe_get_W(net)

      if (is.null(W)) {
        return(DT::datatable(data.frame(Message = "No weight matrix available"), rownames = FALSE))
      }

      edges_df <- extract_edges(W)

      if (is.null(edges_df)) {
        return(DT::datatable(data.frame(Message = "No edges found"), rownames = FALSE))
      }

      DT::datatable(
        edges_df,
        options = list(
          pageLength = 25,
          searching = TRUE,
          ordering = TRUE,
          order = list(list(2, 'desc')),  # Sort by Weight descending
          dom = 'ftp'
        ),
        rownames = FALSE,
        class = 'cell-border stripe',
        escape = FALSE
      ) %>%
        DT::formatRound(columns = 'Weight', digits = 4)
    }
  })

  output$centPlot <- renderPlot({
    res <- est_result_thr()
    shiny::validate(shiny::need(!is.null(res), "Estimate a network first."))
    if (res$split) {
      n_groups <- length(res$groups)
      if (n_groups == 0) return()
      cent_data <- list()
      for (group_name in names(res$groups)) {
        net <- res$groups[[group_name]]
        if (!is.null(net)) {
          W <- safe_get_W(net)
          if (!is.null(W)) {
            W <- as.matrix(W)
            storage.mode(W) <- "double"
            W[!is.finite(W)] <- 0
            df <- prepped_data()
            if (is.null(colnames(W)) && !is.null(df)) {
              colnames(W) <- rownames(W) <- colnames(df)
            }
            if (sum(abs(W)) > 0) {
              ct <- compute_centrality_measures(W, measures = input$cent_include)
              ct$Group <- group_name
              cent_data[[group_name]] <- ct
            }
          }
        }
      }
      if (length(cent_data) == 0) {
        plot.new()
        text(0.5, 0.5, "No valid centrality data", cex = 2)
        return()
      }
      all_cent <- do.call(rbind, cent_data)
      measures <- setdiff(names(all_cent), c("Node", "Group"))
      nm <- length(measures)
      if (nm == 0) return()
      ncol <- min(2, nm)
      nrow <- ceiling(nm / ncol)
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      lbl_chars <- max(nchar(unique(all_cent$Node)), na.rm = TRUE)
      left_mar  <- max(8, ceiling(lbl_chars * 0.6))
      par(mfrow = c(nrow, ncol), mar = c(4, left_mar, 3, 1))
      group_colors <- c("#5B9BD5", "#ED7D31", "#70AD47", "#FFC000", 
                        "#9E7CC1", "#4BACC6", "#F68B8B", "#C6AC65")
      for (m in measures) {
        cent_matrix <- tapply(all_cent[[m]], list(all_cent$Node, all_cent$Group), mean)
        cent_matrix[is.na(cent_matrix)] <- 0
        n_groups <- ncol(cent_matrix)
        colors_to_use <- group_colors[1:min(n_groups, length(group_colors))]
        barplot(t(cent_matrix),
                beside = TRUE,
                horiz = TRUE,
                las = 1,
                main = m,
                legend.text = colnames(cent_matrix),
                args.legend = list(x = "topright", bty = "n", cex = 0.8),
                col = colors_to_use)
      }
    } else {
      net <- res$network
      W <- safe_get_W(net)
      shiny::validate(shiny::need(!is.null(W), "Centrality needs a weight matrix."))
      W <- as.matrix(W)
      storage.mode(W) <- "double"
      W[!is.finite(W)] <- 0
      df <- prepped_data()
      if (is.null(colnames(W)) && !is.null(df)) {
        colnames(W) <- rownames(W) <- colnames(df)
      }
      shiny::validate(shiny::need(!is.null(colnames(W)), "No variable names."))
      shiny::validate(shiny::need(sum(abs(W)) > 0, "All edges are zero."))
      include_vec <- input$cent_include
      ct <- compute_centrality_measures(W, measures = include_vec)
      node_col <- "Node"
      nodes <- ct[[node_col]]
      measures <- setdiff(names(ct), node_col)
      ord_measure <- input$cent_orderBy
      if (!is.null(ord_measure) && ord_measure %in% measures) {
        ord_idx <- order(ct[[ord_measure]], decreasing = isTRUE(input$cent_decreasing))
        ct <- ct[ord_idx, , drop = FALSE]
        nodes <- nodes[ord_idx]
      }
      nm <- length(measures)
      shiny::validate(shiny::need(nm > 0, "No centrality measures selected."))
      ncol <- min(2, nm)
      nrow <- ceiling(nm / ncol)
      op <- par(no.readonly = TRUE)
      on.exit(par(op), add = TRUE)
      lbl_chars <- max(nchar(nodes), na.rm = TRUE)
      left_mar  <- max(8, ceiling(lbl_chars * 0.6))
      par(mfrow = c(nrow, ncol), mar = c(4, left_mar, 3, 1))
      for (m in measures) {
        vals <- ct[[m]]
        barplot(height = vals,
                names.arg = nodes,
                horiz = TRUE,
                las = 1,
                cex.names = max(0.4, min(1.2, input$label_cex * 0.4)),
                main = m,
                border = NA,
                col = "#5B9BD5",  # Soft blue color
                xlim = c(0, max(vals) * 1.1))
        box(bty = "n")
      }
    }
  })
  
  output$centWarn <- renderPrint({
    if (!identical(input$estimator, "relimp") &&
        any(c("InDegree","OutDegree","InStrength","OutStrength") %in% input$cent_include)) {
      cat("Note: In/Out metrics are for directed networks.\n")
    }
  })
  
  output$downloadRDS <- downloadHandler(
    filename = function(){
      paste0("bootnet_results_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".rds")
    },
    content = function(file){
      # Capture reactive values safely
      tryCatch({
        message("=== Download RDS: Starting ===")

        # Get estimation results
        est_res <- isolate(est_result_thr())
        if (is.null(est_res)) {
          showNotification("No network estimated yet. Please run an estimation first.",
                          type = "error", duration = 5)
          return(NULL)
        }

        message("Estimation result captured")

        # Build output object
        obj <- list(
          raw_dat = isolate(raw_dat()),
          dat = isolate(dat()),
          classification = list(
            continuous = isolate(cont_vars()),
            discrete = isolate(disc_vars()),
            count = isolate(count_vars()),
            split = isolate(split_var())
          ),
          estimator = isolate(input$estimator),
          estimate = est_res,
          predictability = isolate(predictability_vec()),
          viz = list(
            layout = isolate(input$layout),
            theme = isolate(input$theme),
            legend = isolate(input$show_legend),
            labels = isolate(input$show_labels),
            vsize = isolate(input$vsize),
            esize = isolate(input$esize),
            asize = isolate(input$asize),
            label_cex = isolate(input$label_cex),
            label_scale_equal = isolate(input$label_scale_equal),
            nodeSizeByMean = isolate(input$nodeSizeByMean),
            nodeSizeRange = isolate(input$nodeSizeRange),
            edgeLabels = isolate(input$edgeLabels),
            edgeLabelCex = isolate(input$edgeLabelCex),
            edgeScaleMin = isolate(input$edgeScaleMin),
            edgeScaleMax = isolate(input$edgeScaleMax),
            minEdgeVis = isolate(input$minEdgeVis),
            pieColorChoice = isolate(input$pieColorChoice),
            pieBorder = isolate(input$pieBorder),
            nodeColor = isolate(input$nodeColor),
            nodeBorderColor = isolate(input$nodeBorderColor)
          ),
          centrality = list(
            include = isolate(input$cent_include),
            orderBy = isolate(input$cent_orderBy),
            decreasing = isolate(input$cent_decreasing)
          ),
          bridge = isolate(tryCatch(bridge_result(), error = function(e) NULL)),
          bridge_stability = isolate(tryCatch(bridge_stability_result(), error = function(e) NULL)),
          relimp = list(
            smart_preprocess = isolate(input$relimp_smart_preprocess),
            dummy = isolate(input$relimp_dummy),
            drop_nonfinite = isolate(input$relimp_drop_nonfinite),
            drop_zerovar = isolate(input$relimp_drop_zerovar),
            threshold = isolate(input$relimpEstMin)
          ),
          timestamp = Sys.time(),
          session_info = list(
            R_version = R.version.string,
            bootnet_version = tryCatch(as.character(packageVersion("bootnet")), error = function(e) "unknown")
          )
        )

        message(paste("Object size:", object.size(obj), "bytes"))
        message("Saving to file...")

        saveRDS(obj, file)

        message("=== Download RDS: Complete ===")

        showNotification("Results saved successfully!", type = "message", duration = 3)

      }, error = function(e) {
        message(paste("ERROR in downloadRDS:", e$message))
        showNotification(paste("Error saving results:", e$message),
                        type = "error", duration = 10)
      })
    }
  )

  output$downloadRscript <- downloadHandler(
    filename = function() {
      paste0("BootRIN_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R")
    },
    content = function(file) {
      tryCatch({
        res <- isolate(est_result_thr())
        if (is.null(res)) {
          showNotification("No network estimated yet. Please run estimation first.",
                           type = "error", duration = 5)
          return(NULL)
        }
        code <- generate_r_code(
          input      = isolate(input),
          res        = res,
          boot_res   = isolate(boot_sig_result()),
          df_prepped = isolate(prepped_data()),
          df_raw     = isolate(dat())
        )
        writeLines(code, file)
        showNotification("R script exported successfully!", type = "message", duration = 3)
      }, error = function(e) {
        message(paste("ERROR in downloadRscript:", e$message))
        showNotification(paste("Error exporting script:", e$message),
                         type = "error", duration = 10)
      })
    }
  )

  # ============================================================
  # NETWORK COMPARISON TEST (NCT)
  # ============================================================

  # Reactive to store NCT results
  nct_results <- reactiveVal(NULL)

  # Dynamic UI for subject ID selector
  output$nct_subject_id_selector <- renderUI({
    req(dat())
    full_df <- dat()

    # Get all column names from full dataset
    all_cols <- colnames(full_df)

    # Exclude split variable and network variables from choices
    exclude_cols <- c(split_var(), cont_vars(), disc_vars(), count_vars())
    id_choices <- setdiff(all_cols, exclude_cols)

    selectInput("nct_subject_id", "Subject ID column (optional)",
               choices = c("None (data pre-sorted)" = "", id_choices),
               selected = "")
  })

  # Output to control UI visibility
  output$nct_available <- reactive({
    !is.null(nct_results())
  })
  outputOptions(output, "nct_available", suspendWhenHidden = FALSE)

  # Run NCT when button is clicked
  observeEvent(input$run_nct, {
    req(dat())

    # Check if split variable exists
    split <- split_var()
    if (is.null(split) || split == "") {
      showNotification("Please select a split variable first", type = "error", duration = 5)
      return(NULL)
    }

    # Get estimation results
    est_res <- est_result_thr()
    if (is.null(est_res) || !est_res$split) {
      showNotification("Please run network estimation with a split variable first",
                      type = "error", duration = 5)
      return(NULL)
    }

    # Check estimator compatibility
    current_estimator <- input$estimator
    compatible_estimators <- c("EBICglasso", "IsingFit", "cor", "pcor")

    if (!current_estimator %in% compatible_estimators) {
      showNotification(
        paste0("NCT is not compatible with the '", current_estimator, "' estimator. ",
              "Compatible estimators: EBICglasso, IsingFit, Correlation, Partial Correlation. ",
              "Please re-estimate your network with a compatible estimator."),
        type = "warning", duration = 10
      )
      return(NULL)
    }

    # Check if we have exactly 2 groups
    groups <- names(est_res$groups)
    if (length(groups) != 2) {
      showNotification(
        paste("Network Comparison Test requires exactly 2 groups.",
             "Your split variable has", length(groups), "groups.",
             "Please use a binary split variable."),
        type = "error", duration = 8
      )
      return(NULL)
    }

    # Show progress
    showNotification("Running Network Comparison Test... This may take a while.",
                    type = "message", duration = NULL, id = "nct_progress")

    tryCatch({
      # Get data for both groups
      full_df <- dat()
      group1_name <- groups[1]
      group2_name <- groups[2]

      group1_idx <- which(full_df[[split]] == group1_name)
      group2_idx <- which(full_df[[split]] == group2_name)

      # Get variable names from the estimated network (not all columns!)
      # Use the network matrix column names which represent the actual variables used
      group1_network <- est_res$groups[[group1_name]]
      varnames <- colnames(group1_network$graph)

      if (is.null(varnames) || length(varnames) == 0) {
        removeNotification("nct_progress")
        showNotification("Cannot extract variable names from estimated network",
                        type = "error", duration = 5)
        return(NULL)
      }

      message(paste("Using", length(varnames), "variables from estimated network"))
      message(paste("Variables:", paste(varnames, collapse = ", ")))

      # Determine if paired comparison
      is_paired <- !is.null(input$nct_test_type) && input$nct_test_type == "paired"

      # Handle paired comparison with optional subject ID matching
      if (is_paired) {
        subject_id_col <- input$nct_subject_id

        if (!is.null(subject_id_col) && subject_id_col != "") {
          # Subject ID provided - match subjects across groups
          message(paste("\n=== Paired NCT with Subject ID matching ==="))
          message(paste("Subject ID column:", subject_id_col))

          # Get subject IDs for both groups
          data1_full <- full_df[group1_idx, c(subject_id_col, varnames), drop = FALSE]
          data2_full <- full_df[group2_idx, c(subject_id_col, varnames), drop = FALSE]

          # Find common subjects
          subjects1 <- data1_full[[subject_id_col]]
          subjects2 <- data2_full[[subject_id_col]]
          common_subjects <- intersect(subjects1, subjects2)

          if (length(common_subjects) == 0) {
            removeNotification("nct_progress")
            showNotification(
              "No matching subjects found between groups. Check your subject ID column.",
              type = "error", duration = 10
            )
            return(NULL)
          }

          message(paste("Found", length(common_subjects), "matching subjects"))
          message(paste("Group 1 had", length(subjects1), "total subjects"))
          message(paste("Group 2 had", length(subjects2), "total subjects"))

          # Match and sort by subject ID
          data1_matched <- data1_full[match(common_subjects, subjects1), varnames, drop = FALSE]
          data2_matched <- data2_full[match(common_subjects, subjects2), varnames, drop = FALSE]

          data1 <- data1_matched
          data2 <- data2_matched

        } else {
          # No subject ID - assume data is pre-sorted
          message(paste("\n=== Paired NCT (pre-sorted data) ==="))

          data1 <- full_df[group1_idx, varnames, drop = FALSE]
          data2 <- full_df[group2_idx, varnames, drop = FALSE]

          if (nrow(data1) != nrow(data2)) {
            removeNotification("nct_progress")
            showNotification(
              paste0("Paired NCT requires equal sample sizes. ",
                    "Group 1: ", nrow(data1), " observations. ",
                    "Group 2: ", nrow(data2), " observations. ",
                    "Either select a Subject ID column or ensure data is pre-sorted with equal sizes."),
              type = "error", duration = 10
            )
            return(NULL)
          }

          message(paste("⚠️ Assuming rows are pre-sorted and matched across groups"))
        }
      } else {
        # Unpaired - no matching needed
        data1 <- full_df[group1_idx, varnames, drop = FALSE]
        data2 <- full_df[group2_idx, varnames, drop = FALSE]
      }

      # Data is already numeric from network estimation, just ensure it's in the right format
      # Convert to matrix for NCT
      data1_numeric <- as.data.frame(data1)
      data2_numeric <- as.data.frame(data2)

      message("\n=== Running NCT ===")
      message(paste("Group 1 (", group1_name, "):", nrow(data1_numeric), "observations"))
      message(paste("Group 2 (", group2_name, "):", nrow(data2_numeric), "observations"))
      message(paste("Variables:", length(varnames)))
      message(paste("Iterations:", input$nct_iterations))
      message(paste("Test edges:", input$nct_test_edges))
      message(paste("Test centrality:", input$nct_test_centrality))

      # Determine if Ising network (binary data)
      is_ising <- current_estimator == "IsingFit"

      # Get gamma parameter from user's estimation settings
      nct_gamma <- if (current_estimator == "EBICglasso") {
        input$tuning %||% 0.5
      } else if (current_estimator == "IsingFit") {
        input$ising_gamma %||% 0.25
      } else {
        0.5  # Default for cor/pcor
      }

      # Get AND rule for Ising
      nct_and <- if (is_ising) {
        input$ising_and %||% TRUE
      } else {
        TRUE
      }

      message(paste("Estimator:", current_estimator))
      message(paste("NCT comparison type:", if(is_paired) "Paired" else "Independent"))
      message(paste("Binary data (Ising):", is_ising))
      message(paste("EBIC gamma:", nct_gamma))
      message(paste("AND-rule:", nct_and))

      # Map user-selected centrality measures to NCT format
      # NCT accepts: "strength", "closeness", "betweenness", "expectedInfluence"
      user_centrality <- input$cent_include %||% c("Strength")
      nct_centrality <- c()

      # Map from app names to NCT names
      if ("Strength" %in% user_centrality) nct_centrality <- c(nct_centrality, "strength")
      if ("Closeness" %in% user_centrality) nct_centrality <- c(nct_centrality, "closeness")
      if ("Betweenness" %in% user_centrality) nct_centrality <- c(nct_centrality, "betweenness")
      if ("ExpectedInfluence" %in% user_centrality) nct_centrality <- c(nct_centrality, "expectedInfluence")

      # If no valid centrality measures selected or centrality test disabled, use default
      if (length(nct_centrality) == 0 || !input$nct_test_centrality) {
        nct_centrality <- c("strength", "expectedInfluence")
      }

      message(paste("Centrality measures for NCT:", paste(nct_centrality, collapse = ", ")))

      # Run NCT
      # Note: NCT always tests network invariance and global strength
      # test.edges and test.centrality are optional additional tests
      nct_result <- NetworkComparisonTest::NCT(
        data1 = data1_numeric,
        data2 = data2_numeric,
        gamma = nct_gamma,  # Use gamma from user's estimation settings
        it = input$nct_iterations %||% 1000,
        binary.data = is_ising,  # TRUE for IsingFit, FALSE otherwise
        paired = is_paired,
        weighted = TRUE,
        AND = nct_and,  # Use AND from Ising settings if applicable
        test.edges = input$nct_test_edges %||% FALSE,
        edges = "all",
        progressbar = FALSE,
        test.centrality = input$nct_test_centrality %||% FALSE,
        centrality = nct_centrality,  # Use user-selected centrality measures
        verbose = TRUE
      )

      # Debug: Show NCT result structure
      message("\n=== NCT Result Structure ===")
      message(paste("All NCT fields:", paste(names(nct_result), collapse = ", ")))
      message(paste("nw.pval class:", class(nct_result$nw.pval)))
      message(paste("nw.pval value:", nct_result$nw.pval))
      message(paste("glstrinv.real class:", class(nct_result$glstrinv.real)))
      message(paste("glstrinv.real length:", length(nct_result$glstrinv.real)))
      message(paste("glstrinv.real values:", paste(nct_result$glstrinv.real, collapse = ", ")))
      message(paste("glstrinv.pval class:", class(nct_result$glstrinv.pval)))
      message(paste("glstrinv.pval value:", nct_result$glstrinv.pval))
      if (!is.null(nct_result$einv.pvals)) {
        message("Edge differences (einv.pvals) found!")
      }
      if (!is.null(nct_result$diffcen.pval)) {
        message("Centrality differences (diffcen.pval) found!")
      }

      # Store results with parameters used
      nct_results(list(
        result = nct_result,
        group1 = group1_name,
        group2 = group2_name,
        n1 = nrow(data1),
        n2 = nrow(data2),
        vars = varnames,
        paired = is_paired,
        estimator = current_estimator,
        is_ising = is_ising,
        iterations = input$nct_iterations %||% 1000,
        gamma = nct_gamma,
        and_rule = nct_and
      ))

      message("=== NCT Complete ===\n")

      removeNotification("nct_progress")
      showNotification("Network Comparison Test completed!", type = "message", duration = 5)

    }, error = function(e) {
      removeNotification("nct_progress")
      message(paste("ERROR in NCT:", e$message))
      showNotification(paste("Error running NCT:", e$message),
                      type = "error", duration = 10)
    })
  })

  # NCT Summary output
  output$nct_summary <- renderPrint({
    req(nct_results())
    nct <- nct_results()

    cat("===================================\n")
    cat("NETWORK COMPARISON TEST RESULTS\n")
    cat("===================================\n\n")

    cat("Test Parameters:\n")
    cat(sprintf("  Estimator: %s\n", nct$estimator))
    if (!is.null(nct$is_ising) && nct$is_ising) {
      cat(sprintf("  Network type: Ising (binary data)\n"))
      cat(sprintf("  AND-rule: %s\n", nct$and_rule))
    } else {
      cat(sprintf("  Network type: Gaussian Graphical Model\n"))
    }
    cat(sprintf("  Comparison type: %s\n", if(nct$paired) "Paired" else "Independent (unpaired)"))
    cat(sprintf("  EBIC gamma: %.2f\n", nct$gamma))
    cat(sprintf("  Permutations: %d\n\n", nct$iterations))

    cat("Groups Compared:\n")
    cat(sprintf("  Group 1: %s (n = %d)\n", nct$group1, nct$n1))
    cat(sprintf("  Group 2: %s (n = %d)\n", nct$group2, nct$n2))
    cat(sprintf("  Variables: %d\n\n", length(nct$vars)))

    # Network structure invariance test (uses nwinv.pval, not nw.pval)
    cat("-----------------------------------\n")
    cat("NETWORK STRUCTURE INVARIANCE TEST\n")
    cat("-----------------------------------\n")

    if (!is.null(nct$result$nwinv.pval)) {
      # Display maximum difference
      if (!is.null(nct$result$nwinv.real)) {
        cat(sprintf("Maximum edge difference: %.4f\n", nct$result$nwinv.real))
      }

      # Get p-value
      pval <- if (is.numeric(nct$result$nwinv.pval) && length(nct$result$nwinv.pval) == 1) {
        nct$result$nwinv.pval
      } else if (is.list(nct$result$nwinv.pval)) {
        nct$result$nwinv.pval[["pval"]] %||% nct$result$nwinv.pval[[1]]
      } else {
        NA
      }

      if (!is.na(pval)) {
        cat(sprintf("p-value: %.4f\n\n", pval))
        if (pval < 0.05) {
          cat("*** Network structures are significantly different (p < 0.05) ***\n\n")
        } else {
          cat("Network structures are NOT significantly different (p >= 0.05)\n\n")
        }
      } else {
        cat("Unable to extract p-value\n\n")
      }
    } else {
      cat("ℹ️  Network structure invariance test requires 'Test individual edges' option\n\n")
    }
  })

  # NCT Global Strength Table
  output$nct_global_table <- renderTable({
    req(nct_results())
    nct <- nct_results()

    # Check if global strength test was performed
    if (is.null(nct$result$glstrinv.pval)) {
      return(data.frame(Message = "Global strength test not performed"))
    }

    tryCatch({
      # NCT returns glstrinv.real as the DIFFERENCE in global strength (single value)
      # NOT the individual group values
      # glstrinv.pval is the p-value for this difference

      diff_val <- if (!is.null(nct$result$glstrinv.real) && length(nct$result$glstrinv.real) == 1)
        nct$result$glstrinv.real[1] else NA

      pval <- if (is.numeric(nct$result$glstrinv.pval) && length(nct$result$glstrinv.pval) == 1) {
        nct$result$glstrinv.pval
      } else if (is.list(nct$result$glstrinv.pval)) {
        nct$result$glstrinv.pval[["pval"]]
      } else {
        NA
      }

      # Interpretation
      sig_text <- if (!is.na(pval)) {
        if (pval < 0.05) "Significantly different" else "Not significantly different"
      } else {
        "N/A"
      }

      data.frame(
        Measure = c("Global Strength Difference", "p-value", "Interpretation"),
        Value = c(
          if (!is.na(diff_val)) sprintf("%.4f", diff_val) else "N/A",
          if (!is.na(pval)) sprintf("%.4f", pval) else "N/A",
          sig_text
        ),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(Message = paste("Error displaying results:", e$message))
    })
  }, rownames = FALSE, colnames = TRUE, spacing = "m", width = "100%")

  # NCT Edge Differences Table
  output$nct_edges_table <- DT::renderDataTable({
    req(nct_results())
    req(input$nct_test_edges)

    nct <- nct_results()
    if (is.null(nct$result$einv.pvals)) return(NULL)

    # Get edge differences - structure may vary
    edges_df <- as.data.frame(nct$result$einv.pvals)

    if (nrow(edges_df) == 0) {
      return(DT::datatable(data.frame(Message = "No edge differences found")))
    }

    # Add significant column if p-value column exists
    pval_col <- which(tolower(colnames(edges_df)) == "p-value")
    if (length(pval_col) > 0) {
      edges_df$significant <- edges_df[[pval_col]] < 0.05
    }

    # Create datatable
    dt <- DT::datatable(
      edges_df,
      options = list(
        pageLength = 10,
        order = list(list(0, 'asc'))
      ),
      rownames = FALSE
    )

    # Add formatting if columns exist
    if ("significant" %in% colnames(edges_df)) {
      dt <- dt %>%
        DT::formatStyle(
          'significant',
          target = 'row',
          backgroundColor = DT::styleEqual(c(TRUE, FALSE), c('#ffcccc', 'white'))
        )
    }

    # Format numeric columns
    numeric_cols <- names(edges_df)[sapply(edges_df, is.numeric)]
    if (length(numeric_cols) > 0) {
      dt <- dt %>% DT::formatRound(columns = numeric_cols, digits = 4)
    }

    dt
  })

  # NCT Centrality Differences Table
  output$nct_centrality_table <- DT::renderDataTable({
    req(nct_results())
    req(input$nct_test_centrality)

    nct <- nct_results()
    if (is.null(nct$result$diffcen.pval)) return(NULL)

    tryCatch({
      # Handle different possible structures
      diff_real <- nct$result$diffcen.real
      diff_pval <- nct$result$diffcen.pval

      if (is.null(diff_real) || is.null(diff_pval)) {
        return(DT::datatable(data.frame(Message = "Centrality test results not available")))
      }

      # Convert to data frame
      if (is.matrix(diff_real)) {
        cent_df <- as.data.frame(diff_real)
        cent_df$Node <- rownames(diff_real)
      } else if (is.vector(diff_real)) {
        cent_df <- data.frame(
          Node = names(diff_real),
          Difference = as.numeric(diff_real),
          stringsAsFactors = FALSE
        )
      } else {
        return(DT::datatable(data.frame(Message = "Unexpected centrality data structure")))
      }

      # Add p-values
      if (is.matrix(diff_pval)) {
        pval_df <- as.data.frame(diff_pval)
        cent_df <- cbind(cent_df, pval_df)
      } else if (is.vector(diff_pval)) {
        cent_df$p_value <- as.numeric(diff_pval)
      }

      # Add significance column if p_value exists
      if ("p_value" %in% colnames(cent_df)) {
        cent_df$significant <- cent_df$p_value < 0.05
      }

      # Reorder columns to put Node first
      if ("Node" %in% colnames(cent_df)) {
        other_cols <- setdiff(colnames(cent_df), "Node")
        cent_df <- cent_df[, c("Node", other_cols)]
      }

      # Create datatable
      dt <- DT::datatable(
        cent_df,
        options = list(
          pageLength = 10,
          order = list(list(0, 'asc'))
        ),
        rownames = FALSE
      )

      # Add formatting
      if ("significant" %in% colnames(cent_df)) {
        dt <- dt %>%
          DT::formatStyle(
            'significant',
            target = 'row',
            backgroundColor = DT::styleEqual(c(TRUE, FALSE), c('#ffcccc', 'white'))
          )
      }

      # Format numeric columns
      numeric_cols <- names(cent_df)[sapply(cent_df, is.numeric)]
      if (length(numeric_cols) > 0) {
        dt <- dt %>% DT::formatRound(columns = numeric_cols, digits = 4)
      }

      dt
    }, error = function(e) {
      DT::datatable(data.frame(Message = paste("Error displaying centrality:", e$message)))
    })
  })

  # ── Bridge Symptoms ──────────────────────────────────────────────────────────

  # Custom bridge plot: bypasses networktools::plot.bridge which breaks on
  # newer ggplot2 versions (object 'Group' not found in geom_path aesthetics).
  plot_bridge_custom <- function(bridge_obj,
                                 include = "Bridge Expected Influence (1-step)",
                                 order   = "value",
                                 main    = "") {
    vals <- bridge_obj[[include]]
    if (is.null(vals) || length(vals) == 0) {
      message(sprintf("[FALLBACK] plot_bridge_custom: metric '%s' not found in bridge object", include))
      available <- names(bridge_obj)[sapply(bridge_obj, is.numeric)]
      message(sprintf("           Available metrics: %s", paste(available, collapse = ", ")))
      plot.new()
      title(paste("No data for:", include))
      return(invisible(NULL))
    }

    df <- data.frame(
      node  = names(vals) %||% as.character(seq_along(vals)),
      value = as.numeric(vals),
      stringsAsFactors = FALSE
    )
    if (order == "value") df <- df[order(df$value, decreasing = FALSE), ]
    df$node <- factor(df$node, levels = df$node)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$node)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$value,
                                         y = .data$node, yend = .data$node),
                            color = "grey70") +
      ggplot2::geom_point(size = 3, color = "steelblue") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(x = include, y = NULL, title = main) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
    print(p)
  }

  bridge_result <- eventReactive(input$run_bridge, {
    res <- est_result_thr()
    shiny::validate(shiny::need(!is.null(res), "Estimate a network first."))

    run_bridge_for_net <- function(net, comm_method, algo, manual_text) {
      W <- safe_get_W(net)
      if (is.null(W)) return(NULL)
      W <- as.matrix(W); storage.mode(W) <- "double"; W[!is.finite(W)] <- 0

      comm_vec <- if (comm_method == "manual" && nzchar(trimws(manual_text %||% ""))) {
        raw <- trimws(unlist(strsplit(manual_text, ",")))
        if (length(raw) != ncol(W)) {
          showNotification(
            paste0("Community vector length (", length(raw), ") \u2260 number of nodes (", ncol(W), ")"),
            type = "error", duration = 10)
          return(NULL)
        }
        raw
      } else {
        detect_bridge_communities(W, algo)
      }
      if (is.null(comm_vec)) return(NULL)

      bridge_res <- tryCatch(
        networktools::bridge(W, communities = comm_vec),
        error = function(e) {
          showNotification(paste("bridge() error:", e$message), type = "error", duration = 10)
          NULL
        }
      )
      list(bridge = bridge_res, communities = comm_vec, W = W)
    }

    if (res$split) {
      results <- list()
      for (g in names(res$groups)) {
        results[[g]] <- run_bridge_for_net(res$groups[[g]],
                                            input$bridge_comm_method,
                                            input$bridge_auto_algo,
                                            input$bridge_comm_manual)
      }
      list(split = TRUE, groups = results)
    } else {
      r <- run_bridge_for_net(res$network,
                               input$bridge_comm_method,
                               input$bridge_auto_algo,
                               input$bridge_comm_manual)
      list(split = FALSE, result = r)
    }
  }, ignoreInit = TRUE)

  output$bridge_comm_info <- renderPrint({
    br <- bridge_result()
    shiny::validate(shiny::need(!is.null(br), "Run bridge analysis first."))
    if (br$split) {
      for (g in names(br$groups)) {
        r <- br$groups[[g]]
        if (is.null(r)) next
        cv <- r$communities
        cat(sprintf("=== Group: %s ===\n", g))
        cat(sprintf("Communities detected: %d\n", length(unique(cv))))
        node_comm <- data.frame(Node = names(cv) %||% seq_along(cv), Community = cv)
        print(node_comm, row.names = FALSE)
        cat("\n")
      }
    } else {
      r <- br$result
      shiny::validate(shiny::need(!is.null(r), "Bridge analysis failed — check community settings."))
      cv <- r$communities
      cat(sprintf("Communities detected: %d\n", length(unique(cv))))
      node_comm <- data.frame(Node = names(cv) %||% seq_along(cv), Community = cv)
      print(node_comm, row.names = FALSE)
    }
  })

  output$bridge_plot <- renderPlot({
    br <- bridge_result()
    shiny::validate(shiny::need(!is.null(br), "Run bridge analysis first."))
    if (br$split) {
      grps  <- br$groups
      plots <- lapply(names(grps), function(g) {
        r <- grps[[g]]
        if (is.null(r) || is.null(r$bridge)) return(NULL)
        # Capture the ggplot object from the custom function
        vals <- r$bridge[["Bridge Expected Influence (1-step)"]]
        if (is.null(vals)) {
          message(sprintf("[FALLBACK] bridge plot group '%s': metric not found", g))
          return(NULL)
        }
        df <- data.frame(
          node  = names(vals) %||% as.character(seq_along(vals)),
          value = as.numeric(vals),
          stringsAsFactors = FALSE
        )
        df <- df[order(df$value), ]
        df$node <- factor(df$node, levels = df$node)
        ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$node)) +
          ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$value,
                                             y = .data$node, yend = .data$node),
                                color = "grey70") +
          ggplot2::geom_point(size = 3, color = "steelblue") +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
          ggplot2::labs(x = "Bridge Expected Influence (1-step)", y = NULL,
                        title = paste("Group:", g)) +
          ggplot2::theme_minimal(base_size = 12) +
          ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
      })
      plots <- Filter(Negate(is.null), plots)
      if (length(plots) > 0) {
        gridExtra::grid.arrange(grobs = plots, ncol = length(plots))
      }
    } else {
      r <- br$result
      shiny::validate(shiny::need(!is.null(r) && !is.null(r$bridge),
                                   "Bridge analysis failed — check community settings."))
      plot_bridge_custom(r$bridge, include = "Bridge Expected Influence (1-step)",
                         order = "value")
    }
  })

  output$bridge_table <- DT::renderDataTable({
    br <- bridge_result()
    shiny::validate(shiny::need(!is.null(br), "Run bridge analysis first."))

    extract_bridge_df <- function(r, group = NULL) {
      b <- r$bridge
      if (is.null(b)) return(NULL)
      metric_names <- c("Bridge Strength", "Bridge Expected Influence (1-step)",
                        "Bridge Expected Influence (2-step)",
                        "Bridge Betweenness", "Bridge Closeness")
      dfs <- lapply(metric_names, function(m) {
        vals <- b[[m]]
        if (is.null(vals)) return(NULL)
        data.frame(Node = names(vals) %||% seq_along(vals),
                   Metric = m,
                   Value = round(as.numeric(vals), 3),
                   stringsAsFactors = FALSE)
      })
      df <- do.call(rbind, Filter(Negate(is.null), dfs))
      if (!is.null(group)) df <- cbind(Group = group, df)
      df
    }

    if (br$split) {
      dfs <- lapply(names(br$groups), function(g) extract_bridge_df(br$groups[[g]], g))
      df <- do.call(rbind, Filter(Negate(is.null), dfs))
    } else {
      df <- extract_bridge_df(br$result)
    }
    shiny::validate(shiny::need(!is.null(df) && nrow(df) > 0, "No bridge metrics available."))
    DT::datatable(df, options = list(pageLength = 15), rownames = FALSE)
  })

  # ── Bridge Stability ─────────────────────────────────────────────────────────

  bridge_stability_result <- eventReactive(input$run_bridge_stability, {
    br <- bridge_result()
    shiny::validate(shiny::need(!is.null(br), "Run bridge analysis first."))
    res <- est_result_thr()

    run_stability <- function(net, comm_vec) {
      tryCatch(
        bootnet::bootnet(net, boots = input$bridge_boots, type = "case",
                         statistics = "bridgeExpectedInfluence",
                         communities = comm_vec,
                         nCores = input$bridge_ncores),
        error = function(e) {
          showNotification(paste("Stability error:", e$message), type = "error", duration = 15)
          NULL
        }
      )
    }

    if (!br$split) {
      boot_res <- run_stability(res$network, br$result$communities)
      list(split = FALSE, result = boot_res)
    } else {
      results <- list()
      for (g in names(res$groups)) {
        results[[g]] <- run_stability(res$groups[[g]], br$groups[[g]]$communities)
      }
      list(split = TRUE, groups = results)
    }
  }, ignoreInit = TRUE)

  output$bridge_stability_plot <- renderPlot({
    bs <- bridge_stability_result()
    shiny::validate(shiny::need(!is.null(bs), "Run bridge stability first."))
    if (bs$split) {
      n <- length(bs$groups)
      par(mfrow = c(1, n))
      for (g in names(bs$groups)) {
        b <- bs$groups[[g]]
        if (!is.null(b))
          plot(b, statistics = "bridgeExpectedInfluence", main = paste("Group:", g))
      }
    } else {
      shiny::validate(shiny::need(!is.null(bs$result), "Stability bootstrap failed."))
      plot(bs$result, statistics = "bridgeExpectedInfluence")
    }
  })

  output$bridge_cs_output <- renderPrint({
    bs <- bridge_stability_result()
    shiny::validate(shiny::need(!is.null(bs), "Run bridge stability first."))
    if (bs$split) {
      for (g in names(bs$groups)) {
        b <- bs$groups[[g]]
        if (!is.null(b)) {
          cat(sprintf("=== Group: %s ===\n", g))
          print(bootnet::corStability(b))
          cat("\n")
        }
      }
    } else {
      shiny::validate(shiny::need(!is.null(bs$result), "Stability bootstrap failed."))
      print(bootnet::corStability(bs$result))
    }
  })

  # ── NSON: Nested Specificity-Oriented Network ─────────────────────────────

  nson_result <- eventReactive(input$run_nson, {
    res <- est_result_thr()
    shiny::validate(shiny::need(!is.null(res),
      "Estimate a network first (Estimation & Visualization tab)."))
    df <- prepped_data()
    shiny::validate(shiny::need(!is.null(df), "No usable data."))

    data_type <- if (identical(input$estimator, "IsingFit")) "binary" else "continuous"
    z_thr     <- input$nson_z_threshold %||% 0
    tol       <- input$nson_tolerance   %||% 0
    miss_pol  <- input$missing_policy   %||% "listwise"

    # Re-use already-computed predictability from main module (no re-estimation needed)
    pred_all <- tryCatch(predictability_vec(), error = function(e) NULL)

    run_nson_for_net <- function(net, group_data, pred_for_group = NULL, raw_group_data = NULL) {
      W <- safe_get_W(net)
      if (is.null(W)) return(NULL)
      W       <- as.matrix(W); storage.mode(W) <- "double"; W[is.na(W)] <- 0
      var_nms <- colnames(W)
      gdf     <- group_data[, intersect(var_nms, colnames(group_data)), drop = FALSE]

      gnss_tbl <- tryCatch(
        compute_gnss(gdf, type = data_type, z_threshold = z_thr, missing = miss_pol),
        error = function(e) { message("GNSS error: ", e$message); NULL }
      )
      if (is.null(gnss_tbl)) return(NULL)

      # Mean levels — use raw (unstandardized) data when available;
      # prepped_data() may be z-scored, giving means near 0 for continuous variables
      mean_src <- if (!is.null(raw_group_data)) {
        raw_group_data[, intersect(var_nms, colnames(raw_group_data)), drop = FALSE]
      } else { gdf }
      mean_lv <- tryCatch(
        setNames(vapply(mean_src, function(x) mean(as.numeric(x), na.rm = TRUE), numeric(1)),
                 colnames(mean_src)),
        error = function(e) setNames(rep(NA_real_, length(var_nms)), var_nms)
      )
      gnss_tbl$mean_level <- round(mean_lv[gnss_tbl$variable], 4)

      # Log-odds (logit) of activation rate — informative for both binary and thresholded continuous
      p_clamp <- pmax(0.001, pmin(0.999, gnss_tbl$activation_rate))
      gnss_tbl$logit <- round(log(p_clamp / (1 - p_clamp)), 4)

      # Predictability from main module's already-computed values
      if (!is.null(pred_for_group) && length(pred_for_group) > 0) {
        pv <- pred_for_group[gnss_tbl$variable]
        pv[!is.finite(pv)] <- NA_real_
        gnss_tbl$predictability <- round(pv, 4)
      }

      edge_tbl <- orient_network_by_gnss(W, gnss_tbl, tolerance = tol)
      list(W = W, gnss = gnss_tbl, edge_table = edge_tbl, group_df = gdf)
    }

    if (res$split) {
      full_df <- dat(); split <- split_var()
      results <- list()
      for (g in names(res$groups)) {
        g_data <- if (!is.null(split) && split %in% colnames(full_df)) {
          df[which(full_df[[split]] == g), , drop = FALSE]
        } else { df }
        raw_g <- if (!is.null(split) && split %in% colnames(full_df)) {
          full_df[which(full_df[[split]] == g), , drop = FALSE]
        } else { full_df }
        pred_g <- if (is.list(pred_all)) pred_all[[as.character(g)]] else pred_all
        results[[g]] <- run_nson_for_net(res$groups[[g]], g_data, pred_g, raw_g)
      }
      list(split = TRUE, groups = results, data_type = data_type)
    } else {
      full_raw <- tryCatch(dat(), error = function(e) NULL)
      pred_g   <- if (!is.list(pred_all)) pred_all else NULL
      r <- run_nson_for_net(res$network, df, pred_g, full_raw)
      shiny::validate(shiny::need(!is.null(r),
        "NSON computation failed. Check console for details."))
      c(r, list(split = FALSE, data_type = data_type))
    }
  }, ignoreInit = TRUE)

  output$nson_gnss_table <- DT::renderDataTable({
    nson <- nson_result()
    shiny::validate(shiny::need(!is.null(nson), "Click 'Compute NSON' first."))
    if (nson$split) {
      tbls <- lapply(names(nson$groups), function(g) {
        t <- nson$groups[[g]]$gnss
        if (!is.null(t)) cbind(Group = g, t) else NULL
      })
      tbl <- do.call(rbind, Filter(Negate(is.null), tbls))
    } else {
      tbl <- nson$gnss
    }
    shiny::validate(shiny::need(!is.null(tbl) && nrow(tbl) > 0, "No GNSS data."))
    num_cols <- intersect(c("activation_rate", "gnss", "mean_level", "logit", "predictability"), colnames(tbl))
    DT::datatable(tbl, options = list(pageLength = 20), rownames = FALSE) %>%
      DT::formatRound(columns = num_cols, digits = 4)
  })

  # Read and sanitize shared appearance inputs for NSON plots
  nson_appearance <- reactive({
    nc  <- input$nson_node_color    %||% "#FFFFFF"
    bc  <- input$nson_border_color  %||% "#8F8F8F"
    if (!is_valid_color(nc))  nc  <- "#FFFFFF"
    if (!is_valid_color(bc))  bc  <- "#8F8F8F"
    vs  <- input$nson_vsize         %||% 6;   if (!is.finite(vs)  || vs  <= 0) vs  <- 6
    es  <- input$nson_esize         %||% 3;   if (!is.finite(es)  || es  <= 0) es  <- 3
    as_ <- input$nson_asize         %||% 0.5; if (!is.finite(as_) || as_ <= 0) as_ <- 0.5
    me  <- input$nson_min_edge      %||% 0;   if (!is.finite(me)  || me  <  0) me  <- 0
    sm  <- input$nson_edge_scale_min %||% 0;  if (!is.finite(sm)  || sm  <  0) sm  <- 0
    sx  <- input$nson_edge_scale_max %||% 1;  if (!is.finite(sx)  || sx  <= sm) sx <- sm + 1
    lc  <- input$nson_label_cex     %||% 4.5; if (!is.finite(lc)  || lc  <= 0) lc  <- 4.5
    th  <- input$nson_theme %||% "classic"
    list(
      node_color   = nc,
      border_color = bc,
      vsize        = vs,
      esize        = es,
      asize        = as_,
      min_edge     = me,
      scl_min      = sm,
      scl_max      = sx,
      label_cex    = lc,
      label_bold   = isTRUE(input$nson_label_bold),
      show_labels  = isTRUE(input$nson_show_labels),
      size_by_mean = isTRUE(input$nson_size_by_mean),
      size_range   = { sr <- input$nson_size_range; if (is.null(sr) || length(sr)!=2 || anyNA(sr)) c(4,10) else sr },
      rev_scaling  = isTRUE(input$nson_reverse_scaling),
      theme        = th
    )
  })

  nson_comm_colors <- function(n) {
    pal <- c("#72AFD3","#FF6B92","#FFCC00","#74D6B7","#FFB07C",
             "#B09FCA","#85C1E9","#F1948A","#A9DFBF","#F9E79F")
    if (n <= length(pal)) pal[seq_len(n)] else grDevices::rainbow(n)
  }

  nson_vsize_vec <- function(var_names, ap, means_vec, global_from_range = NULL) {
    if (ap$size_by_mean && !is.null(means_vec) && !is.list(means_vec)) {
      m   <- means_vec[var_names]; m[!is.finite(m)] <- NA
      idx <- which(!is.na(m))
      sv  <- rep(ap$vsize, length(var_names))
      if (length(idx) > 0) {
        if (!is.null(global_from_range) && diff(global_from_range) > 0) {
          mn  <- global_from_range[1]; mx <- global_from_range[2]
          sc  <- pmax(0, pmin(1, (m[idx] - mn) / (mx - mn)))
          if (ap$rev_scaling) sc <- 1 - sc
          rng <- ap$size_range
          sv[idx] <- sc * (rng[2] - rng[1]) + rng[1]
        } else {
          sv[idx] <- scale_to_range(m[idx], rng = ap$size_range, reverse = ap$rev_scaling)
        }
      }
      sv
    } else {
      ap$vsize
    }
  }

  nson_layout_mat <- reactive({
    nson <- nson_result()
    req(!is.null(nson))
    lt  <- input$nson_layout          %||% "spring"
    sep <- input$nson_ega_separation  %||% 3

    compute_lm <- function(W, gdf) {
      p       <- ncol(W)
      var_nms <- colnames(W)
      spring_fallback <- function() {
        lm2 <- qgraph::qgraph(W, DoNotPlot = TRUE, layout = "spring")$layout
        rownames(lm2) <- var_nms
        list(layout = lm2, communities = NULL)
      }
      tryCatch({
        if (lt == "circle") {
          angles <- seq(0, 2 * pi, length.out = p + 1)[-(p + 1)]
          lm2    <- cbind(cos(angles), sin(angles))
          rownames(lm2) <- var_nms
          list(layout = lm2, communities = NULL)
        } else if (lt == "ega") {
          res  <- compute_ega_layout(W, data = gdf, return_communities = TRUE, separation = sep)
          lm2  <- if (is.list(res) && is.matrix(res$layout)) res$layout
                  else qgraph::qgraph(W, DoNotPlot = TRUE, layout = "spring")$layout
          rownames(lm2) <- var_nms
          comm <- if (is.list(res)) res$communities else NULL
          # k-means fallback when EGA finds ≤1 community
          if (is.null(comm) || length(unique(comm)) <= 1) {
            k    <- if (p >= 6) 3L else 2L
            comm <- tryCatch({
              km <- kmeans(lm2, centers = min(k, p - 1), nstart = 20)
              setNames(km$cluster, var_nms)
            }, error = function(e) NULL)
          }
          list(layout = lm2, communities = comm)
        } else if (lt == "pca") {
          lm2 <- compute_pca_layout(W, data = gdf)
          if (is.character(lm2)) lm2 <- qgraph::qgraph(W, DoNotPlot = TRUE, layout = "spring")$layout
          rownames(lm2) <- var_nms
          # Color by k-means on the PCA coordinates — geometrically meaningful,
          # entirely independent of graph topology (unlike EGA/Louvain)
          comm <- tryCatch({
            k   <- if (p >= 6) 3L else 2L
            km  <- kmeans(lm2, centers = k, nstart = 20)
            setNames(km$cluster, var_nms)
          }, error = function(e) NULL)
          list(layout = lm2, communities = comm)
        } else if (lt == "main") {
          stored <- main_net_layout()
          lm2 <- if (!is.null(stored) && is.matrix(stored) && nrow(stored) == p) {
            rownames(stored) <- var_nms; stored
          } else {
            lm_f <- qgraph::qgraph(W, DoNotPlot = TRUE, layout = "spring")$layout
            rownames(lm_f) <- var_nms; lm_f
          }
          list(layout = lm2, communities = NULL)
        } else {
          lm2 <- qgraph::qgraph(W, DoNotPlot = TRUE, layout = "spring")$layout
          rownames(lm2) <- var_nms
          list(layout = lm2, communities = NULL)
        }
      }, error = function(e) spring_fallback())
    }

    if (nson$split) {
      valid_groups <- Filter(Negate(is.null), nson$groups)
      if (length(valid_groups) == 0)
        return(list(split = TRUE, layouts = list(), communities = list()))

      # Compute layout ONCE from the first group so all groups share the same
      # coordinate system (comparable node positions across groups)
      first_res   <- compute_lm(valid_groups[[1]]$W, valid_groups[[1]]$group_df)
      shared_lm   <- first_res$layout

      # Shared communities: computed ONCE from first group, applied to ALL groups
      # so node colors are consistent across groups for direct comparison.
      # first_res$communities already has EGA + k-means fallback applied (via compute_lm).
      shared_ega_comm <- if (lt == "ega")  first_res$communities else NULL
      shared_pca_comm <- if (lt == "pca" && !is.null(shared_lm)) {
        tryCatch({
          k  <- if (ncol(valid_groups[[1]]$W) >= 6) 3L else 2L
          km <- kmeans(shared_lm, centers = k, nstart = 20)
          setNames(km$cluster, rownames(shared_lm))
        }, error = function(e) NULL)
      } else NULL

      # Per-group: reuse shared layout and shared communities
      results <- lapply(nson$groups, function(r) {
        if (is.null(r)) return(list(layout = shared_lm, communities = NULL))
        comm <- if      (lt == "ega") shared_ega_comm
                else if (lt == "pca") shared_pca_comm
                else                  NULL
        list(layout = shared_lm, communities = comm)
      })

      list(split      = TRUE,
           layouts     = lapply(results, function(x) x$layout),
           communities = lapply(results, function(x) x$communities))
    } else {
      res <- compute_lm(nson$W, nson$group_df)
      list(split = FALSE, layout = res$layout, communities = res$communities)
    }
  })

  nson_plot_helper <- function(W, gnss_tbl, layout_mat, tolerance,
                                title, show_gnss_labels, ap, vsize_vec,
                                communities = NULL, orient_by = "gnss") {
    var_names   <- colnames(W)
    p           <- ncol(W)
    gnss_vals   <- setNames(gnss_tbl$gnss, gnss_tbl$variable)
    act_rates   <- setNames(gnss_tbl$activation_rate, gnss_tbl$variable)
    orient_vals <- if (orient_by == "activation_rate") act_rates else gnss_vals
    min_edge    <- ap$min_edge

    ef <- character(0); et <- character(0)
    ew_abs <- numeric(0); e_sign <- character(0); ea <- logical(0)

    for (i in seq_len(p - 1)) {
      for (j in (i + 1):p) {
        w <- W[i, j]
        if (is.na(w) || abs(w) < 1e-10 || abs(w) < min_edge) next
        vi <- var_names[i]; vj <- var_names[j]
        gi <- orient_vals[vi]; gj <- orient_vals[vj]
        sgn <- ifelse(w > 0, "positive", "negative")
        aw  <- abs(w)
        if (is.na(gi) || is.na(gj) || abs(gi - gj) <= tolerance) {
          ef <- c(ef, vi); et <- c(et, vj)
          ew_abs <- c(ew_abs, aw); e_sign <- c(e_sign, sgn); ea <- c(ea, FALSE)
        } else if (gi < gj) {
          ef <- c(ef, vi); et <- c(et, vj)
          ew_abs <- c(ew_abs, aw); e_sign <- c(e_sign, sgn); ea <- c(ea, TRUE)
        } else {
          ef <- c(ef, vj); et <- c(et, vi)
          ew_abs <- c(ew_abs, aw); e_sign <- c(e_sign, sgn); ea <- c(ea, TRUE)
        }
      }
    }

    if (length(ef) == 0) {
      return(ggplot2::ggplot() +
               ggplot2::annotate("text", x=0.5, y=0.5, label="No edges to display", size=5) +
               ggplot2::theme_void() + ggplot2::ggtitle(title))
    }

    scl_min <- ap$scl_min; scl_max <- ap$scl_max
    ew_scaled <- (ew_abs - scl_min) / (scl_max - scl_min)
    ew_scaled[ew_scaled < 0] <- 0; ew_scaled[ew_scaled > 1] <- 1
    ewidths <- ew_scaled * ap$esize + 0.2

    # Alpha proportional to width: thin edges = transparent, thick edges = opaque
    ew_min <- min(ewidths); ew_max <- max(ewidths)
    ealpha <- if (ew_max > ew_min) {
      (ewidths - ew_min) / (ew_max - ew_min) * 0.65 + 0.35
    } else { rep(0.85, length(ewidths)) }

    theme_cols <- switch(ap$theme,
      colorblind = c(positive = "#009E73", negative = "#D55E00"),
      gray       = c(positive = "#888888", negative = "#333333"),
      Borkulo    = c(positive = "#009900", negative = "#FF0000"),
      Smurf      = c(positive = "#4A7EC5", negative = "#283269"),
                    c(positive = "#009900", negative = "#CC0000")
    )

    edf <- data.frame(
      from        = ef,
      to          = et,
      e_sign      = e_sign,
      ewidth      = ewidths,
      ealpha      = ealpha,
      is_directed = ea,
      stringsAsFactors = FALSE
    )
    vdf <- data.frame(name = var_names, stringsAsFactors = FALSE)
    g   <- igraph::graph_from_data_frame(edf, directed = TRUE, vertices = vdf)

    vorder <- igraph::V(g)$name

    gg_scale   <- 3.5   # ggraph size units are smaller than qgraph vsize units
    node_sizes <- if (length(vsize_vec) > 1) {
      sv <- vsize_vec[match(vorder, var_names)]
      sv[!is.finite(sv)] <- ap$vsize
      sv * gg_scale
    } else {
      rep(ap$vsize * gg_scale, length(vorder))
    }
    igraph::V(g)$node_size <- node_sizes

    # Node fill colors — palette when communities provided (NULL = user color)
    if (!is.null(communities) && length(communities) >= p) {
      comm_ids  <- communities[var_names]
      comm_ids[is.na(comm_ids)] <- max(comm_ids, na.rm = TRUE) + 1L
      u_comms   <- sort(unique(comm_ids))
      pal       <- nson_comm_colors(length(u_comms))
      col_map   <- setNames(pal, as.character(u_comms))
      fill_vec  <- unname(col_map[as.character(comm_ids)])
      igraph::V(g)$node_fill <- fill_vec[match(vorder, var_names)]
    } else {
      igraph::V(g)$node_fill <- rep(ap$node_color, length(vorder))
    }

    vlabels <- if (!ap$show_labels && !show_gnss_labels) {
      rep(NA_character_, length(vorder))
    } else if (isTRUE(show_gnss_labels)) {
      gv <- gnss_vals[vorder]
      paste0(vorder, "\n(", ifelse(is.na(gv), "NA", round(gv, 3)), ")")
    } else {
      vorder
    }
    igraph::V(g)$vlabel <- vlabels

    lmat_ord <- layout_mat[match(vorder, var_names), , drop = FALSE]
    lay      <- ggraph::create_layout(g, layout = "manual",
                                      x = lmat_ord[, 1], y = lmat_ord[, 2])

    label_size_gg <- ap$label_cex / 4.5 * 3.5
    if (!is.finite(label_size_gg) || label_size_gg <= 0) label_size_gg <- 3.5
    label_font    <- if (ap$label_bold) "bold" else "plain"
    cap_mm        <- mean(node_sizes, na.rm = TRUE) * 0.5

    ggraph::ggraph(lay) +
      ggraph::geom_edge_link(
        ggplot2::aes(colour = e_sign, width = ewidth, alpha = ealpha, filter = !is_directed),
        show.legend = FALSE
      ) +
      ggraph::geom_edge_link(
        ggplot2::aes(colour = e_sign, width = ewidth, alpha = ealpha, filter = is_directed),
        arrow     = grid::arrow(type = "closed", angle = 20,
                                length = grid::unit(ap$asize * 0.5, "cm")),
        end_cap   = ggraph::circle(cap_mm, "mm"),
        show.legend = FALSE
      ) +
      ggraph::scale_edge_colour_manual(values = theme_cols) +
      ggraph::scale_edge_width_identity() +
      ggraph::scale_edge_alpha_identity() +
      ggraph::geom_node_point(
        ggplot2::aes(size = node_size, fill = node_fill),
        colour = ap$border_color,
        shape  = 21,
        stroke = 0.8
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_fill_identity() +
      ggraph::geom_node_text(
        ggplot2::aes(label = vlabel),
        size     = label_size_gg,
        fontface = label_font,
        family   = "sans"
      ) +
      ggplot2::scale_x_continuous(expand = c(0.2, 0)) +
      ggplot2::scale_y_continuous(expand = c(0.2, 0)) +
      ggplot2::ggtitle(title) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 12))
  }

  output$nson_classic_plot_ui <- renderUI({
    h <- round(600 * (input$nson_zoom %||% 100) / 100)
    plotOutput("nson_classic_plot", height = paste0(h, "px"))
  })
  output$nson_directed_plot_ui <- renderUI({
    h <- round(600 * (input$nson_zoom %||% 100) / 100)
    plotOutput("nson_directed_plot", height = paste0(h, "px"))
  })

  output$nson_classic_plot <- renderPlot({
    nson <- nson_result()
    shiny::validate(shiny::need(!is.null(nson), "Click 'Compute NSON' first."))
    ap          <- nson_appearance()
    means_all <- tryCatch(means_for_nodes(), error = function(e) NULL)

    # Global mean range — only network variables, for cross-group comparability
    net_vars_c <- if (nson$split) {
      vg <- Filter(Negate(is.null), nson$groups)
      if (length(vg) > 0) colnames(vg[[1]]$W) else character(0)
    } else colnames(nson$W)
    classic_global_range <- if (ap$size_by_mean && is.list(means_all)) {
      all_vals <- unlist(lapply(means_all, function(m) {
        v <- m[intersect(names(m), net_vars_c)]; v[is.finite(v)]
      }), use.names = FALSE)
      if (length(all_vals) > 1) c(min(all_vals), max(all_vals)) else NULL
    } else NULL

    qgraph_theme <- switch(ap$theme, colorblind = "colorblind", gray = "gray", "classic")
    smurf_cols   <- if (ap$theme == "Smurf") list(posCol = "#4A7EC5", negCol = "#283269") else list()

    node_fill_vec <- function(W, communities) {
      vn <- colnames(W)
      if (!is.null(communities) && length(communities) >= ncol(W)) {
        comm_ids <- communities[vn]
        comm_ids[is.na(comm_ids)] <- max(comm_ids, na.rm = TRUE) + 1L
        u_comms  <- sort(unique(comm_ids))
        pal      <- nson_comm_colors(length(u_comms))
        col_map  <- setNames(pal, as.character(u_comms))
        unname(col_map[as.character(comm_ids)])
      } else {
        ap$node_color
      }
    }

    draw_classic <- function(W, layout_mat, title_str, communities, grp_means, gfr) {
      vn     <- colnames(W)
      W_viz  <- W
      if (ap$min_edge > 0) W_viz[abs(W_viz) < ap$min_edge] <- 0
      if (ap$scl_max > ap$scl_min)
        W_viz  <- rescale_edge_weights(W_viz, ap$scl_min, ap$scl_max)
      vsz    <- nson_vsize_vec(vn, ap, grp_means, gfr)
      lbls   <- if (!ap$show_labels) FALSE else vn
      nc     <- node_fill_vec(W, communities)
      do.call(suppressWarnings, list(qgraph::qgraph(
        W_viz,
        layout            = layout_mat %||% "spring",
        directed          = FALSE,
        theme             = qgraph_theme,
        title             = title_str,
        labels            = lbls,
        label.cex         = ap$label_cex * 0.7,
        label.font        = if (ap$label_bold) 2 else 1,
        label.scale.equal = TRUE,
        color             = nc,
        border.color      = ap$border_color,
        vsize             = vsz,
        esize             = ap$esize,
        minimum           = ap$min_edge,
        fade              = TRUE,
        posCol            = smurf_cols$posCol %||% "#009900",
        negCol            = smurf_cols$negCol %||% "#CC0000"
      )))
    }

    lm_data <- nson_layout_mat()
    if (nson$split) {
      grps <- nson$groups; n <- length(grps)
      par(mfrow = c(1, min(n, 3)))
      for (g in names(grps)) {
        r    <- grps[[g]]; if (is.null(r)) next
        lm   <- lm_data$layouts[[g]]
        comm <- lm_data$communities[[g]]
        grp_means <- if (is.list(means_all)) means_all[[as.character(g)]] else means_all
        draw_classic(r$W, lm, paste("Classic –", g), comm, grp_means, classic_global_range)
      }
    } else {
      means_vec <- if (!is.list(means_all)) means_all else NULL
      draw_classic(nson$W, lm_data$layout, "Classic Undirected Network",
                   lm_data$communities, means_vec, NULL)
    }
  })

  output$nson_directed_plot <- renderPlot({
    plot_mode <- input$nson_plot_mode %||% "gnss"
    lm_data   <- nson_layout_mat()
    ap        <- nson_appearance()

    if (plot_mode == "cp") {
      cp_res    <- cp_result()
      min_p     <- input$cp_min_prob      %||% 0.1
      slbl      <- isTRUE(input$cp_edge_labels)
      curv      <- input$cp_curvature     %||% 0.35
      reverse   <- isTRUE(input$cp_reverse_arrows)
      mode      <- cp_res$mode
      means_all <- tryCatch(means_for_nodes(), error = function(e) NULL)

      cp_global_range <- if (ap$size_by_mean && is.list(means_all)) {
        net_vars <- if (cp_res$split) {
          vg <- Filter(Negate(is.null), cp_res$groups)
          if (length(vg) > 0) rownames(vg[[1]]$cp) else character(0)
        } else rownames(cp_res$cp)
        all_vals <- unlist(lapply(means_all, function(m) {
          v <- m[intersect(names(m), net_vars)]; v[is.finite(v)]
        }), use.names = FALSE)
        if (length(all_vals) > 1) c(min(all_vals), max(all_vals)) else NULL
      } else NULL

      if (cp_res$split) {
        grps  <- cp_res$groups; n <- length(grps)
        plots <- lapply(grps, function(gr) {
          lm        <- lm_data$layouts[[gr$label]] %||% lm_data$layout
          shiny::validate(shiny::need(!is.null(lm), "No layout."))
          grp_means <- if (is.list(means_all)) means_all[[as.character(gr$label)]] else means_all
          vs        <- nson_vsize_vec(rownames(gr$cp), ap, grp_means, cp_global_range)
          cp_draw_graph(gr$cp, gr$W, lm, mode, min_p, slbl, ap,
                        paste("Conditional Probability —", gr$label), curv, reverse, vs)
        })
        gridExtra::grid.arrange(grobs = Filter(Negate(is.null), plots),
                                ncol = min(n, 2))
      } else {
        lm        <- lm_data$layout
        shiny::validate(shiny::need(!is.null(lm), "No layout."))
        means_vec <- if (!is.list(means_all)) means_all else NULL
        vs        <- nson_vsize_vec(rownames(cp_res$cp), ap, means_vec)
        print(cp_draw_graph(cp_res$cp, cp_res$W, lm, mode, min_p, slbl, ap,
                            "Conditional Probability Graph", curv, reverse, vs))
      }
    } else {
      nson <- nson_result()
      shiny::validate(shiny::need(!is.null(nson), "Click 'Compute NSON' first."))
      tol       <- input$nson_tolerance  %||% 0
      slbl      <- isTRUE(input$nson_show_gnss_labels)
      orient_by <- input$nson_orient_by  %||% "gnss"
      means_all <- tryCatch(means_for_nodes(), error = function(e) NULL)

      net_vars_d <- if (nson$split) {
        vg <- Filter(Negate(is.null), nson$groups)
        if (length(vg) > 0) colnames(vg[[1]]$W) else character(0)
      } else colnames(nson$W)
      directed_global_range <- if (ap$size_by_mean && is.list(means_all)) {
        all_vals <- unlist(lapply(means_all, function(m) {
          v <- m[intersect(names(m), net_vars_d)]; v[is.finite(v)]
        }), use.names = FALSE)
        if (length(all_vals) > 1) c(min(all_vals), max(all_vals)) else NULL
      } else NULL

      if (nson$split) {
        grps  <- nson$groups; n <- length(grps)
        plots <- lapply(names(grps), function(g) {
          r    <- grps[[g]]; if (is.null(r)) return(NULL)
          lm   <- lm_data$layouts[[g]];   if (is.null(lm)) return(NULL)
          comm <- lm_data$communities[[g]]
          grp_means <- if (is.list(means_all)) means_all[[as.character(g)]] else means_all
          vs   <- nson_vsize_vec(colnames(r$W), ap, grp_means, directed_global_range)
          nson_plot_helper(r$W, r$gnss, lm, tol,
                           paste0("NSON – Group: ", g), slbl, ap, vs, comm, orient_by)
        })
        plots <- Filter(Negate(is.null), plots)
        gridExtra::grid.arrange(grobs = plots, ncol = min(n, 3))
      } else {
        lm <- lm_data$layout
        shiny::validate(shiny::need(!is.null(lm), "Layout computation failed."))
        means_vec <- if (!is.list(means_all)) means_all else NULL
        vs <- nson_vsize_vec(colnames(nson$W), ap, means_vec)
        print(nson_plot_helper(nson$W, nson$gnss, lm, tol,
                               "Nested Specificity-Oriented Network (NSON)", slbl, ap, vs,
                               lm_data$communities, orient_by))
      }
    }
  })

  output$nson_edge_table <- DT::renderDataTable({
    nson <- nson_result()
    shiny::validate(shiny::need(!is.null(nson), "Click 'Compute NSON' first."))
    if (nson$split) {
      tbls <- lapply(names(nson$groups), function(g) {
        t <- nson$groups[[g]]$edge_table
        if (!is.null(t) && nrow(t) > 0) cbind(Group = g, t) else NULL
      })
      tbl <- do.call(rbind, Filter(Negate(is.null), tbls))
    } else {
      tbl <- nson$edge_table
    }
    shiny::validate(shiny::need(!is.null(tbl) && nrow(tbl) > 0, "No edges in network."))
    num_cols <- intersect(c("weight", "gnss_from", "gnss_to", "abs_gnss_diff"),
                          colnames(tbl))
    DT::datatable(tbl, options = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE) %>%
      DT::formatRound(columns = num_cols, digits = 4)
  })

  output$download_nson_edges <- downloadHandler(
    filename = function() paste0("nson_edges_", Sys.Date(), ".csv"),
    content  = function(file) {
      nson <- nson_result()
      if (is.null(nson)) return(NULL)
      if (nson$split) {
        tbls <- lapply(names(nson$groups), function(g) {
          t <- nson$groups[[g]]$edge_table
          if (!is.null(t) && nrow(t) > 0) cbind(Group = g, t) else NULL
        })
        tbl <- do.call(rbind, Filter(Negate(is.null), tbls))
      } else {
        tbl <- nson$edge_table
      }
      if (!is.null(tbl)) write.csv(tbl, file, row.names = FALSE)
    }
  )

  output$download_nson_activation <- downloadHandler(
    filename = function() paste0("nson_activation_", Sys.Date(), ".csv"),
    content  = function(file) {
      nson <- nson_result()
      if (is.null(nson) || nson$split) return(NULL)
      df_raw  <- nson$group_df
      dtype   <- nson$data_type
      z_thr   <- input$nson_z_threshold %||% 0
      act_mat <- if (dtype == "continuous") {
        m <- as.data.frame((scale(df_raw) > z_thr) * 1.0)
        colnames(m) <- colnames(df_raw); m
      } else { df_raw }
      write.csv(act_mat, file, row.names = FALSE)
    }
  )

  # ── Conditional Probability Graph ────────────────────────────────────────

  cp_result <- eventReactive(input$run_cp_graph, {
    nson  <- nson_result()
    shiny::validate(shiny::need(!is.null(nson),
      "Compute NSON first."))
    dtype <- nson$data_type
    z_thr <- input$nson_z_threshold %||% 0
    miss  <- input$missing_policy   %||% "listwise"
    mode  <- input$cp_mode          %||% "sparse"

    run_for_group <- function(gdf, W, label) {
      cp <- tryCatch(
        compute_cond_prob_matrix(gdf, type = dtype, z_threshold = z_thr, missing = miss),
        error = function(e) { message("CP error [", label, "]: ", e$message); NULL })
      if (is.null(cp)) return(NULL)
      list(cp = cp, W = W, label = label)
    }

    if (nson$split) {
      parts <- lapply(names(nson$groups), function(g) {
        r <- nson$groups[[g]]; if (is.null(r)) return(NULL)
        run_for_group(r$group_df, r$W, as.character(g))
      })
      list(split = TRUE, groups = Filter(Negate(is.null), parts),
           mode = mode, dtype = dtype)
    } else {
      r <- run_for_group(nson$group_df, nson$W, "All")
      shiny::validate(shiny::need(!is.null(r), "CP computation failed."))
      list(split = FALSE, cp = r$cp, W = r$W, mode = mode, dtype = dtype)
    }
  }, ignoreInit = TRUE)

  cp_build_edges <- function(cp_mat, W, mode, min_prob, reverse = FALSE) {
    p         <- nrow(cp_mat)
    var_names <- rownames(cp_mat)
    rows      <- vector("list", p * (p - 1))
    ri        <- 0L
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        if (i == j) next
        # Default (conventional): X→Y width = P(Y|X) = cp_mat[i,j]
        # Reversed:               X→Y width = P(X|Y) = cp_mat[j,i]
        prob <- if (reverse) cp_mat[j, i] else cp_mat[i, j]
        if (is.na(prob) || prob < min_prob) next
        if (mode == "sparse" && !is.null(W)) {
          if (abs(W[i,j]) < 1e-10 && abs(W[j,i]) < 1e-10) next
        }
        ri <- ri + 1L
        rows[[ri]] <- data.frame(from=var_names[i], to=var_names[j],
                                 prob=round(prob,4), stringsAsFactors=FALSE)
      }
    }
    if (ri == 0L) return(data.frame(from=character(), to=character(), prob=numeric()))
    do.call(rbind, rows[seq_len(ri)])
  }

  cp_draw_graph <- function(cp_mat, W, layout_mat, mode, min_prob,
                             show_labels, ap, title_str, curvature = 0.35,
                             reverse = FALSE, vsize_vec = NULL) {
    edf <- cp_build_edges(cp_mat, W, mode, min_prob, reverse)
    shiny::validate(shiny::need(nrow(edf) > 0,
      "No edges above the probability threshold. Lower the threshold."))

    var_names <- rownames(cp_mat)
    vdf <- data.frame(name = var_names, stringsAsFactors = FALSE)
    g   <- igraph::graph_from_data_frame(edf, directed = TRUE, vertices = vdf)

    vorder <- igraph::V(g)$name

    gg_scale   <- 3.5
    node_sizes <- if (!is.null(vsize_vec) && length(vsize_vec) > 1) {
      sv <- vsize_vec[match(vorder, var_names)]
      sv[!is.finite(sv)] <- ap$vsize
      sv * gg_scale
    } else {
      rep(ap$vsize * gg_scale, length(vorder))
    }
    igraph::V(g)$node_size <- node_sizes

    lmat_ord <- layout_mat[match(vorder, var_names), , drop = FALSE]
    lay      <- ggraph::create_layout(g, layout = "manual",
                                      x = lmat_ord[,1], y = lmat_ord[,2])

    cap_r    <- mean(node_sizes, na.rm = TRUE) * 0.5
    edge_aes <- if (show_labels) {
      ggplot2::aes(width = prob, alpha = prob, colour = prob,
                   label = round(prob, 2))
    } else {
      ggplot2::aes(width = prob, alpha = prob, colour = prob)
    }
    edge_layer <- ggraph::geom_edge_arc(
      edge_aes,
      curvature    = curvature,
      fold         = FALSE,
      arrow        = grid::arrow(type = "closed", angle = 18,
                                 length = grid::unit(0.35, "cm")),
      start_cap    = ggraph::circle(cap_r, "mm"),
      end_cap      = ggraph::circle(cap_r, "mm"),
      angle_calc   = "along",
      label_size   = 2.6,
      label_colour = "gray20",
      show.legend  = TRUE
    )

    ggraph::ggraph(lay) +
      edge_layer +
      ggraph::scale_edge_width_continuous(range = c(0.3, 3.5), guide = "none") +
      ggraph::scale_edge_alpha_continuous(range = c(0.2, 0.95), guide = "none") +
      ggraph::scale_edge_colour_gradient(
        low  = "#AED6F1", high = "#1A5276",
        name = if (reverse) "P(X=1|Y=1)" else "P(Y=1|X=1)") +
      ggraph::geom_node_point(
        ggplot2::aes(size = node_size),
        fill   = ap$node_color,
        colour = ap$border_color,
        shape  = 21, stroke = 0.8) +
      ggplot2::scale_size_identity() +
      ggraph::geom_node_text(
        ggplot2::aes(label = name),
        size     = ap$label_cex / 4.5 * 3.5, family = "sans",
        fontface = if (ap$label_bold) "bold" else "plain") +
      ggplot2::scale_x_continuous(expand = c(0.22, 0)) +
      ggplot2::scale_y_continuous(expand = c(0.22, 0)) +
      ggplot2::ggtitle(title_str) +
      ggraph::theme_graph(base_family = "sans") +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust=0.5, face="bold", size=13),
        legend.position = "top"
      )
  }

  output$download_cp_matrix <- downloadHandler(
    filename = function() paste0("cond_prob_matrix_", Sys.Date(), ".csv"),
    content  = function(file) {
      cp_res <- tryCatch(cp_result(), error=function(e) NULL)
      if (is.null(cp_res)) return(NULL)
      mat <- if (cp_res$split && length(cp_res$groups) > 0) cp_res$groups[[1]]$cp
             else cp_res$cp
      if (!is.null(mat)) {
        df_out <- as.data.frame(mat)
        df_out <- cbind(from=rownames(df_out), df_out)
        write.csv(df_out, file, row.names=FALSE)
      }
    }
  )

  # ── k-way Syndromic Depth & Centrality ───────────────────────────────────

  kway_result <- eventReactive(input$run_kway, {
    nson  <- nson_result()
    shiny::validate(shiny::need(!is.null(nson),
      "Compute NSON first (Estimation & Visualization tab → NSON tab → Compute NSON)."))
    dtype <- nson$data_type
    z_thr <- input$nson_z_threshold %||% 0
    miss  <- input$missing_policy   %||% "listwise"
    max_k <- as.integer(input$kway_max_k %||% 4L)

    run_for_data <- function(gdf, group_label = NULL) {
      res <- tryCatch(
        compute_kway_gnss(gdf, type = dtype, z_threshold = z_thr,
                          max_k = max_k, missing = miss),
        error = function(e) { message("k-way error [", group_label, "]: ", e$message); NULL }
      )
      if (!is.null(res) && !is.null(group_label)) res$group <- group_label
      res
    }

    if (nson$split) {
      parts <- lapply(names(nson$groups), function(g) {
        r <- nson$groups[[g]]; if (is.null(r)) return(NULL)
        run_for_data(r$group_df, group_label = as.character(g))
      })
      out <- do.call(rbind, Filter(Negate(is.null), parts))
    } else {
      out <- run_for_data(nson$group_df)
      if (!is.null(out)) out$group <- "All"
    }
    shiny::validate(shiny::need(!is.null(out) && nrow(out) > 0,
      "k-way computation returned no results. Check network size and max k."))
    out
  }, ignoreInit = TRUE)

  output$kway_max_k_ui <- renderUI({
    nson  <- tryCatch(nson_result(), error = function(e) NULL)
    n_nodes <- if (!is.null(nson)) {
      if (nson$split) {
        vg <- Filter(Negate(is.null), nson$groups)
        if (length(vg) > 0) ncol(vg[[1]]$W) else 10L
      } else ncol(nson$W)
    } else 10L

    abs_max <- max(2L, n_nodes - 1L)
    n_combos_warn <- if (abs_max >= 6) {
      HTML(sprintf(
        "<span style='color:#E67E22;font-size:11px;'>&#9888; Network has <strong>%d nodes</strong> (max k = %d). High k values compute C(%d, k−1) combinations per node — computation may take <strong>several minutes</strong> for k &gt; 5 with large networks.</span>",
        n_nodes, abs_max, n_nodes - 1L))
    } else NULL

    tagList(
      sliderInput("kway_max_k", "Maximum k",
                  min = 2, max = abs_max, value = min(4L, abs_max), step = 1),
      if (!is.null(n_combos_warn)) p(n_combos_warn) else
        helpText(sprintf("Max k = %d (number of nodes − 1).", abs_max))
    )
  })

  output$kway_group_selector_ui <- renderUI({
    nson <- tryCatch(nson_result(), error = function(e) NULL)
    if (is.null(nson) || !nson$split) return(NULL)
    grp_names <- names(nson$groups)
    shinyWidgets::pickerInput("kway_selected_group", "Display group",
                              choices = grp_names, selected = grp_names[1])
  })

  kway_plot_data <- reactive({
    df  <- kway_result()
    grp <- input$kway_selected_group %||% "All"
    if ("group" %in% colnames(df)) df <- df[df$group == grp, ]
    df
  })

  kway_node_order <- reactive({
    df   <- kway_plot_data()
    sort <- input$kway_sort_by %||% "k2_and"
    all_vars <- unique(df$variable)

    if (sort == "alpha") return(sort(all_vars))

    if (sort == "k2_and") {
      sub <- df[df$k == 2 & !is.na(df$and_gnss), c("variable","and_gnss")]
      sub <- sub[order(-sub$and_gnss), ]; return(unique(sub$variable))
    }
    if (sort == "k2_or") {
      sub <- df[df$k == 2 & !is.na(df$or_gnss), c("variable","or_gnss")]
      sub <- sub[order(-sub$or_gnss), ]; return(unique(sub$variable))
    }
    if (sort == "kmax_and") {
      mk  <- max(df$k, na.rm = TRUE)
      sub <- df[df$k == mk & !is.na(df$and_gnss), c("variable","and_gnss")]
      sub <- sub[order(-sub$and_gnss), ]; return(unique(sub$variable))
    }
    if (sort == "kmax_or") {
      mk  <- max(df$k, na.rm = TRUE)
      sub <- df[df$k == mk & !is.na(df$or_gnss), c("variable","or_gnss")]
      sub <- sub[order(-sub$or_gnss), ]; return(unique(sub$variable))
    }
    # fallback
    sub <- df[df$k == 2 & !is.na(df$and_gnss), c("variable","and_gnss")]
    sub <- sub[order(-sub$and_gnss), ]; unique(sub$variable)
  })

  make_kway_plot <- function(df, metric, style, node_order, text_sz = 9) {
    y_col   <- if (metric == "and") "and_gnss" else "or_gnss"
    k_range <- sort(unique(df$k))
    k_labs  <- paste0("k=", k_range)

    if (metric == "and") {
      pal <- setNames(colorRampPalette(c("#1A5276","#AED6F1"))(length(k_range)), k_labs)
      ttl <- if (style == "grouped") "AND Rule — Syndromic Decay"
             else                    "AND Rule — Total Syndromic Footprint"
      sub <- "P(all k−1 peers active | Z=1)"
    } else {
      pal <- setNames(colorRampPalette(c("#7B241C","#FDEBD0"))(length(k_range)), k_labs)
      ttl <- if (style == "grouped") "OR Rule — Vulnerability Spread"
             else                    "OR Rule — Total Comorbidity Footprint"
      sub <- "P(≥1 of k−1 peers active | Z=1)"
    }

    df_p          <- df[!is.na(df[[y_col]]), ]
    df_p$k_lab    <- factor(paste0("k=", df_p$k), levels = k_labs)
    df_p$variable <- factor(df_p$variable, levels = node_order)
    df_p$y_val    <- df_p[[y_col]]

    shared_theme <- ggplot2::theme_minimal(base_size = text_sz + 2) +
      ggplot2::theme(
        axis.text.x        = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5,
                                                   size = text_sz),
        axis.text.y        = ggplot2::element_text(size = text_sz),
        axis.title.y       = ggplot2::element_text(size = text_sz + 1),
        legend.text        = ggplot2::element_text(size = text_sz),
        legend.position    = "top",
        plot.title         = ggplot2::element_text(face = "bold", size = text_sz + 3),
        plot.subtitle      = ggplot2::element_text(size = text_sz - 1, color = "gray50"),
        plot.margin        = ggplot2::margin(5, 10, max(20, text_sz * 5), 5),
        panel.grid.major.x = ggplot2::element_blank()
      )

    y_scale <- ggplot2::scale_y_continuous(
      limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05)))

    if (style == "grouped") {
      # Line chart with filled bullet points per k level
      ggplot2::ggplot(df_p, ggplot2::aes(x = variable, y = y_val,
                                          color = k_lab, group = k_lab)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_point(size = 2.8, shape = 21,
                            ggplot2::aes(fill = k_lab), color = "white", stroke = 0.6) +
        ggplot2::scale_color_manual(values = pal, name = NULL) +
        ggplot2::scale_fill_manual(values  = pal, name = NULL) +
        ggplot2::labs(title = ttl, subtitle = sub, x = NULL, y = "Probability") +
        y_scale + shared_theme
    } else {
      # Stacked bar chart (unchanged)
      ggplot2::ggplot(df_p, ggplot2::aes(x = variable, y = y_val, fill = k_lab)) +
        ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.65) +
        ggplot2::scale_fill_manual(values = pal, name = NULL) +
        ggplot2::labs(title = ttl, subtitle = sub, x = NULL, y = "Probability") +
        y_scale + shared_theme
    }
  }

  output$kway_main_plot <- renderPlot({
    df    <- kway_plot_data()
    ord   <- kway_node_order()
    shiny::validate(shiny::need(nrow(df) > 0 && length(ord) > 0, "No data."))
    choice  <- input$kway_chart_choice %||% "and_grouped"
    parts   <- strsplit(choice, "_", fixed = TRUE)[[1]]
    txt_sz  <- input$kway_text_size %||% 9
    make_kway_plot(df, parts[1], parts[2], ord, text_sz = txt_sz)
  })

  output$download_kway <- downloadHandler(
    filename = function() paste0("kway_gnss_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- tryCatch(kway_result(), error = function(e) NULL)
      if (!is.null(df)) write.csv(df, file, row.names = FALSE)
    }
  )

}


shinyApp(ui, server)