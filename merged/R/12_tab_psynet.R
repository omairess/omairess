# ============================================================================
# 12_tab_psynet.R — psychonetrics tab + NCT comparison (from PsychoNetrix.R,
#                   NCT relocated here from BootSON per the merge plan)
#
# CORRECTNESS-CRITICAL LOGIC (do not alter casually):
#   * estimator RESOLVED to a concrete value before fitting — the string
#     "default" (psychonetrics' silent ML->FIML auto-switch) never reaches
#     the model call, the log, or the exported script
#   * ordinal-data guard: declaring ordered/dichotomous data with an
#     ML-family estimator triggers a visible warning (Rhemtulla et al. 2012);
#     when the data type is ordered/dichotomous the variables are passed via
#     `ordered =` so psychonetrics uses polychoric machinery
#   * NCT: paired = TRUE/FALSE is an explicit required choice; paired rows
#     are aligned by a subject-ID column (or equal-n + explicit row-order
#     assumption as fallback); edge-level p-values get multiple-comparison
#     correction (Holm default)
#
# DECISION: model families ported from PsychoNetrix.R:1365-1409 verbatim in
#   behaviour: lvm/lnm/rnm/lrnm via one builder switch (rnm & lrnm get
#   omega_epsilon = "full"), then runmodel -> optional prune(alpha, adjust)
#   -> stepup -> modelsearch, in that order.
# DECISION: the latent Relative-Importance network uses the EXACT Shapley/
#   LMG decomposition from PsychoNetrix.R:2135-2175 (ported verbatim), and
#   the exported fragment embeds the SAME function via deparse() so app and
#   script can never drift.
# DECISION: optimizer pinned to nlminb (recorded), not left to whatever
#   psychonetrics' current default happens to be.
# ============================================================================

# --- Estimator resolution: no "default" ever reaches runmodel ---------------
resolve_psynet_estimator <- function(choice, data_has_na) {
  if (choice != "auto") return(choice)
  if (data_has_na) "FIML" else "ML"
}

# --- Multi-group matrix collapse (PsychoNetrix.R:2056) ----------------------
collapse_mg <- function(m, grp = 1L) {
  if (is.list(m) && !is.data.frame(m)) return(as.matrix(m[[grp]]))
  if (is.array(m) && length(dim(m)) == 3) return(as.matrix(m[, , grp]))
  as.matrix(m)
}

# --- Latent Relative-Importance network (PsychoNetrix.R:2135, verbatim) -----
# Johnson/LMG relative weights via exact Shapley decomposition on the latent
# correlation matrix: RI[j, i] = average marginal R^2 contribution of latent
# j to outcome latent i over all predictor subsets. colSums(RI)[i] = total
# R^2 of outcome i on the other latents.
johnson_rw_from_cor <- function(R_lat) {
  p  <- ncol(R_lat)
  nm <- colnames(R_lat)
  if (is.null(nm)) nm <- paste0("L", seq_len(p))
  RI <- matrix(0, p, p, dimnames = list(nm, nm))

  r2_sub <- function(sub, rxy, Rxx) {
    if (length(sub) == 0L) return(0)
    rys <- rxy[sub]
    Rxs <- Rxx[sub, sub, drop = FALSE]
    tryCatch(max(0, drop(t(rys) %*% solve(Rxs, rys))),
             error = function(e) max(0, drop(t(rys) %*% (MASS::ginv(Rxs) %*% rys))))
  }

  for (i in seq_len(p)) {
    pred <- seq_len(p)[-i]
    q    <- length(pred)
    if (q < 1L) next
    rxy  <- R_lat[pred, i]
    Rxx  <- R_lat[pred, pred, drop = FALSE]
    if (q == 1L) { RI[pred, i] <- rxy^2; next }

    contrib <- numeric(q)
    for (j in seq_len(q)) {
      others <- seq_len(q)[-j]
      n_oth  <- length(others)
      for (mask in seq_len(2L^n_oth) - 1L) {
        S  <- if (n_oth == 0L) integer(0) else
                others[as.logical(intToBits(mask)[seq_len(n_oth)])]
        s  <- length(S)
        wt <- factorial(s) * factorial(q - 1L - s) / factorial(q)
        contrib[j] <- contrib[j] +
          wt * (r2_sub(c(S, j), rxy, Rxx) - r2_sub(S, rxy, Rxx))
      }
    }
    RI[pred, i] <- contrib
  }
  RI
}

# The exported script embeds the SAME function (deparse = zero drift).
JOHNSON_RW_FRAGMENT <- paste0(
  "johnson_rw_from_cor <- ",
  paste(deparse(johnson_rw_from_cor), collapse = "\n"))

# Bake a 0/1 lambda matrix as a runnable literal for the exported script.
lambda_literal <- function(lam) {
  sprintf(paste("lambda <- matrix(c(%s), nrow = %d,",
                "  dimnames = list(%s,",
                "                  %s))", sep = "\n"),
          paste(as.integer(lam), collapse = ", "), nrow(lam),
          vars_literal(rownames(lam)), vars_literal(colnames(lam)))
}

# --- Mean factor scores for latent node scaling (PsychoNetrix.R:2189-2212) --
# Regression-method factor scores, uncentered: colMeans(X %*% S^-1 %*% lambda
# %*% sigma_zeta) — the latent "severity level" each factor sits at. Latent
# node size scales by |this|, NOT by observed column means.
mean_fscores_psynet <- function(mod, lambda, dat) {
  tryCatch({
    sigma_z <- tryCatch(collapse_mg(psychonetrics::getmatrix(mod, "sigma_zeta"), 1L),
                        error = function(e) NULL)
    if (is.null(sigma_z)) sigma_z <- diag(ncol(lambda))
    keep <- intersect(rownames(lambda), colnames(dat))
    if (length(keep) != nrow(lambda)) stop("lambda/data name mismatch")
    dm <- as.matrix(dat[, keep, drop = FALSE])
    dm <- dm[stats::complete.cases(dm), , drop = FALSE]
    if (nrow(dm) < 2) stop("too few complete rows")
    S_inv <- tryCatch(solve(stats::cov(dm)),
                      error = function(e) MASS::ginv(stats::cov(dm)))
    colMeans(dm %*% S_inv %*% lambda %*% as.matrix(sigma_z))
  }, error = function(e) NULL)
}

# --- Combined factor-structure plot spec (from PsychoNetrix.R:2439-2588) ----
# Builds the (latents + indicators) super-network: latent-latent block =
# omega_zeta (undirected) or the Johnson/LMG RI network (directed), plus
# directed loading arrows latent -> indicator. Latents are circles on an
# inner ring; indicators are squares fanned around their dominant latent.
# Labels are always the FULL variable/factor names.
build_factor_spec <- function(mod, lambda_hint, type = c("latent", "ri"),
                              ri_norm = FALSE) {
  type    <- match.arg(type)
  lambda  <- tryCatch(collapse_mg(psychonetrics::getmatrix(mod, "lambda"), 1L),
                      error = function(e) NULL)
  if (is.null(lambda)) return(NULL)
  # psychonetrics matrices often come back with empty dimnames — restore the
  # names from the lambda the user built in the editor.
  if ((is.null(rownames(lambda)) || all(!nzchar(rownames(lambda)))) &&
      !is.null(lambda_hint) && nrow(lambda_hint) == nrow(lambda))
    rownames(lambda) <- rownames(lambda_hint)
  if ((is.null(colnames(lambda)) || all(!nzchar(colnames(lambda)))) &&
      !is.null(lambda_hint) && ncol(lambda_hint) == ncol(lambda))
    colnames(lambda) <- colnames(lambda_hint)
  lat_names <- colnames(lambda) %||% paste0("L", seq_len(ncol(lambda)))
  obs_labs  <- rownames(lambda) %||% paste0("V", seq_len(nrow(lambda)))

  n_lat <- ncol(lambda); n_obs <- nrow(lambda); n_total <- n_lat + n_obs

  lat_block <- if (type == "ri") {
    sigma_z <- tryCatch(collapse_mg(psychonetrics::getmatrix(mod, "sigma_zeta"), 1L),
                        error = function(e) NULL)
    if (is.null(sigma_z)) return(NULL)
    RI <- johnson_rw_from_cor(stats::cov2cor(sigma_z))
    if (ri_norm) { cs <- colSums(RI); cs[cs < 1e-10] <- 1; RI <- sweep(RI, 2, cs, "/") }
    RI
  } else {
    tryCatch(collapse_mg(psychonetrics::getmatrix(mod, "omega_zeta"), 1L),
             error = function(e) NULL)
  }
  if (is.null(lat_block)) return(NULL)

  W <- matrix(0, n_total, n_total)
  W[1:n_lat, 1:n_lat]             <- lat_block
  W[1:n_lat, (n_lat + 1):n_total] <- t(lambda)

  # loadings always directed; latent block directed only for the RI network
  dir_mat <- matrix(FALSE, n_total, n_total)
  dir_mat[1:n_lat, (n_lat + 1):n_total] <- TRUE
  if (type == "ri") dir_mat[1:n_lat, 1:n_lat] <- TRUE

  # boost near-zero loadings for DISPLAY so cross-loadings stay visible;
  # edge labels still show the true fitted values
  lambda_w <- lambda
  max_load <- if (any(lambda_w != 0)) max(abs(lambda_w[lambda_w != 0])) else 0
  if (max_load > 1e-6) {
    min_vis <- max_load * 0.20
    nz <- lambda_w != 0 & abs(lambda_w) < min_vis
    lambda_w[nz] <- sign(lambda_w[nz]) * min_vis
  }
  W_disp <- W
  W_disp[1:n_lat, (n_lat + 1):n_total] <- t(lambda_w)

  edge_lab <- matrix("", n_total, n_total)
  edge_lab[1:n_lat, 1:n_lat] <-
    ifelse(lat_block != 0, as.character(round(lat_block, 2)), "")
  edge_lab[1:n_lat, (n_lat + 1):n_total] <-
    ifelse(t(lambda) != 0, as.character(round(t(lambda), 2)), "")

  # latents on an inner ring; indicators fanned around their dominant latent
  dominant   <- apply(abs(lambda), 1, which.max)
  angles_lat <- seq(0, 2 * pi, length.out = n_lat + 1)[seq_len(n_lat)]
  lat_pos    <- cbind(cos(angles_lat), sin(angles_lat)) * 0.38
  obs_pos    <- matrix(0, n_obs, 2)
  for (k in seq_len(n_lat)) {
    idx <- which(dominant == k); n_k <- length(idx)
    if (n_k == 0) next
    half <- (pi / n_lat) * 0.85
    angs <- if (n_k == 1) angles_lat[k] else
              seq(angles_lat[k] - half, angles_lat[k] + half, length.out = n_k)
    obs_pos[idx, 1] <- cos(angs); obs_pos[idx, 2] <- sin(angs)
  }

  pal <- house_group_colors(n_lat)
  list(W = W_disp, dir_mat = dir_mat, edge_lab = edge_lab,
       labels = c(lat_names, obs_labs),
       shapes = c(rep("circle", n_lat), rep("square", n_obs)),
       node_cols = c(pal, pal[dominant]),
       layout = rbind(lat_pos, obs_pos),
       n_lat = n_lat,
       what = if (type == "ri")
                sprintf("factor structure + RI latent network (%s)",
                        if (ri_norm) "normalized" else "raw")
              else "factor structure + latent network (omega_zeta)")
}

# --- Fit-index interpretation (ported from PsychoNetrix.R:1535-1707) --------
# Thresholds + a plain-text "Model Fit Summary" that labels each index
# Excellent/Good/Acceptable/Poor with references. Behaviour unchanged.
PSYNET_FIT_THRESHOLDS <- list(
  CFI  = list(dir="high", cuts=c(.90,.95,.97), labs=c("Poor (<.90)","Acceptable (.90-.94)","Good (.95-.96)","Excellent (>=.97)")),
  TLI  = list(dir="high", cuts=c(.90,.95,.97), labs=c("Poor (<.90)","Acceptable (.90-.94)","Good (.95-.96)","Excellent (>=.97)")),
  NNFI = list(dir="high", cuts=c(.90,.95,.97), labs=c("Poor (<.90)","Acceptable (.90-.94)","Good (.95-.96)","Excellent (>=.97)")),
  NFI  = list(dir="high", cuts=c(.90,.95,.97), labs=c("Poor (<.90)","Acceptable (.90-.94)","Good (.95-.96)","Excellent (>=.97)")),
  IFI  = list(dir="high", cuts=c(.90,.95,.97), labs=c("Poor (<.90)","Acceptable (.90-.94)","Good (.95-.96)","Excellent (>=.97)")),
  RMSEA= list(dir="low",  cuts=c(.05,.08,.10), labs=c("Excellent (<=.05)","Good (.05-.08)","Acceptable (.08-.10)","Poor (>.10)")),
  SRMR = list(dir="low",  cuts=c(.05,.08,.10), labs=c("Excellent (<=.05)","Good (.05-.08)","Acceptable (.08-.10)","Poor (>.10)")),
  GFI  = list(dir="high", cuts=c(.85,.90,.95), labs=c("Poor (<.85)","Acceptable (.85-.89)","Good (.90-.94)","Excellent (>=.95)"))
)

print_psynet_fit_guide <- function(ft) {
  if (is.null(ft) || !is.data.frame(ft)) { print(ft); return(invisible(NULL)) }
  get_val <- function(nm) {
    idx <- which(toupper(ft$Measure) == toupper(nm))
    if (!length(idx)) return(NA_real_)
    suppressWarnings(as.numeric(ft$Value[idx[1]]))
  }
  classify <- function(val, spec) {
    if (is.na(val)) return(NA_character_)
    cu <- spec$cuts; la <- spec$labs
    if (spec$dir == "high")
      if (val >= cu[3]) la[4] else if (val >= cu[2]) la[3] else if (val >= cu[1]) la[2] else la[1]
    else
      if (val <= cu[1]) la[1] else if (val <= cu[2]) la[2] else if (val <= cu[3]) la[3] else la[4]
  }
  verdict <- function(val, nm) {
    spec <- PSYNET_FIT_THRESHOLDS[[toupper(nm)]]
    if (is.null(spec) || is.na(val)) return("")
    lab <- classify(val, spec); if (is.na(lab)) return("")
    tag <- if (grepl("^Excellent", lab)) "[++]" else if (grepl("^Good", lab)) "[+] "
           else if (grepl("^Acceptable", lab)) "[~] " else "[x] "
    paste0("  ", tag, " ", lab)
  }
  row1 <- function(nm, label = nm) {
    v <- get_val(nm); if (is.na(v)) return(invisible(NULL))
    cat(sprintf("  %-16s %7.3f%s\n", label, v, verdict(v, nm)))
  }
  SEP <- "  --------------------------------------------------------------\n"
  cat("==================================================================\n")
  cat("  Model Fit Summary\n")
  cat("==================================================================\n")
  info <- c(if (!is.na(get_val("nvar"))) sprintf("Variables: %g", get_val("nvar")),
            if (!is.na(get_val("npar"))) sprintf("Free parameters: %g", get_val("npar")),
            if (!is.na(get_val("df")))   sprintf("df: %g", get_val("df")))
  if (length(info)) cat("  ", paste(info, collapse = "   "), "\n", sep = "")
  chisq <- get_val("chisq"); pv <- get_val("pvalue")
  if (!is.na(chisq)) {
    cat("\n  Chi-square test of exact fit\n"); cat(SEP)
    cat(sprintf("  %-16s %7.2f%s\n", sprintf("chi2 (df=%g)", get_val("df")), chisq,
                if (is.na(pv)) "" else if (pv < .001) "  p < .001" else sprintf("  p = %.3f", pv)))
    cat("  (Significant chi2 is expected for large N; rely on indices below)\n")
  }
  if (any(!is.na(c(get_val("rmsea"), get_val("srmr"), get_val("gfi"))))) {
    cat("\n  Absolute fit\n"); cat(SEP)
    vr <- get_val("rmsea")
    if (!is.na(vr)) {
      lb <- get_val("rmsea.ci.lower"); ub <- get_val("rmsea.ci.upper")
      ci <- if (!is.na(lb) && !is.na(ub)) sprintf("  90%% CI [%.3f, %.3f]", lb, ub) else ""
      cat(sprintf("  %-16s %7.3f%s%s\n", "RMSEA", vr, ci, verdict(vr, "RMSEA")))
    }
    row1("srmr", "SRMR"); row1("gfi", "GFI")
  }
  if (any(vapply(c("cfi","tli","nnfi","nfi","ifi"), function(n) !is.na(get_val(n)), TRUE))) {
    cat("\n  Incremental fit\n"); cat(SEP)
    for (nm in c("cfi","tli","nnfi","nfi","ifi")) row1(nm, toupper(nm))
  }
  ic <- Filter(function(x) !is.na(x$v),
    lapply(list(c("aic","AIC"), c("bic","BIC"), c("ebic.5","eBIC(.5)")),
           function(p) list(lab = p[2], v = get_val(p[1]))))
  if (length(ic)) {
    cat("\n  Information criteria  (lower = better; for model comparison)\n"); cat(SEP)
    for (x in ic) cat(sprintf("  %-16s %10.2f\n", x$lab, x$v))
  }
  cat("\n"); cat(SEP)
  cat("  Legend: [++] Excellent  [+] Good  [~] Acceptable  [x] Poor\n")
  cat("  Refs: Hu & Bentler (1999); Browne & Cudeck (1992);\n")
  cat("        Schermelleh-Engel et al. (2003); Kline (2016).\n")
  cat("==================================================================\n")
}

psynetTabUI <- function(id) {
  ns <- shiny::NS(id)
  latent_family <- sprintf("['lvm','lnm','rnm','lrnm'].includes(input['%s'])",
                           ns("family"))
  shiny::tagList(
    varselectUI(ns("vars"), "Indicators (psychonetrics)"),
    shiny::selectInput(ns("data_type"), "Variable type",
                       c("Continuous" = "continuous",
                         "Ordered categorical" = "ordered",
                         "Dichotomous (0/1)" = "dichotomous")),
    shiny::selectInput(ns("estimator"), "Estimator",
                       c("Auto (ML; FIML if missing data)" = "auto",
                         "ML" = "ML", "FIML" = "FIML", "DWLS (ordinal)" = "DWLS",
                         "ULS" = "ULS", "WLS" = "WLS")),
    shiny::uiOutput(ns("ordinal_warning")),
    shiny::radioButtons(ns("family"), "Model family",
                        c("GGM (no latents)" = "ggm",
                          "CFA / SEM (lvm)" = "lvm",
                          "Latent Network Model (LNM)" = "lnm",
                          "Residual Network Model (RNM)" = "rnm",
                          "Latent + Residual Network (LRNM)" = "lrnm")),

    # -- Lambda (loading structure) editor, latent families only -------------
    shiny::conditionalPanel(latent_family,
      shiny::h5("Factor-loading structure (lambda)"),
      shiny::textInput(ns("latent_names"),
                       "Latent variable names (comma-separated)",
                       placeholder = "e.g. Prequel, Originals, Sequel"),
      shiny::actionButton(ns("gen_lambda"), "Generate lambda editor",
                          class = "btn-default btn-sm"),
      shiny::actionButton(ns("simple_structure"), "Auto simple structure",
                          class = "btn-default btn-sm"),
      shiny::helpText("Tick = free loading; untick = fixed to zero.",
                      "Each column is one latent factor."),
      shiny::uiOutput(ns("lambda_editor")),
      shiny::radioButtons(ns("identification"), "Identification",
                          c("Variance (fix latent variance)" = "variance",
                            "Loadings (fix first loading)"   = "loadings"),
                          selected = "variance")
    ),

    # -- Search & pruning (order: prune -> stepup -> modelsearch) ------------
    shiny::h5("Search & pruning"),
    shiny::checkboxInput(ns("do_prune"), "Prune non-significant edges", TRUE),
    shiny::conditionalPanel(sprintf("input['%s']", ns("do_prune")),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("prune_alpha"), "Prune alpha",
                                             value = 0.01, min = 0.001,
                                             max = 0.10, step = 0.005)),
        shiny::column(6, shiny::selectInput(ns("prune_adjust"), "Adjust p-values",
                                            c("none", "fdr", "bonferroni",
                                              "holm", "BY")))
      )),
    shiny::checkboxInput(ns("do_stepup"), "Step-up search (add edges)", FALSE),
    shiny::checkboxInput(ns("do_modelsearch"), "Full model search (slow)", FALSE),

    shiny::actionButton(ns("run"), "Fit model", class = "btn-primary"),

    # -- Results (fit + interpretation, parameters, MIs, matrices) -----------
    shiny::hr(),
    shiny::tabsetPanel(
      shiny::tabPanel("Fit indices",
        shiny::verbatimTextOutput(ns("fit_summary"))),
      shiny::tabPanel("Parameters",
        DT::dataTableOutput(ns("param_table"))),
      shiny::tabPanel("Modification indices",
        shiny::helpText("Largest suggested additions first. High MI = adding",
                        "that parameter would improve fit; corroborate with theory."),
        DT::dataTableOutput(ns("mi_table"))),
      shiny::tabPanel("Matrices",
        shiny::selectInput(ns("matrix_select"), "Matrix",
          c("omega", "omega_zeta", "sigma_zeta", "lambda",
            "omega_epsilon", "beta", "sigma", "kappa")),
        shiny::verbatimTextOutput(ns("matrix_output")))
    ),

    # -- Plots ----------------------------------------------------------------
    shiny::hr(),
    shiny::fluidRow(
      shiny::column(8, shiny::uiOutput(ns("psynet_plot_ui"))),
      shiny::column(4,
        shiny::selectInput(ns("plot_type"), "Network to plot",
          c("Factor structure + latent network"    = "factor_latent",
            "Factor structure + RI latent network" = "factor_ri",
            "Latent network (omega_zeta)"          = "latent",
            "RI latent network (raw)"              = "ri_raw",
            "RI latent network (normalized)"       = "ri_norm",
            "Residual network (omega_epsilon)"     = "residual",
            "Observed GGM (omega)"                 = "omega",
            "Path diagram (SEM, semPlot)"          = "semplot")),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'factor_ri'", ns("plot_type")),
          shiny::checkboxInput(ns("factor_ri_norm"),
                               "Normalize RI edges (column = share of R-squared)",
                               value = FALSE)),
        shiny::selectInput(ns("psynet_layout"), "Layout (matrix plots)",
                           c("Spring" = "spring", "Circle" = "circle"),
                           selected = "spring"),
        shiny::helpText("RI = Relative-Importance latent network",
                        "(Johnson/LMG): a DIRECTED network where the arrow",
                        "j -> i shows latent j's share of the R-squared of",
                        "latent i. Multi-group fits draw one panel per group",
                        "on a SHARED layout (matrix plots); composite/SEM",
                        "views show group 1."),
        # edge labels ON by default so the RI network shows its weights
        # as numbers on the arrows rather than hiding them (request #4).
        appearanceUI(ns("look"), signed = TRUE, default_edge_labels = TRUE,
                     scale_label = "Scale node size by mean latent score (latents) / column mean (indicators)")
      )
    )

    # TODO(port): panel families (dlvm1 / panelgvar / ri_clpm) + wave
    # detection + beta-matrix editor (PsychoNetrix.R:334-399, 1112-1250,
    # 1411-1468); Ising family; factor scores; data-transform pipeline.
  )
}

psynetTabServer <- function(id, data_bus, rec) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "psynet")
    rv  <- shiny::reactiveValues(model = NULL, family = NULL,
                                 estimator_used = NULL, lambda_used = NULL,
                                 lambda_template = NULL, plot_fn = NULL)

    # --- Ordinal/ML mismatch guard ------------------------------------------
    output$ordinal_warning <- shiny::renderUI({
      ml_family <- input$estimator %in% c("auto", "ML", "FIML")
      if (input$data_type %in% c("ordered", "dichotomous") && ml_family) {
        shiny::div(class = "alert alert-warning",
          "Ordinal/dichotomous items with an ML-family estimator treat categories
           as continuous. With < ~5-7 categories this biases estimates
           (Rhemtulla et al., 2012) — consider DWLS or ULS. The variables WILL
           be passed as ordered so psychonetrics can use polychoric input.")
      }
    })

    # --- Lambda editor (ported pattern from PsychoNetrix.R:1017-1110) -------
    shiny::observeEvent(input$gen_lambda, {
      shiny::req(sel$vars(), input$latent_names)
      latents <- trimws(strsplit(input$latent_names, ",")[[1]])
      latents <- latents[latents != ""]
      if (!length(latents)) {
        shiny::showNotification("Enter at least one latent name.",
                                type = "warning"); return()
      }
      rv$lambda_template <- matrix(0L, length(sel$vars()), length(latents),
                                   dimnames = list(sel$vars(), latents))
    })

    shiny::observeEvent(input$simple_structure, {
      shiny::req(sel$vars(), input$latent_names)
      latents <- trimws(strsplit(input$latent_names, ",")[[1]])
      latents <- latents[latents != ""]
      shiny::req(length(latents) >= 1)
      vars <- sel$vars()
      lam  <- matrix(0L, length(vars), length(latents),
                     dimnames = list(vars, latents))
      per  <- ceiling(length(vars) / length(latents))
      for (j in seq_along(latents)) {
        rows <- ((j - 1) * per + 1):min(j * per, length(vars))
        lam[rows, j] <- 1L
      }
      rv$lambda_template <- lam
      shiny::showNotification("Simple structure applied.", type = "message")
    })

    output$lambda_editor <- shiny::renderUI({
      shiny::req(rv$lambda_template)
      ns   <- session$ns
      lam  <- rv$lambda_template
      vars <- rownames(lam); latents <- colnames(lam)
      header <- shiny::tags$tr(
        shiny::tags$th("Variable"),
        lapply(latents, function(l)
          shiny::tags$th(l, style = "text-align:center; padding:2px 8px;")))
      rows <- lapply(seq_along(vars), function(i) shiny::tags$tr(
        shiny::tags$td(shiny::tags$b(vars[i]),
                       style = "padding:1px 8px; white-space:nowrap;"),
        lapply(seq_along(latents), function(j) shiny::tags$td(
          shiny::checkboxInput(ns(sprintf("lam_%d_%d", i, j)), NULL,
                               value = as.logical(lam[i, j]),
                               width = "24px"),
          style = "text-align:center; padding:0 8px; width:1%;"))))
      shiny::tags$div(style = "overflow-x:auto; display:inline-block;",
        # compact: kill checkbox form-group margins, shrink the table to
        # its content instead of the full panel width
        shiny::tags$style(shiny::HTML(
          ".lambda-tight table {width:auto !important; margin-bottom:4px;}
           .lambda-tight .checkbox {margin:0; min-height:0;}
           .lambda-tight .form-group {margin:0;}
           .lambda-tight .shiny-input-container {margin:0; padding:0; width:24px !important;}")),
        shiny::tags$div(class = "lambda-tight",
          shiny::tags$table(class = "table table-condensed table-bordered",
                            shiny::tags$thead(header),
                            shiny::tags$tbody(rows))))
    })

    get_lambda <- shiny::reactive({
      shiny::req(rv$lambda_template)
      lam <- rv$lambda_template
      for (i in seq_len(nrow(lam))) for (j in seq_len(ncol(lam))) {
        val <- input[[sprintf("lam_%d_%d", i, j)]]
        if (!is.null(val)) lam[i, j] <- if (isTRUE(val)) 1L else 0L
      }
      lam
    })

    # --- Fit (family dispatch ported from PsychoNetrix.R:1315-1409) ----------
    shiny::observeEvent(input$run, {
      vars <- sel$vars()
      dat  <- data_bus$wide()[, vars, drop = FALSE]
      # Per-module transform: after variable selection, before fitting.
      trans <- sel$transform()
      dat   <- apply_house_transform(dat, trans)
      fam  <- input$family
      est  <- resolve_psynet_estimator(input$estimator, anyNA(dat))
      ordered_vars <- if (input$data_type %in% c("ordered", "dichotomous"))
                        vars else character(0)
      grp  <- sel$group_var()
      if (!is.null(grp)) {
        gv  <- data_bus$wide()[[grp]]
        dat[[grp]] <- gv
        dat <- dat[!is.na(dat[[grp]]), , drop = FALSE]   # NA group crashes optimizer
      }

      lambda <- NULL
      if (fam %in% c("lvm", "lnm", "rnm", "lrnm")) {
        lambda <- get_lambda()
        shiny::validate(shiny::need(!is.null(lambda) && any(lambda == 1L),
          "Generate the lambda editor and free at least one loading first."))
        shiny::validate(shiny::need(setequal(rownames(lambda), vars),
          "Lambda rows no longer match the selected variables — regenerate the editor."))
      }

      mod <- tryCatch({
        shiny::withProgress(message = "Fitting psychonetrics model", value = 0.2, {
          m <- if (fam == "ggm") {
            psychonetrics::ggm(dat, vars = vars, omega = "full",
                               ordered = ordered_vars, estimator = est,
                               groupvar = grp)
          } else {
            builder <- switch(fam,
              lvm  = function(...) psychonetrics::lvm(...),
              lnm  = function(...) psychonetrics::lnm(...),
              rnm  = function(...) psychonetrics::rnm(..., omega_epsilon = "full"),
              lrnm = function(...) psychonetrics::lrnm(..., omega_epsilon = "full"))
            builder(data = dat, lambda = lambda, vars = rownames(lambda),
                    latents = colnames(lambda),
                    identification = input$identification,
                    ordered = ordered_vars, estimator = est, groupvar = grp)
          }
          # optimizer pinned; then prune -> stepup -> modelsearch (orig. order)
          m <- psychonetrics::runmodel(psychonetrics::setoptimizer(m, "nlminb"))
          if (isTRUE(input$do_prune))
            m <- psychonetrics::prune(m, alpha = input$prune_alpha,
                                      adjust = input$prune_adjust)
          if (isTRUE(input$do_stepup))      m <- psychonetrics::stepup(m)
          if (isTRUE(input$do_modelsearch)) m <- psychonetrics::modelsearch(m)
          m
        })
      }, error = function(e) {
        shiny::showNotification(paste("Model error:", conditionMessage(e)),
                                type = "error", duration = 10); NULL
      })
      shiny::req(mod)
      rv$model <- mod; rv$family <- fam; rv$dat <- dat
      rv$estimator_used <- est; rv$lambda_used <- lambda; rv$transform <- trans

      # -- recorder fragment: every setting explicit, lambda baked literal ---
      fit_line <- if (fam == "ggm") {
        sprintf(paste(
          'mod <- psychonetrics::ggm(dat_psynet, vars = psynet_vars,',
          '  omega = "full", ordered = %s, estimator = "%s"%s)', sep = "\n"),
          if (length(ordered_vars)) "psynet_vars" else "character(0)", est,
          if (!is.null(grp)) sprintf(', groupvar = "%s"', grp) else "")
      } else {
        sprintf(paste(
          "%s",
          'mod <- psychonetrics::%s(dat_psynet, lambda = lambda,',
          "  vars = rownames(lambda), latents = colnames(lambda),",
          '  identification = "%s", ordered = %s, estimator = "%s"%s%s)',
          sep = "\n"),
          lambda_literal(lambda), fam, input$identification,
          if (length(ordered_vars)) "psynet_vars" else "character(0)", est,
          if (fam %in% c("rnm", "lrnm")) ', omega_epsilon = "full"' else "",
          if (!is.null(grp)) sprintf(', groupvar = "%s"', grp) else "")
      }
      chain <- c(
        'mod <- psychonetrics::runmodel(psychonetrics::setoptimizer(mod, "nlminb"))',
        if (isTRUE(input$do_prune))
          sprintf('mod <- psychonetrics::prune(mod, alpha = %s, adjust = "%s")',
                  input$prune_alpha, input$prune_adjust),
        if (isTRUE(input$do_stepup))      "mod <- psychonetrics::stepup(mod)",
        if (isTRUE(input$do_modelsearch)) "mod <- psychonetrics::modelsearch(mod)")

      rec_upsert(
        rec, "psynet_analysis", "analysis",
        description = sprintf(
          "[psychonetrics] Fitted %s on %d variables (n = %d): estimator = %s (resolved explicitly%s), optimizer = nlminb, data type = %s, transform: %s%s%s%s.",
          toupper(fam), length(vars), nrow(dat), est,
          if (input$estimator == "auto") " from 'auto'" else "",
          input$data_type,
          names(TRANSFORM_LABELS)[TRANSFORM_LABELS == trans],
          if (fam != "ggm") sprintf(", identification = %s, %d latents",
                                    input$identification, ncol(lambda)) else "",
          if (isTRUE(input$do_prune))
            sprintf("; pruned at alpha = %s (adjust = %s)",
                    input$prune_alpha, input$prune_adjust) else "",
          if (isTRUE(input$do_stepup)) "; step-up search" else ""),
        code = paste(c(
          sprintf("psynet_vars <- %s", vars_literal(vars)),
          "dat_psynet  <- dat_wide[, psynet_vars]",
          transform_code_fragment(trans, "dat_psynet"),
          fit_line, chain), collapse = "\n")
      )
    })

    output$fit_summary <- shiny::renderPrint({
      shiny::req(rv$model)
      ft <- tryCatch(psychonetrics::fit(rv$model), error = function(e) NULL)
      print_psynet_fit_guide(ft)
    })

    output$param_table <- DT::renderDataTable({
      shiny::req(rv$model)
      pr <- tryCatch(rv$model@parameters, error = function(e) NULL)
      shiny::validate(shiny::need(!is.null(pr), "No parameters available."))
      keep <- intersect(c("var1", "op", "var2", "matrix", "row", "col",
                          "est", "se", "p", "std"), names(pr))
      tb <- as.data.frame(pr)[, keep, drop = FALSE]
      num <- intersect(c("est", "se", "p", "std"), names(tb))
      DT::datatable(tb, options = list(pageLength = 25, scrollX = TRUE)) |>
        DT::formatRound(num, 3)
    })

    output$mi_table <- DT::renderDataTable({
      shiny::req(rv$model)
      mi <- tryCatch(
        psychonetrics::MIs(rv$model, matrices = "omega", type = "free"),
        error = function(e) tryCatch(psychonetrics::MIs(rv$model),
                                     error = function(e2) NULL))
      shiny::validate(shiny::need(!is.null(mi) && NROW(mi) > 0,
        "No modification indices available for this model."))
      mi <- as.data.frame(mi)
      mic <- intersect(c("mi", "epc", "matrix", "row", "col"), names(mi))
      mi <- mi[, mic, drop = FALSE]
      if ("mi" %in% names(mi)) mi <- mi[order(-mi$mi), ]
      DT::datatable(mi, options = list(pageLength = 25, scrollX = TRUE)) |>
        DT::formatRound(intersect(c("mi", "epc"), names(mi)), 3)
    })

    output$matrix_output <- shiny::renderPrint({
      shiny::req(rv$model)
      m <- tryCatch(psychonetrics::getmatrix(rv$model, input$matrix_select),
                    error = function(e) NULL)
      if (is.null(m)) {
        cat(sprintf("Matrix '%s' is not available for this model.\n",
                    input$matrix_select))
      } else {
        print(round(collapse_mg(m, 1L), 4))
      }
    })

    # --- Plots: omega / latent / RI / residual -------------------------------
    look <- appearanceServer("look", plot_closure = shiny::reactive(rv$plot_fn))

    # Extract + prepare the matrix/matrices for the chosen plot type.
    # psychonetrics matrices often come back with EMPTY dimnames — without
    # restoring them qgraph would label nodes 1..p. Full names, always (#6).
    # Multi-group fits return W_list with one matrix PER GROUP (#9).
    plot_matrix <- shiny::reactive({
      shiny::req(rv$model)
      pt <- input$plot_type
      lat_names <- colnames(rv$lambda_used)
      obs_names <- if (!is.null(rv$lambda_used)) rownames(rv$lambda_used)
                   else tryCatch(sel$vars(), error = function(e) NULL)
      get_mats <- function(name, nms = NULL) {
        raw <- tryCatch(psychonetrics::getmatrix(rv$model, name),
                        error = function(e) NULL)
        if (is.null(raw)) return(NULL)
        mats <- if (is.list(raw) && !is.data.frame(raw))
                  lapply(raw, as.matrix) else list(as.matrix(raw))
        lapply(mats, function(m) {
          if ((is.null(colnames(m)) || all(!nzchar(colnames(m)))) &&
              !is.null(nms) && length(nms) == ncol(m))
            dimnames(m) <- list(nms, nms)
          m
        })
      }
      pack <- function(mats, directed, what) {
        if (is.null(mats)) return(list(W = NULL))
        list(W = mats[[1]], W_list = mats, directed = directed, what = what)
      }
      if (pt == "omega")
        return(pack(get_mats("omega", obs_names), FALSE, "observed GGM (omega)"))
      if (pt == "latent")
        return(pack(get_mats("omega_zeta", lat_names), FALSE,
                    "latent network (omega_zeta)"))
      if (pt == "residual")
        return(pack(get_mats("omega_epsilon", obs_names), FALSE,
                    "residual network (omega_epsilon)"))
      # RI networks from sigma_zeta -> latent correlations -> Johnson/LMG,
      # computed per group for multi-group fits
      szs <- get_mats("sigma_zeta", lat_names)
      if (is.null(szs)) return(list(W = NULL))
      ri_one <- function(sz) {
        R_lat <- tryCatch(stats::cov2cor(sz), error = function(e) NULL)
        if (is.null(R_lat)) return(NULL)
        RI <- johnson_rw_from_cor(R_lat)
        if (pt == "ri_norm") {
          cs <- colSums(RI); cs[cs < 1e-10] <- 1
          RI <- sweep(RI, 2, cs, "/")
        }
        RI
      }
      ris <- lapply(szs, ri_one)
      if (any(vapply(ris, is.null, TRUE))) return(list(W = NULL))
      pack(ris, TRUE, sprintf("RI latent network (%s)",
                              if (pt == "ri_norm") "normalized" else "raw"))
    })

    # Resizable plot window: height follows the appearance slider.
    output$psynet_plot_ui <- shiny::renderUI({
      shiny::plotOutput(session$ns("psynet_plot"),
                        height = sprintf("%dpx", look()$plot_height))
    })

    output$psynet_plot <- shiny::renderPlot({
      s  <- look()
      pt <- input$plot_type

      # --- SEM path diagram (semPlot), ported from PsychoNetrix.R:3078-3117 --
      if (pt == "semplot") {
        shiny::req(rv$model)
        shiny::validate(shiny::need(requireNamespace("semPlot", quietly = TRUE),
          "Install the 'semPlot' package for path diagrams."))
        lam <- tryCatch(collapse_mg(psychonetrics::getmatrix(rv$model, "lambda"), 1L),
                        error = function(e) NULL)
        shiny::validate(shiny::need(!is.null(lam),
          "The SEM path diagram needs a fitted CFA/LNM/RNM/LRNM model (lambda)."))
        psi <- tryCatch(collapse_mg(psychonetrics::getmatrix(rv$model, "sigma_zeta"), 1L),
                        error = function(e) diag(ncol(lam)))
        tht <- tryCatch(collapse_mg(psychonetrics::getmatrix(rv$model, "sigma_epsilon"), 1L),
                        error = function(e) diag(nrow(lam)))
        # psychonetrics matrices have empty dimnames — restore full names
        lat_l <- colnames(rv$lambda_used) %||% paste0("L", seq_len(ncol(lam)))
        obs_l <- rownames(rv$lambda_used) %||% paste0("V", seq_len(nrow(lam)))
        dimnames(lam) <- list(obs_l, lat_l)
        dimnames(psi) <- list(lat_l, lat_l)
        dimnames(tht) <- list(obs_l, obs_l)
        fn <- function() semPlot::semPaths(
          semPlot::lisrelModel(LY = lam, PS = psi, TE = tht),
          what = "std", whatLabels = "est", layout = "tree2",
          sizeLat = s$vsize + 2, sizeMan = s$vsize,
          edge.label.cex = s$edge_label_cex, mar = c(6, 1, 6, 1),
          nCharNodes = 0, nCharEdges = 0,   # do NOT abbreviate labels
          label.font = if (isTRUE(s$label_bold)) 2 else 1,
          pastel = TRUE, borders = TRUE)
        rv$plot_fn <- fn
        rec_upsert(
          rec, "psynet_plot", "plot",
          description = "[psychonetrics] Plotted the SEM path diagram (semPlot::semPaths, standardized, tree2 layout).",
          code = paste(
            'lam <- psychonetrics::getmatrix(mod, "lambda")',
            "if (is.list(lam)) lam <- lam[[1]]",
            'psi <- psychonetrics::getmatrix(mod, "sigma_zeta"); if (is.list(psi)) psi <- psi[[1]]',
            'tht <- psychonetrics::getmatrix(mod, "sigma_epsilon"); if (is.list(tht)) tht <- tht[[1]]',
            "dimnames(lam) <- list(rownames(lambda), colnames(lambda))",
            "dimnames(psi) <- list(colnames(lambda), colnames(lambda))",
            "dimnames(tht) <- list(rownames(lambda), rownames(lambda))",
            "semPlot::semPaths(semPlot::lisrelModel(LY = lam, PS = psi, TE = tht),",
            '  what = "std", whatLabels = "est", layout = "tree2", pastel = TRUE,',
            "  nCharNodes = 0, nCharEdges = 0)  # full labels, no abbreviation",
            sep = "\n")
        )
        fn(); return(invisible(NULL))
      }

      # --- Composite: factor structure + latent / RI network ----------------
      if (pt %in% c("factor_latent", "factor_ri")) {
        shiny::req(rv$model)
        spec <- build_factor_spec(
          rv$model, rv$lambda_used,
          type = if (pt == "factor_ri") "ri" else "latent",
          ri_norm = isTRUE(input$factor_ri_norm))
        shiny::validate(shiny::need(!is.null(spec),
          "Factor-structure plots need a fitted CFA/LNM/RNM/LRNM model (lambda + sigma_zeta)."))

        # Latent node size: scale by |mean factor score| (the latent severity
        # level), NOT observed column means (request #3).
        n_obs <- length(spec$labels) - spec$n_lat
        lat_vsize <- rep(s$vsize * 1.6, spec$n_lat)
        if (isTRUE(s$scale_nodes)) {
          fs <- mean_fscores_psynet(rv$model, rv$lambda_used, rv$dat)
          if (!is.null(fs) && length(fs) == spec$n_lat)
            lat_vsize <- scale_vsize_by_mean(abs(fs), s$vsize_min * 1.3,
                                             s$vsize_max * 1.3)
        }
        # Predictability rings on LATENTS only (analytic R^2 from sigma_zeta);
        # observed squares get no ring (NULL slots).
        pie_arg <- if (isTRUE(s$show_pred)) {
          sz <- tryCatch(collapse_mg(psychonetrics::getmatrix(rv$model, "sigma_zeta"), 1L),
                         error = function(e) NULL)
          r2l <- if (!is.null(sz)) latent_predictability_r2(sz) else NULL
          if (!is.null(r2l) && length(r2l) == spec$n_lat)
            c(as.list(r2l), vector("list", n_obs)) else NULL
        } else NULL

        fn <- function() do.call(qgraph::qgraph, c(
          list(spec$W,
          directed = spec$dir_mat, layout = spec$layout,
          labels = spec$labels, label.scale = FALSE,       # full names, equal size
          label.cex = s$label_cex,
          label.font = if (isTRUE(s$label_bold)) 2 else 1,
          shape = spec$shapes,
          vsize = c(lat_vsize, rep(s$vsize, n_obs)),
          color = spec$node_cols, border.color = s$node_border,
          posCol = s$pos_edge, negCol = s$neg_edge,
          esize = s$esize, asize = s$asize, minimum = s$min_edge,
          edge.labels = if (isTRUE(s$show_edge_labels)) spec$edge_lab else FALSE,
          edge.label.cex = s$edge_label_cex),
          house_pie_args(pie_arg, s)))
        rv$plot_fn <- fn
        rec_upsert(
          rec, "psynet_plot", "plot",
          description = sprintf(
            "[psychonetrics] Plotted the %s: latents as circles on an inner ring (one colour per factor), indicators as squares coloured by their dominant factor, loading arrows latent -> indicator.",
            spec$what),
          code = paste(c(
            'lambda  <- psychonetrics::getmatrix(mod, "lambda")',
            "if (is.list(lambda)) lambda <- lambda[[1]]",
            if (pt == "factor_ri") c(
              'sigma_z <- psychonetrics::getmatrix(mod, "sigma_zeta")',
              "if (is.list(sigma_z)) sigma_z <- sigma_z[[1]]",
              JOHNSON_RW_FRAGMENT,
              "lat_block <- johnson_rw_from_cor(cov2cor(as.matrix(sigma_z)))",
              if (isTRUE(input$factor_ri_norm))
                "cs <- colSums(lat_block); cs[cs < 1e-10] <- 1; lat_block <- sweep(lat_block, 2, cs, \"/\")")
            else c(
              'lat_block <- psychonetrics::getmatrix(mod, "omega_zeta")',
              "if (is.list(lat_block)) lat_block <- lat_block[[1]]"),
            "n_lat <- ncol(lambda); n_obs <- nrow(lambda); n <- n_lat + n_obs",
            "W <- matrix(0, n, n)",
            "W[1:n_lat, 1:n_lat] <- as.matrix(lat_block)",
            "W[1:n_lat, (n_lat + 1):n] <- t(as.matrix(lambda))",
            "dir_mat <- matrix(FALSE, n, n)",
            "dir_mat[1:n_lat, (n_lat + 1):n] <- TRUE",
            if (pt == "factor_ri") "dir_mat[1:n_lat, 1:n_lat] <- TRUE",
            "# latents = circles (inner ring), indicators = squares near their",
            "# dominant factor; labels are the full variable/factor names.",
            'qgraph::qgraph(W, directed = dir_mat, layout = "spring",',
            "  labels = c(colnames(lambda), rownames(lambda)),",
            sprintf("  label.scale = FALSE, label.cex = %s,", s$label_cex),
            sprintf('  shape = c(rep("circle", n_lat), rep("square", n_obs)),'),
            sprintf("  esize = %s, minimum = %s)", s$esize, s$min_edge)),
            collapse = "\n")
        )
        fn(); return(invisible(NULL))
      }

      # --- Simple single-matrix plots ----------------------------------------
      pm <- plot_matrix()
      shiny::validate(shiny::need(!is.null(pm$W), paste(
        "This matrix is not available for the fitted model —",
        "latent/RI plots need an LNM/RNM/LRNM/CFA fit;",
        "omega needs a GGM.")))

      is_latent_net <- pt %in% c("latent", "ri_raw", "ri_norm")
      W_list <- pm$W_list %||% list(pm$W)
      multi  <- length(W_list) > 1
      lay_choice <- input$psynet_layout %||% "spring"

      # Node size: latent networks scale by |mean factor score| (request #3);
      # observed networks scale by column means. (Single-group only — group-
      # specific scores/means would differ per panel.)
      vsize_arg <- NULL
      if (isTRUE(s$scale_nodes) && !multi) {
        vsize_arg <- if (is_latent_net) {
          fs <- mean_fscores_psynet(rv$model, rv$lambda_used, rv$dat)
          if (!is.null(fs) && length(fs) == ncol(pm$W))
            scale_vsize_by_mean(abs(fs), s$vsize_min, s$vsize_max) else NULL
        } else if (!is.null(rv$dat)) {
          means <- vapply(rv$dat[, intersect(colnames(pm$W), names(rv$dat)),
                                 drop = FALSE],
                          function(x) mean(x, na.rm = TRUE), numeric(1))
          if (length(means) == ncol(pm$W))
            scale_vsize_by_mean(means, s$vsize_min, s$vsize_max) else NULL
        }
      }

      # Predictability rings (single-group only, same reason as above).
      r2 <- if (isTRUE(s$show_pred) && !multi) {
        if (is_latent_net) {
          sz <- tryCatch(collapse_mg(psychonetrics::getmatrix(rv$model, "sigma_zeta"), 1L),
                         error = function(e) NULL)
          if (!is.null(sz)) latent_predictability_r2(sz) else NULL
        } else if (!is.null(rv$dat)) {
          keep <- intersect(colnames(pm$W), names(rv$dat))
          if (length(keep) == ncol(pm$W))
            node_predictability_r2(rv$dat[, keep, drop = FALSE]) else NULL
        }
      } else NULL

      # means under labels (single-group): factor scores for latent nets,
      # column means for observed nets.
      nvals <- if (!multi) {
        if (is_latent_net) {
          fs <- mean_fscores_psynet(rv$model, rv$lambda_used, rv$dat)
          if (!is.null(fs) && length(fs) == ncol(pm$W))
            stats::setNames(fs, colnames(pm$W)) else NULL
        } else if (!is.null(rv$dat)) {
          keep <- intersect(colnames(pm$W), names(rv$dat))
          if (length(keep) == ncol(pm$W))
            vapply(rv$dat[, keep, drop = FALSE], mean, numeric(1), na.rm = TRUE)
          else NULL
        }
      } else NULL

      args <- house_qgraph_args(pm$W, s, directed = pm$directed,
                                vsize = vsize_arg, pie = r2,
                                node_values = nvals)
      grp_labels <- {
        gv <- sel$group_var()
        lv <- if (!is.null(gv)) sort(unique(stats::na.omit(
                data_bus$wide()[[gv]]))) else NULL
        if (!is.null(lv) && length(lv) == length(W_list))
          as.character(lv) else paste("Group", seq_along(W_list))
      }
      fn <- if (multi) {
        # Side-by-side group networks on ONE shared layout so node positions
        # match across panels (request #9); groups named in panel titles.
        L <- if (lay_choice == "circle") "circle"
             else do.call(qgraph::averageLayout, W_list)
        function() {
          op <- graphics::par(mfrow = c(1, length(W_list)))
          on.exit(graphics::par(op), add = TRUE)
          for (k in seq_along(W_list))
            do.call(qgraph::qgraph,
                    c(list(W_list[[k]], layout = L, title = grp_labels[k]),
                      args))
        }
      } else {
        function() do.call(qgraph::qgraph,
                           c(list(pm$W, layout = lay_choice), args))
      }
      rv$plot_fn <- fn

      extract_code <- switch(pt,
        omega    = 'W <- psychonetrics::getmatrix(mod, "omega")',
        latent   = 'W <- psychonetrics::getmatrix(mod, "omega_zeta")',
        residual = 'W <- psychonetrics::getmatrix(mod, "omega_epsilon")',
        paste(c(
          'sigma_z <- psychonetrics::getmatrix(mod, "sigma_zeta")',
          "if (is.list(sigma_z) && !is.data.frame(sigma_z)) sigma_z <- sigma_z[[1]]",
          "R_lat <- cov2cor(as.matrix(sigma_z))",
          JOHNSON_RW_FRAGMENT,
          "W <- johnson_rw_from_cor(R_lat)",
          if (pt == "ri_norm") paste(
            "cs <- colSums(W); cs[cs < 1e-10] <- 1",
            'W <- sweep(W, 2, cs, "/")  # each column sums to 1', sep = "\n")),
          collapse = "\n"))

      rec_upsert(
        rec, "psynet_plot", "plot",
        description = sprintf(
          "[psychonetrics] Plotted the %s (qgraph, %s layout%s%s).", pm$what,
          lay_choice,
          if (pm$directed) "; directed: arrow j -> i = j's share of i's R-squared"
          else "",
          if (multi) sprintf("; %d groups side by side on one shared layout (%s)",
                             length(W_list), paste(grp_labels, collapse = ", "))
          else ""),
        code = paste(c(
          extract_code,
          if (multi) c(
            "W_list <- if (is.list(W) && !is.data.frame(W)) lapply(W, as.matrix) else list(as.matrix(W))",
            "L <- do.call(qgraph::averageLayout, W_list)  # shared layout across groups",
            "op <- par(mfrow = c(1, length(W_list))); on.exit(par(op), add = TRUE)",
            "for (k in seq_along(W_list)) qgraph::qgraph(W_list[[k]], layout = L,",
            sprintf('  title = %s[k],', vars_literal(grp_labels)),
            house_qgraph_args_code(s, pm$directed, wobj = "W_list[[k]]"),
            ")")
          else c(
            "if (is.list(W) && !is.data.frame(W)) W <- W[[1]]",
            "W <- as.matrix(W)",
            sprintf('qgraph::qgraph(W, layout = "%s",', lay_choice),
            house_qgraph_args_code(s, pm$directed, wobj = "W"),
            ")")), collapse = "\n")
      )
      fn()
    })

    list(model = shiny::reactive(rv$model))
  })
}

# ============================================================================
# NCT sub-module — paired handling + subject-ID matching + p correction
# ============================================================================

nctUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(4,
      # show_group = FALSE: NCT has its own required grouping selector below,
      # so the module's optional one would be redundant.
      varselectUI(ns("vars"), "Nodes to compare", show_group = FALSE),
      shiny::selectInput(ns("group_col"),
                         "Grouping variable (defines the two groups)",
                         choices = NULL),
      shiny::radioButtons(ns("design"), "Design (required — determines paired)",
                          c("Independent groups (paired = FALSE)" = "independent",
                            "Pre-post / repeated measures (paired = TRUE)" = "paired")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'paired'", ns("design")),
        shiny::selectInput(ns("id_col"), "Subject-ID column (aligns the pairs)",
                           choices = c("(none — assume matching row order)" = "")),
        shiny::helpText("With an ID column, rows are matched subject-by-subject",
                        "and unmatched subjects dropped (reported). Without one,",
                        "equal group sizes and matching row order are ASSUMED.")
      ),
      shiny::selectInput(ns("nct_est"), "Network estimator",
                         c("Auto (IsingFit if binary, else EBICglasso)" = "auto",
                           "EBICglasso" = "EBICglasso",
                           "ggmModSelect" = "ggmModSelect",
                           "pcor (partial correlations)" = "pcor",
                           "cor (marginal correlations)" = "cor",
                           "IsingFit (binary)" = "IsingFit"),
                         selected = "auto"),
      shiny::selectInput(ns("adjust"), "Edge-test p-value correction",
                         c("Holm (recommended)" = "holm", "BH / FDR" = "BH",
                           "Bonferroni" = "bonferroni",
                           "None (NOT recommended)" = "none"),
                         selected = "holm"),
      shiny::numericInput(ns("iterations"), "Permutation iterations",
                          value = 1000, min = 100, step = 100),
      shiny::actionButton(ns("run"), "Run NCT", class = "btn-primary"),
      shiny::helpText("Both group networks use the SAME chosen estimator.",
                      "A null NCT is not proof of equality.")
    ),
    shiny::column(8,
      shiny::h4("Group networks (same layout, so they are directly comparable)"),
      shiny::plotOutput(ns("group_plots"), height = "360px"),
      shiny::hr(),
      shiny::verbatimTextOutput(ns("nct_summary")),
      DT::dataTableOutput(ns("edge_table"))
    )
  )
}

# Estimate one network from a data matrix with any of the supported
# estimators (request #4). `est` in {auto, EBICglasso, ggmModSelect, pcor,
# cor, IsingFit}; "auto" -> IsingFit if binary else EBICglasso.
nct_estimate_one <- function(m, est, gamma, binary = FALSE) {
  if (identical(est, "auto")) est <- if (binary) "IsingFit" else "EBICglasso"
  net <- switch(est,
    IsingFit     = bootnet::estimateNetwork(m, default = "IsingFit", tuning = gamma),
    EBICglasso   = bootnet::estimateNetwork(m, default = "EBICglasso",
                                            corMethod = "cor_auto", tuning = gamma),
    ggmModSelect = bootnet::estimateNetwork(m, default = "ggmModSelect",
                                            corMethod = "cor_auto", tuning = gamma),
    pcor         = bootnet::estimateNetwork(m, default = "pcor",
                                            corMethod = "cor_auto"),
    cor          = bootnet::estimateNetwork(m, default = "cor",
                                            corMethod = "cor_auto"),
    bootnet::estimateNetwork(m, default = "EBICglasso",
                             corMethod = "cor_auto", tuning = gamma))
  qgraph::getWmat(net)
}

# Self-contained permutation Network Comparison Test — the fallback used when
# the installed NetworkComparisonTest package errors (e.g. the version-
# specific "comparison (==)..." bug). Tests, by permuting group labels and
# re-estimating with the chosen estimator:
#   * global strength invariance  (sum of |edge weights|)
#   * network structure invariance (max |edge difference|)
#   * EVERY individual edge (|difference| per edge, request #5) — with the
#     requested multiple-comparison correction.
manual_nct <- function(m1, m2, it, est, gamma, binary, paired, seed,
                       adjust = "holm") {
  set.seed(seed)
  W1 <- nct_estimate_one(m1, est, gamma, binary)
  W2 <- nct_estimate_one(m2, est, gamma, binary)
  ut    <- upper.tri(W1)
  glstr <- function(W) sum(abs(W[ut]))
  s_obs <- abs(glstr(W1) - glstr(W2))
  d_obs <- max(abs(W1 - W2))
  e_obs <- abs(W1[ut] - W2[ut])                 # observed |diff| per edge
  nm    <- colnames(W1)
  pairs <- which(ut, arr.ind = TRUE)
  pooled <- rbind(as.matrix(m1), as.matrix(m2))
  n1 <- nrow(m1); N <- nrow(pooled)
  s_perm <- numeric(it); d_perm <- numeric(it)
  e_ge   <- numeric(length(e_obs)); e_n <- 0L    # count perms with |diff|>=obs
  for (b in seq_len(it)) {
    if (paired && nrow(m1) == nrow(m2)) {
      swap <- stats::rbinom(n1, 1, 0.5) == 1
      pm1 <- rbind(as.matrix(m1)[!swap, , drop = FALSE],
                   as.matrix(m2)[swap, , drop = FALSE])
      pm2 <- rbind(as.matrix(m2)[!swap, , drop = FALSE],
                   as.matrix(m1)[swap, , drop = FALSE])
    } else {
      idx <- sample.int(N)
      pm1 <- pooled[idx[seq_len(n1)], , drop = FALSE]
      pm2 <- pooled[idx[(n1 + 1):N], , drop = FALSE]
    }
    PW1 <- tryCatch(nct_estimate_one(pm1, est, gamma, binary), error = function(e) NULL)
    PW2 <- tryCatch(nct_estimate_one(pm2, est, gamma, binary), error = function(e) NULL)
    if (is.null(PW1) || is.null(PW2)) { s_perm[b] <- NA; d_perm[b] <- NA; next }
    s_perm[b] <- abs(glstr(PW1) - glstr(PW2))
    d_perm[b] <- max(abs(PW1 - PW2))
    e_ge <- e_ge + (abs(PW1[ut] - PW2[ut]) >= e_obs)
    e_n  <- e_n + 1L
  }
  edge_p <- if (e_n > 0) (e_ge + 1) / (e_n + 1) else rep(NA_real_, length(e_obs))
  edge_p <- stats::p.adjust(edge_p, method = adjust)
  einv <- data.frame(
    Var1 = nm[pairs[, 1]], Var2 = nm[pairs[, 2]],
    diff = e_obs, p = edge_p, stringsAsFactors = FALSE)
  einv <- einv[order(einv$p, -einv$diff), ]
  structure(list(
    glstrinv.real = s_obs,
    glstrinv.pval = mean(c(s_perm, s_obs) >= s_obs, na.rm = TRUE),
    nwinv.real    = d_obs,
    nwinv.pval    = mean(c(d_perm, d_obs) >= d_obs, na.rm = TRUE),
    einv.pvals    = einv,
    nperm = sum(!is.na(s_perm))),
    class = "manual_nct")
}

print.manual_nct <- function(x, ...) {
  cat("Fallback permutation NCT\n")
  cat(sprintf("  Global strength invariance:   diff = %.4f, p = %.4f\n",
              x$glstrinv.real, x$glstrinv.pval))
  cat(sprintf("  Network structure invariance: max|diff| = %.4f, p = %.4f\n",
              x$nwinv.real, x$nwinv.pval))
  if (!is.null(x$einv.pvals)) {
    sig <- sum(x$einv.pvals$p < 0.05, na.rm = TRUE)
    cat(sprintf("  Edge-level tests: %d of %d edges differ (corrected p < .05); see the table.\n",
                sig, nrow(x$einv.pvals)))
  }
  cat(sprintf("  (%d valid permutations)\n", x$nperm))
  invisible(x)
}
# Return the object (no printing) so an outer print() renders it exactly once
# — the previous version printed inside summary AND via the outer print,
# which duplicated the output.
summary.manual_nct <- function(object, ...) invisible(object)

# Core NCT call. Prefers the installed package; if it errors (version-specific
# bugs such as the "comparison (==) is possible only for atomic and list
# types" crash), falls back to the self-contained permutation test above.
# Data are coerced to a plain numeric matrix first (data frames with any odd
# column class are a common trigger). No estimator FUNCTION is ever passed.
run_nct_corrected <- function(dat1, dat2, paired, it, adjust, gg, seed,
                              est = "auto") {
  m1 <- data.matrix(dat1); m2 <- data.matrix(dat2)
  all_binary <- all(apply(rbind(m1, m2), 2,
                          function(x) length(unique(stats::na.omit(x))) <= 2))
  use_binary <- all_binary || identical(est, "IsingFit")
  res <- tryCatch({
    set.seed(seed)
    nct_formals <- names(formals(NetworkComparisonTest::NCT))
    args <- list(data1 = as.data.frame(m1), data2 = as.data.frame(m2),
                 it = it, paired = paired, test.edges = TRUE, edges = "all")
    if ("gamma" %in% nct_formals) args$gamma <- gg$tuning
    if (use_binary && "binary.data" %in% nct_formals) args$binary.data <- TRUE
    if ("p.adjust.methods" %in% nct_formals) args$p.adjust.methods <- adjust
    if ("progressbar" %in% nct_formals) args$progressbar <- FALSE  # no console bar
    # swallow any residual console progress printing from NCT/bootnet
    r <- NULL
    invisible(utils::capture.output(suppressMessages(
      r <- do.call(NetworkComparisonTest::NCT, args))))
    if (!("p.adjust.methods" %in% nct_formals) &&
        !identical(adjust, "none") && !is.null(r$einv.pvals)) {
      pcol <- grep("p", names(r$einv.pvals), ignore.case = TRUE, value = TRUE)[1]
      r$einv.pvals[[pcol]] <- stats::p.adjust(r$einv.pvals[[pcol]], method = adjust)
    }
    attr(r, "engine") <- "NetworkComparisonTest"
    r
  }, error = function(e) {
    r <- NULL
    invisible(utils::capture.output(suppressMessages(
      r <- manual_nct(m1, m2, min(it, 500L), est, gg$tuning, use_binary, paired,
                      seed, adjust = adjust))))
    attr(r, "engine") <- "fallback"
    attr(r, "pkg_error") <- conditionMessage(e)
    r
  })
  attr(res, "edge_p_adjust") <- adjust
  attr(res, "binary_data")   <- all_binary
  res
}

# Align two paired data frames by a subject-ID vector. Keeps subjects present
# exactly once in EACH group, in matching order; reports what was dropped.
align_paired_by_id <- function(d1, d2, ids1, ids2) {
  keep1 <- ids1[!duplicated(ids1)]
  keep2 <- ids2[!duplicated(ids2)]
  common <- intersect(keep1, keep2)
  list(
    d1 = d1[match(common, ids1), , drop = FALSE],
    d2 = d2[match(common, ids2), , drop = FALSE],
    n_matched = length(common),
    n_dropped = (nrow(d1) - length(common)) + (nrow(d2) - length(common))
  )
}

nctServer <- function(id, data_bus, rec, gg_settings) {
  shiny::moduleServer(id, function(input, output, session) {

    sel <- varselectServer("vars", data_bus, rec, rec_prefix = "nct")

    shiny::observeEvent(data_bus$wide(), {
      shiny::updateSelectInput(session, "group_col",
                               choices = names(data_bus$wide()))
      shiny::updateSelectInput(session, "id_col",
        choices = c("(none — assume matching row order)" = "",
                    names(data_bus$wide())))
    })

    # prep_r: build the two group data frames — NO NCT here, so the
    # side-by-side group networks always render even if the NCT test errors.
    prep_r <- shiny::eventReactive(input$run, {
      dat <- data_bus$wide()
      shiny::req(input$group_col %in% names(dat))
      g   <- dat[[input$group_col]]
      lv  <- names(sort(table(g), decreasing = TRUE))[1:2]
      shiny::validate(shiny::need(length(lv) == 2 && !anyNA(lv),
        "The grouping variable needs at least two non-missing groups."))
      num <- setdiff(sel$vars(), c(input$group_col, input$id_col))
      num <- intersect(num, names(dat)[vapply(dat, is.numeric, TRUE)])
      shiny::validate(shiny::need(length(num) >= 2,
        "Select at least 2 numeric nodes to compare."))
      trans <- sel$transform()
      dat[, num] <- apply_house_transform(dat[, num, drop = FALSE], trans)
      d1  <- dat[g == lv[1], num, drop = FALSE]
      d2  <- dat[g == lv[2], num, drop = FALSE]
      paired  <- identical(input$design, "paired")
      id_note <- ""
      if (paired) {
        if (!is.null(input$id_col) && nzchar(input$id_col) &&
            input$id_col %in% names(dat)) {
          al <- align_paired_by_id(d1, d2,
                                   dat[[input$id_col]][g == lv[1]],
                                   dat[[input$id_col]][g == lv[2]])
          shiny::validate(shiny::need(al$n_matched >= 10,
            "Fewer than 10 matched subjects across the two conditions."))
          if (al$n_dropped > 0) shiny::showNotification(
            sprintf("Paired matching by '%s': %d matched pairs; %d unmatched rows dropped.",
                    input$id_col, al$n_matched, al$n_dropped),
            type = "warning", duration = 8)
          d1 <- al$d1; d2 <- al$d2
          id_note <- sprintf(" rows matched by ID column '%s' (%d pairs);",
                             input$id_col, al$n_matched)
        } else {
          shiny::validate(shiny::need(nrow(d1) == nrow(d2),
            "Paired design without an ID column needs equal group sizes (row order is assumed aligned) — or pick a subject-ID column."))
          id_note <- " NO ID column: matching row order ASSUMED;"
        }
      }
      binary <- all(apply(rbind(data.matrix(d1), data.matrix(d2)), 2,
                          function(x) length(unique(stats::na.omit(x))) <= 2))
      list(d1 = d1, d2 = d2, groups = lv, num = num, paired = paired,
           id_note = id_note, trans = trans, binary = binary,
           est = input$nct_est %||% "auto",
           gg = gg_settings(), seed = rec_seed(rec))
    })

    res_r <- shiny::eventReactive(input$run, {
      p <- prep_r()
      lv <- p$groups; num <- p$num; trans <- p$trans
      paired <- p$paired; id_note <- p$id_note; gg <- p$gg; seed <- p$seed
      res <- shiny::withProgress(
        message = "Running network comparison",
        detail = "estimating + permuting (this can take a minute)...",
        value = 0.5,
        run_nct_corrected(p$d1, p$d2, paired, input$iterations,
                          input$adjust, gg, seed, est = p$est))
      if (identical(attr(res, "engine"), "fallback"))
        shiny::showNotification(paste0(
          "Your installed NetworkComparisonTest package errored (see the ",
          "results box for the reason), so XS4ALL ran its BUILT-IN permutation ",
          "NCT instead — same method (van Borkulo et al.): global strength, ",
          "global structure, and per-edge tests. Centrality-difference tests ",
          "still need the package. To use the package itself, update it in R: ",
          "remotes::install_github('cvborkulo/NetworkComparisonTest')."),
          type = "warning", duration = 15)

      engine_note <- if (identical(attr(res, "engine"), "fallback"))
        "built-in permutation fallback (package errored)" else
        if (isTRUE(p$binary)) "IsingFit, binary.data = TRUE"
        else "EBICglasso on cor_auto"
      rec_upsert(
        rec, "nct_comparison", "comparison",
        description = sprintf(
          "[NCT] Compared '%s' vs '%s' (%s design, paired = %s):%s %d permutations, seed %d; transform: %s; edge-level p-values corrected with '%s'; estimation = %s, gamma = %s. A null NCT is not proof of equality.",
          lv[1], lv[2], input$design, paired, id_note,
          input$iterations, seed,
          names(TRANSFORM_LABELS)[TRANSFORM_LABELS == trans],
          input$adjust, engine_note, gg$tuning),
        code = paste(c(
          sprintf("nct_vars <- %s", vars_literal(num)),
          "dat_nct <- dat_wide[, nct_vars]",
          transform_code_fragment(trans, "dat_nct"),
          sprintf('g <- dat_wide[["%s"]]', input$group_col),
          sprintf('d1 <- dat_nct[g == "%s", ]', lv[1]),
          sprintf('d2 <- dat_nct[g == "%s", ]', lv[2]),
          if (paired && nzchar(id_note) && grepl("matched by ID", id_note)) c(
            sprintf('ids1 <- dat_wide[["%s"]][g == "%s"]', input$id_col, lv[1]),
            sprintf('ids2 <- dat_wide[["%s"]][g == "%s"]', input$id_col, lv[2]),
            "common <- intersect(ids1[!duplicated(ids1)], ids2[!duplicated(ids2)])",
            "d1 <- d1[match(common, ids1), ]; d2 <- d2[match(common, ids2), ]"),
          sprintf("set.seed(%d)", seed),
          "nct_res <- NetworkComparisonTest::NCT(d1, d2,",
          sprintf("  it = %d, paired = %s, gamma = %s,",
                  input$iterations, paired, gg$tuning),
          if (isTRUE(p$binary))
            "  binary.data = TRUE,  # all nodes binary -> IsingFit path",
          '  test.edges = TRUE, edges = "all",',
          sprintf('  p.adjust.methods = "%s")  # edge-level correction',
                  input$adjust)), collapse = "\n")
      )
      list(res = res, groups = lv)
    })

    # --- Side-by-side group networks, SHARED layout (independent of NCT) ----
    output$group_plots <- shiny::renderPlot({
      p <- prep_r()
      est <- function(d) tryCatch({
        W <- NULL
        invisible(utils::capture.output(suppressMessages(
          W <- nct_estimate_one(data.matrix(d), p$est, p$gg$tuning, p$binary))))
        W
      }, error = function(e) NULL)
      W1 <- est(p$d1); W2 <- est(p$d2)
      shiny::validate(shiny::need(!is.null(W1) && !is.null(W2),
        "Could not estimate one of the group networks (a node may be constant within a group)."))
      pal <- house_pastel()
      L <- qgraph::averageLayout(W1, W2)   # shared layout across panels
      op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
      qgraph::qgraph(W1, layout = L, labels = colnames(W1), label.scale = FALSE,
                     posCol = pal$pos_edge, negCol = pal$neg_edge,
                     color = pal$node_fill, border.color = pal$node_border,
                     title = paste("Group:", p$groups[1]))
      qgraph::qgraph(W2, layout = L, labels = colnames(W2), label.scale = FALSE,
                     posCol = pal$pos_edge, negCol = pal$neg_edge,
                     color = pal$node_fill, border.color = pal$node_border,
                     title = paste("Group:", p$groups[2]))
    })

    output$nct_summary <- shiny::renderPrint({
      r <- res_r()
      cat(sprintf("Comparison: %s vs %s | paired = %s | engine = %s\n",
                  r$groups[1], r$groups[2],
                  identical(input$design, "paired"),
                  attr(r$res, "engine")))
      if (!is.null(attr(r$res, "pkg_error")))
        cat(sprintf("Package error (fell back): %s\n",
                    attr(r$res, "pkg_error")))
      cat("\n")
      print(summary(r$res))
    })

    output$edge_table <- DT::renderDataTable({
      r <- res_r()
      shiny::validate(shiny::need(!is.null(r$res$einv.pvals),
        "Edge-level tests are only available from the NetworkComparisonTest package (not the permutation fallback)."))
      tb <- as.data.frame(r$res$einv.pvals)
      pcol <- grep("p", names(tb), ignore.case = TRUE, value = TRUE)[1]
      tb$significant <- tb[[pcol]] < 0.05   # on CORRECTED p-values
      DT::datatable(tb, caption = sprintf(
        "Edge-level differences; p-values %s-corrected.",
        attr(r$res, "edge_p_adjust")))
    })

    res_r
  })
}
