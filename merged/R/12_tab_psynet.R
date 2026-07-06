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
      shiny::column(8, shiny::plotOutput(ns("psynet_plot"), height = "480px")),
      shiny::column(4,
        shiny::selectInput(ns("plot_type"), "Network to plot",
          c("Observed GGM (omega)"                 = "omega",
            "Latent network (omega_zeta)"          = "latent",
            "RI latent network (raw)"              = "ri_raw",
            "RI latent network (normalized)"       = "ri_norm",
            "Residual network (omega_epsilon)"     = "residual")),
        shiny::helpText("RI = Relative-Importance latent network",
                        "(Johnson/LMG): a DIRECTED network where the arrow",
                        "j -> i shows latent j's share of the R-squared of",
                        "latent i. Multi-group models show group 1."),
        # edge labels ON by default so the RI network shows its weights
        # as numbers on the arrows rather than hiding them (request #4).
        appearanceUI(ns("look"), signed = TRUE, default_edge_labels = TRUE)
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
          shiny::tags$th(l, style = "text-align:center;")))
      rows <- lapply(seq_along(vars), function(i) shiny::tags$tr(
        shiny::tags$td(shiny::tags$b(vars[i])),
        lapply(seq_along(latents), function(j) shiny::tags$td(
          shiny::checkboxInput(ns(sprintf("lam_%d_%d", i, j)), NULL,
                               value = as.logical(lam[i, j])),
          style = "text-align:center; padding:0 6px;"))))
      shiny::tags$div(style = "overflow-x:auto;",
        shiny::tags$table(class = "table table-condensed table-bordered",
                          shiny::tags$thead(header),
                          shiny::tags$tbody(rows)))
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
      rv$model <- mod; rv$family <- fam
      rv$estimator_used <- est; rv$lambda_used <- lambda

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
          "[psychonetrics] Fitted %s on %d variables (n = %d): estimator = %s (resolved explicitly%s), optimizer = nlminb, data type = %s%s%s%s.",
          toupper(fam), length(vars), nrow(dat), est,
          if (input$estimator == "auto") " from 'auto'" else "",
          input$data_type,
          if (fam != "ggm") sprintf(", identification = %s, %d latents",
                                    input$identification, ncol(lambda)) else "",
          if (isTRUE(input$do_prune))
            sprintf("; pruned at alpha = %s (adjust = %s)",
                    input$prune_alpha, input$prune_adjust) else "",
          if (isTRUE(input$do_stepup)) "; step-up search" else ""),
        code = paste(c(
          sprintf("psynet_vars <- %s", vars_literal(vars)),
          "dat_psynet  <- dat_wide[, psynet_vars]",
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

    # Extract + prepare the matrix for the chosen plot type; NULL if absent.
    plot_matrix <- shiny::reactive({
      shiny::req(rv$model)
      pt <- input$plot_type
      get_mat <- function(name) tryCatch(
        collapse_mg(psychonetrics::getmatrix(rv$model, name), 1L),
        error = function(e) NULL)
      if (pt == "omega")    return(list(W = get_mat("omega"), directed = FALSE,
                                        what = "observed GGM (omega)"))
      if (pt == "latent")   return(list(W = get_mat("omega_zeta"), directed = FALSE,
                                        what = "latent network (omega_zeta)"))
      if (pt == "residual") return(list(W = get_mat("omega_epsilon"), directed = FALSE,
                                        what = "residual network (omega_epsilon)"))
      # RI networks from sigma_zeta -> latent correlations -> Johnson/LMG
      sz <- get_mat("sigma_zeta")
      if (is.null(sz)) return(list(W = NULL))
      R_lat <- tryCatch(stats::cov2cor(sz), error = function(e) NULL)
      if (is.null(R_lat)) return(list(W = NULL))
      lat <- colnames(rv$lambda_used)
      if (!is.null(lat) && length(lat) == ncol(R_lat))
        dimnames(R_lat) <- list(lat, lat)
      RI <- johnson_rw_from_cor(R_lat)
      if (pt == "ri_norm") {
        cs <- colSums(RI); cs[cs < 1e-10] <- 1
        RI <- sweep(RI, 2, cs, "/")
      }
      list(W = RI, directed = TRUE,
           what = sprintf("RI latent network (%s)",
                          if (pt == "ri_norm") "normalized" else "raw"))
    })

    output$psynet_plot <- shiny::renderPlot({
      pm <- plot_matrix()
      shiny::validate(shiny::need(!is.null(pm$W), paste(
        "This matrix is not available for the fitted model —",
        "latent/RI plots need an LNM/RNM/LRNM/CFA fit;",
        "omega needs a GGM.")))
      s  <- look()
      pt <- input$plot_type

      args <- house_qgraph_args(pm$W, s, directed = pm$directed)
      fn <- function() do.call(qgraph::qgraph,
                               c(list(pm$W, layout = "spring"), args))
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
          "[psychonetrics] Plotted the %s (qgraph%s).", pm$what,
          if (pm$directed) ", directed: arrow j -> i = j's share of i's R-squared"
          else ""),
        code = paste(
          extract_code,
          "if (is.list(W) && !is.data.frame(W)) W <- W[[1]]  # multi-group: group 1",
          "W <- as.matrix(W)",
          'qgraph::qgraph(W, layout = "spring",',
          house_qgraph_args_code(s, pm$directed, wobj = "W"),
          ")", sep = "\n")
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
      varselectUI(ns("vars"), "Nodes to compare"),
      shiny::selectInput(ns("group_col"), "Grouping column", choices = NULL),
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
      shiny::selectInput(ns("adjust"), "Edge-test p-value correction",
                         c("Holm (recommended)" = "holm", "BH / FDR" = "BH",
                           "Bonferroni" = "bonferroni",
                           "None (NOT recommended)" = "none"),
                         selected = "holm"),
      shiny::numericInput(ns("iterations"), "Permutation iterations",
                          value = 1000, min = 100, step = 100),
      shiny::actionButton(ns("run"), "Run NCT", class = "btn-primary"),
      shiny::helpText("Networks are estimated with the GGM tab's settings",
                      "(default = EBICglasso). A null NCT is not proof of equality.")
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

# Core NCT call. Pins estimation to the shared gg settings and applies the
# chosen edge-level correction natively when NCT supports it, manually if not.
run_nct_corrected <- function(dat1, dat2, paired, it, adjust, gg, seed) {
  set.seed(seed)
  nct_formals <- names(formals(NetworkComparisonTest::NCT))
  args <- list(data1 = dat1, data2 = dat2, it = it, paired = paired,
               test.edges = TRUE, edges = "all")
  # Pin estimation so NCT uses the same settings as the bootnet tab:
  if (all(c("estimator", "estimatorArgs") %in% nct_formals)) {
    args$estimator <- bootnet::estimateNetwork
    args$estimatorArgs <- list(default = gg$default,
                               corMethod = gg$corMethod, tuning = gg$tuning)
  } else if ("gamma" %in% nct_formals) {
    args$gamma <- gg$tuning
  }
  native_adjust <- "p.adjust.methods" %in% nct_formals
  if (native_adjust) args$p.adjust.methods <- adjust
  res <- do.call(NetworkComparisonTest::NCT, args)
  # Manual fallback for NCT builds without p.adjust.methods:
  if (!native_adjust && !identical(adjust, "none") &&
      !is.null(res$einv.pvals)) {
    pcol <- grep("p", names(res$einv.pvals), ignore.case = TRUE, value = TRUE)[1]
    res$einv.pvals[[pcol]] <- stats::p.adjust(res$einv.pvals[[pcol]],
                                              method = adjust)
  }
  attr(res, "edge_p_adjust") <- adjust
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

    res_r <- shiny::eventReactive(input$run, {
      dat <- data_bus$wide()
      shiny::req(input$group_col %in% names(dat))
      g   <- dat[[input$group_col]]
      lv  <- names(sort(table(g), decreasing = TRUE))[1:2]
      # nodes: the user's selection, minus the group/id columns
      num <- setdiff(sel$vars(), c(input$group_col, input$id_col))
      num <- intersect(num, names(dat)[vapply(dat, is.numeric, TRUE)])
      shiny::validate(shiny::need(length(num) >= 2,
        "Select at least 2 numeric nodes to compare."))
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
      gg   <- gg_settings()
      seed <- rec_seed(rec)
      res  <- run_nct_corrected(d1, d2, paired, input$iterations,
                                input$adjust, gg, seed)

      rec_upsert(
        rec, "nct_comparison", "comparison",
        description = sprintf(
          "[NCT] Compared '%s' vs '%s' (%s design, paired = %s):%s %d permutations, seed %d; edge-level p-values corrected with '%s'; estimation pinned to %s / %s / gamma = %s. A null NCT is not proof of equality.",
          lv[1], lv[2], input$design, paired, id_note,
          input$iterations, seed, input$adjust,
          gg$default, gg$corMethod, gg$tuning),
        code = paste(c(
          sprintf("nct_vars <- %s", vars_literal(num)),
          sprintf('g <- dat_wide[["%s"]]', input$group_col),
          sprintf('d1 <- dat_wide[g == "%s", nct_vars]', lv[1]),
          sprintf('d2 <- dat_wide[g == "%s", nct_vars]', lv[2]),
          if (paired && nzchar(id_note) && grepl("matched by ID", id_note)) c(
            sprintf('ids1 <- dat_wide[["%s"]][g == "%s"]', input$id_col, lv[1]),
            sprintf('ids2 <- dat_wide[["%s"]][g == "%s"]', input$id_col, lv[2]),
            "common <- intersect(ids1[!duplicated(ids1)], ids2[!duplicated(ids2)])",
            "d1 <- d1[match(common, ids1), ]; d2 <- d2[match(common, ids2), ]"),
          sprintf("set.seed(%d)", seed),
          "nct_res <- NetworkComparisonTest::NCT(d1, d2,",
          sprintf("  it = %d, paired = %s,", input$iterations, paired),
          '  test.edges = TRUE, edges = "all",',
          sprintf('  p.adjust.methods = "%s")  # edge-level correction',
                  input$adjust)), collapse = "\n")
      )
      list(res = res, groups = lv, d1 = d1, d2 = d2, gg = gg, num = num)
    })

    # --- Side-by-side group networks, SHARED layout (request #10) -----------
    output$group_plots <- shiny::renderPlot({
      r  <- res_r()
      gg <- r$gg
      est <- function(d) tryCatch(qgraph::getWmat(
        bootnet::estimateNetwork(d, default = gg$default,
                                 corMethod = gg$corMethod, tuning = gg$tuning)),
        error = function(e) NULL)
      W1 <- est(r$d1); W2 <- est(r$d2)
      shiny::validate(shiny::need(!is.null(W1) && !is.null(W2),
        "Could not estimate one of the group networks."))
      pal <- house_pastel()
      # one averaged layout so node positions match across the two panels
      L <- qgraph::averageLayout(W1, W2)
      op <- graphics::par(mfrow = c(1, 2)); on.exit(graphics::par(op), add = TRUE)
      qgraph::qgraph(W1, layout = L, labels = colnames(W1), label.scale = FALSE,
                     posCol = pal$pos_edge, negCol = pal$neg_edge,
                     color = pal$node_fill, border.color = pal$node_border,
                     title = paste("Group:", r$groups[1]))
      qgraph::qgraph(W2, layout = L, labels = colnames(W2), label.scale = FALSE,
                     posCol = pal$pos_edge, negCol = pal$neg_edge,
                     color = pal$node_fill, border.color = pal$node_border,
                     title = paste("Group:", r$groups[2]))
    })

    output$nct_summary <- shiny::renderPrint({
      r <- res_r()
      cat(sprintf("NCT: %s vs %s | paired = %s | edge p-adjust = %s\n\n",
                  r$groups[1], r$groups[2],
                  identical(input$design, "paired"),
                  attr(r$res, "edge_p_adjust")))
      print(summary(r$res))
    })

    output$edge_table <- DT::renderDataTable({
      r <- res_r()
      shiny::req(!is.null(r$res$einv.pvals))
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
