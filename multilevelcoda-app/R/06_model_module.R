# ============================================================================
# 06_model_module.R — complr() + brmcoda() (Le et al., 2025)
#
# Composition as predictor: outcome ~ between-person ilrs (bz) + within-person
# ilrs (wz) + covariates + (1 | ID). With one participant (a single patient's
# diary) the model is single-level: outcome ~ z + covariates; nights are then
# treated as independent, which they are not (serial dependence), so results
# are descriptive for that person.
# Composition as outcome: mvbind(z1, z2, ...) ~ predictors + (1 | ID).
# ============================================================================

modelModuleUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      width = 4,
      shiny::radioButtons(ns("role"), "Role of the composition",
        c("Predictor of an outcome (enables substitution)" = "predictor",
          "Outcome (multivariate model of the ilrs)" = "outcome")),
      shiny::selectizeInput(ns("order"), "Part order \u2014 drag to reorder",
        choices = NULL, multiple = TRUE,
        options = list(plugins = list("drag_drop", "remove_button"))),
      shiny::helpText("Default sequential binary partition = pivot balances:",
        "z1 contrasts the FIRST part with all others, z2 the second with the",
        "remaining ones, and so on. Put the part you want to interpret first."),
      shiny::tableOutput(ns("sbp")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'predictor'", ns("role")),
        shiny::selectInput(ns("outcome"), "Outcome", choices = NULL)),
      shinyWidgets::pickerInput(ns("covars"), "Covariates / predictors",
        choices = NULL, multiple = TRUE,
        options = list(`live-search` = TRUE, `actions-box` = TRUE)),
      shiny::uiOutput(ns("level_ui")),
      shiny::selectInput(ns("family"), "Likelihood",
        c("Gaussian" = "gaussian", "Student-t (robust to outliers)" = "student")),
      shiny::h5("Sampler"),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("chains"), "Chains", 4, 1, 8)),
        shiny::column(6, shiny::numericInput(ns("cores"), "Cores",
                                             min(4, max(1, parallel::detectCores() - 1)), 1, 64))),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("iter"), "Iterations", 2000, 200, 50000, 500)),
        shiny::column(6, shiny::numericInput(ns("warmup"), "Warm-up", 1000, 100, 25000, 500))),
      shiny::fluidRow(
        shiny::column(6, shiny::numericInput(ns("adapt_delta"), "adapt_delta", 0.95, 0.8, 0.999, 0.01)),
        shiny::column(6, shiny::numericInput(ns("seed"), "Seed", 2025, 1))),
      shiny::selectInput(ns("backend"), "Stan backend",
        c("rstan" = "rstan", "cmdstanr" = "cmdstanr")),
      shiny::actionButton(ns("fit"), "Fit model", class = "btn-primary"),
      shiny::helpText("The model is compiled by Stan first (about a minute), then",
                      "sampled. The app is busy until it finishes.")
    ),
    shiny::mainPanel(
      width = 8,
      shiny::h4("Model"),
      shiny::verbatimTextOutput(ns("formula")),
      shiny::uiOutput(ns("design_note")),
      shiny::h4("Convergence"),
      shiny::tableOutput(ns("diag")),
      shiny::h4("Summary"),
      shiny::verbatimTextOutput(ns("summary")),
      shiny::h4("Posterior predictive check"),
      shiny::plotOutput(ns("ppc"), height = "320px")
    )
  )
}

modelModuleServer <- function(id, comp, rec) {
  shiny::moduleServer(id, function(input, output, session) {
    ana <- comp$ana

    multilevel <- shiny::reactive(length(unique(ana()$ID)) > 1)

    shiny::observe({
      p <- comp$parts()
      cur <- shiny::isolate(input$order)
      sel <- if (setequal(cur, p)) cur else p
      shiny::updateSelectizeInput(session, "order", choices = p, selected = sel)
    })
    parts_ordered <- shiny::reactive({
      p <- comp$parts()
      o <- input$order
      if (length(o) == length(p) && setequal(o, p)) o else p
    })
    output$sbp <- shiny::renderTable(
      as.data.frame(multilevelcoda::build.sbp(parts_ordered()),
                    row.names = paste0("z", seq_len(length(parts_ordered()) - 1))),
      rownames = TRUE, digits = 0)

    vars <- shiny::reactive({
      a <- ana()
      setdiff(names(a)[vapply(a, is.numeric, logical(1))],
              c("night", comp$parts()))
    })
    shiny::observe({
      v <- vars()
      shiny::updateSelectInput(session, "outcome", choices = v,
        selected = if (isTRUE(shiny::isolate(input$outcome) %in% v)) shiny::isolate(input$outcome) else v[1])
    })
    shiny::observe({
      v <- setdiff(c(vars(), "night"), if (input$role == "predictor") input$outcome)
      shinyWidgets::updatePickerInput(session, "covars", choices = v,
        selected = intersect(shiny::isolate(input$covars), v))
    })

    output$level_ui <- shiny::renderUI({
      if (!multilevel()) {
        return(shiny::helpText(shiny::strong("One participant:"),
          "single-level model on the nights (no random effects); substitution",
          "is at the aggregate level."))
      }
      if (input$role != "predictor") return(NULL)
      shiny::radioButtons(session$ns("level"), "Composition enters as",
        c("Between- and within-person ilrs (bz + wz)" = "bw",
          "Total ilrs, not decomposed (z)" = "aggregate"))
    })

    formula_txt <- shiny::reactive({
      k <- length(parts_ordered()) - 1
      cov <- input$covars
      re <- if (multilevel()) "(1 | ID)" else NULL
      if (input$role == "predictor") {
        shiny::req(input$outcome)
        z <- if (multilevel() && !identical(input$level, "aggregate"))
          c(paste0("bz", seq_len(k), "_1"), paste0("wz", seq_len(k), "_1"))
        else paste0("z", seq_len(k), "_1")
        rhs <- c(z, cov, re)
        sprintf("%s ~ %s", input$outcome, paste(rhs, collapse = " + "))
      } else {
        rhs <- c(cov, re)
        if (!length(rhs)) rhs <- "1"
        sprintf("mvbind(%s) ~ %s", paste0("z", seq_len(k), "_1", collapse = ", "),
                paste(rhs, collapse = " + "))
      }
    })
    model_vars <- shiny::reactive(
      c(if (input$role == "predictor") input$outcome, input$covars))

    output$formula <- shiny::renderText(formula_txt())
    output$design_note <- shiny::renderUI({
      a <- ana()
      n_id <- length(unique(a$ID))
      msg <- if (n_id == 1)
        "Single-person model: nights are treated as independent observations. Use it to describe this patient, not to generalise."
      else if (n_id < 10)
        sprintf("Only %d participants: the between-person (bz) estimates rest on %d means and will be very uncertain; the within-person (wz) estimates use all nights.", n_id, n_id)
      else NULL
      if (is.null(msg)) NULL else shiny::div(class = "alert alert-warning", msg)
    })

    fit <- shiny::reactiveVal(NULL)
    fit_warnings <- shiny::reactiveVal(character(0))
    # a new data set or new choices make an old fit stale
    shiny::observeEvent(list(ana(), formula_txt()), {
      fit(NULL)
      rec$clear_from("complr")
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$fit, {
      p <- parts_ordered()
      mv <- model_vars()
      a <- ana()
      keep <- c("ID", "night", p, setdiff(mv, "night"))
      a <- a[stats::complete.cases(a[, c(p, mv), drop = FALSE]), keep, drop = FALSE]
      if (nrow(a) < 5) {
        shiny::showNotification("Too few complete nights to fit a model.", type = "error")
        return()
      }
      idvar <- if (multilevel()) "ID" else NULL
      fam <- input$family
      warn <- character(0)
      res <- withCallingHandlers(
        tryCatch({
          shiny::withProgress(message = "Compiling and sampling with Stan...", value = 0.2, {
            cx <- multilevelcoda::complr(data = a, parts = p,
                    sbp = multilevelcoda::build.sbp(p), idvar = idvar,
                    total = comp$total(), transform = "ilr")
            m <- multilevelcoda::brmcoda(
              complr = cx, formula = stats::as.formula(formula_txt()),
              family = fam,
              chains = input$chains, iter = input$iter, warmup = input$warmup,
              cores = input$cores, seed = input$seed,
              control = list(adapt_delta = input$adapt_delta),
              backend = input$backend, refresh = 0)
            shiny::incProgress(0.8)
            m
          })
        }, error = function(e) e),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        })
      fit_warnings(unique(warn))
      if (inherits(res, "error")) {
        shiny::showNotification(paste("Model failed:", conditionMessage(res)),
                                type = "error", duration = NULL)
        return()
      }
      fit(res)

      rec$record("complr",
        sprintf("Computed the composition and its ilr coordinates (complr; pivot SBP, part order %s; %s)",
                paste(p, collapse = " > "),
                if (is.null(idvar)) "single participant" else "between/within-person split by ID"),
        c(sprintf("model_vars <- %s", r_lit(setdiff(mv, "night"))),
          sprintf("ana_m <- ana[stats::complete.cases(ana[, c(parts, %s)]),\n  c(\"ID\", \"night\", %s, model_vars)]",
                  r_lit(mv), r_lit(p)),
          sprintf("cx <- complr(data = ana_m, parts = %s,\n  sbp = build.sbp(%s), idvar = %s, total = total, transform = \"ilr\")",
                  r_lit(p), r_lit(p), r_lit(idvar)),
          "print(summary(cx))"))
      rec$record("model",
        sprintf("Fitted brmcoda: %s (%s likelihood, %d chains x %d iterations, %d warm-up, adapt_delta %.3f, seed %d, brms default priors)",
                formula_txt(), fam, input$chains, input$iter, input$warmup,
                input$adapt_delta, as.integer(input$seed)),
        c(sprintf("fit <- brmcoda(complr = cx,\n  formula = %s,", formula_txt()),
          sprintf("  family = \"%s\", chains = %d, iter = %d, warmup = %d, cores = %d,",
                  fam, input$chains, input$iter, input$warmup, input$cores),
          sprintf("  seed = %d, control = list(adapt_delta = %s), backend = %s)",
                  as.integer(input$seed), input$adapt_delta, r_lit(input$backend)),
          "print(summary(fit))",
          "print(brms::prior_summary(fit$model))"))
    })

    output$summary <- shiny::renderPrint({
      f <- fit()
      if (is.null(f)) return(cat("Not fitted yet."))
      print(summary(f))
      if (length(fit_warnings())) {
        cat("\nWarnings during fitting:\n")
        cat(paste("-", fit_warnings()), sep = "\n")
      }
    })
    output$diag <- shiny::renderTable({
      f <- fit()
      shiny::req(f)
      m <- f$model
      rh <- brms::rhat(m)
      ess <- brms::neff_ratio(m)
      np <- brms::nuts_params(m)
      div <- sum(np$Value[np$Parameter == "divergent__"])
      data.frame(
        check = c("Max R-hat (want < 1.01)", "Min ESS ratio (want > 0.1)",
                  "Divergent transitions (want 0)", "Post-warm-up draws"),
        value = c(sprintf("%.3f", max(rh, na.rm = TRUE)),
                  sprintf("%.3f", min(ess, na.rm = TRUE)),
                  as.character(div), as.character(brms::ndraws(m))))
    })
    output$ppc <- shiny::renderPlot({
      f <- fit()
      shiny::req(f)
      if (input$role == "predictor") brms::pp_check(f$model, ndraws = 50)
      else brms::pp_check(f$model, ndraws = 50, resp = "z11")
    })

    list(fit = fit, role = shiny::reactive(input$role),
         multilevel = multilevel, parts = parts_ordered)
  })
}
