# ============================================================
# FaRiMaTo: Fatigue Risk Management Tool  v4.0
#
# Implements the McCauley / Van Dongen biomathematical fatigue model.
# THREE published versions are available and selectable in the app. All
# equations and parameter values were transcribed from the printed tables
# of the source papers (Frontiers via its MathML XML; the other two from
# the author-manuscript PDFs), not recalled from memory.
#
#   "2021"     McCauley ME et al. (2021) Transp Res F 79:94-106.
#              doi:10.1016/j.trf.2021.04.006 -- Table 2.
#              No sleep inertia, no reciprocal feedback.
#
#   "2024jtb"  McCauley ME et al. (2024) J Theor Biol 590:111851.
#              doi:10.1016/j.jtbi.2024.111851 -- Tables 1, 2.
#              Adds sleep inertia h(t); nu = 1.37 and gamma = 0.71 estimated
#              there; eta_w = 0.0074 and mu_w = 0.33 held at 2021 values.
#
#   "2024fe"   McCauley ME et al. (2024) Front Environ Health 3:1362755.
#              doi:10.3389/fenvh.2024.1362755 -- Tables 1, 2.
#              Adds reciprocal feedback zeta = 1.31 and REFITS eta_w to
#              0.0126 and mu_w to 0.466 (PVT) / 0.82 (KSS).
#
# The two 2024 papers are NOT interchangeable: they publish different values
# for eta_w and mu_w. Mixing them would be an error, so each version's
# parameters are loaded as a set and the app records which one is in use.
#
# Verification performed (see the Provenance tab):
#   - 2009 Eq. (22) with its Eq. (27) estimates reproduces W_c = 20.2 h, and
#     alpha_22 = 0.00743 / sigma_22 = -0.0397 correspond to eta_w = 0.0074
#     and -eta_s, consistent with the 2021 and 2024 tables.
#     [McCauley P et al. (2009) J Theor Biol 256:227-239]
#   - The 2021 g_s reformulation is algebraically identical to Table 1 of the
#     2024 papers with h(t) = 0 and gamma = 0.
#   - PVT impairment accumulates across days of sleep restriction while KSS
#     levels off after a few days, in all three versions.
#   - The 2024fe version predicts greater build-up across 14 days of 4 h sleep
#     than 2024jtb (PVT day-14 ~20.8 vs ~17.3 lapses), which is the improvement
#     the Frontiers paper reports making.
#   - Impairment stabilises below W_c and escalates above it.
#   - Sleep inertia raises the first post-awakening value in both 2024
#     versions and is absent in the 2021 version.
#
# NOT FROM ANY PAPER, and flagged as such throughout the app:
#   - the composite risk score and its zone cut-points
#   - the VAS-to-KSS rescaling offered for empirical VAS data
#   - the Euler integration step (no paper states an integration scheme)
#
# Earlier versions of this file (FaRiMa v1.x, FaRiMaTo v2.0) used an unsourced
# parameterisation and carried several incorrect citations. Both the equations
# and the values have been replaced.
# ============================================================

# ---- 0. Package management (House Rule 0) ----------------------------
required_pkgs <- c(
  "shiny", "shinydashboard", "shinyWidgets", "DT",
  "plotly", "ggplot2", "dplyr", "tidyr",
  "readr", "readxl", "haven", "lubridate", "RColorBrewer"
)
miss <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(miss) > 0) {
  message("Installing: ", paste(miss, collapse = ", "))
  install.packages(miss, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(shiny); library(shinydashboard); library(shinyWidgets)
  library(DT); library(plotly); library(ggplot2)
  library(dplyr); library(tidyr)
  library(readr); library(readxl); library(haven)
  library(lubridate); library(RColorBrewer)
})
pkg_versions <- sapply(required_pkgs, function(p)
  tryCatch(as.character(packageVersion(p)), error = function(e) "N/A"))


# ============================================================
# 1. STEP RECORDER  (House Rules 6 + 8)
# ============================================================
new_recorder <- function() {
  steps <- list()
  list(
    record = function(id, description, code) {
      steps[[id]] <<- list(id = id, description = description,
                           code = code, ts = Sys.time())
    },
    build_script = function() {
      header <- paste(c(
        "# FaRiMa -- Reproducible Fatigue Risk Analysis",
        paste("# Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        "#",
        "# Package versions:",
        paste0("#  ", names(pkg_versions), ": ", pkg_versions),
        "",
        "library(dplyr); library(ggplot2); library(lubridate)",
        ""
      ), collapse = "\n")
      code_blocks <- paste(sapply(steps, `[[`, "code"), collapse = "\n\n")
      paste0(header, "\n", code_blocks)
    },
    build_log = function() {
      if (length(steps) == 0) return("No steps recorded yet.")
      paste(mapply(function(i, s) sprintf(
        "[%s] Step %d (%s): %s",
        format(s$ts, "%H:%M:%S"), i, s$id, s$description),
        seq_along(steps), steps), collapse = "\n")
    }
  )
}


# ============================================================
# 2. MODEL CORE  --  McCauley et al. (2024) biomathematical model
# ============================================================
#
# SOURCE (verified against the published article, equations and
# parameter values transcribed from the journal's NLM XML):
#
#   McCauley ME, McCauley P, Kalachev LV, Banks S, Dinges DF,
#   Van Dongen HPA (2024). The dynamics of neurobehavioral impairment
#   and recovery sleep: improved biomathematical modeling for fatigue
#   risk management in operational settings.
#   Front. Environ. Health 3:1362755. doi:10.3389/fenvh.2024.1362755
#   Open access (CC BY 4.0).
#
#   Equations  : Table 1 (T1-1 ... T1-6), plus Equations 4 and 9.
#   Parameters : Table 2.
#
# STATE VARIABLES
#   p(t)      neurobehavioral impairment during wakefulness (or the
#             nominal propensity thereof during sleep). This IS the
#             predicted outcome: PVT lapses, or KSS score.
#   u(t)      allostatic process (days to weeks)
#   h(t)      sleep inertia
#   kappa(t)  time-varying amplitude of the circadian impact
#
# The PVT/KSS distinction is handled INSIDE the equations through
# metric-specific p_f, alpha_w, alpha_s, W_c, mu_w and xi. There is
# no link function and no separate "output arm": the model is run
# once per metric with that metric's parameter set.

#' Published parameter values for the three verified model versions.
#'
#' All three share the same equation family; they differ in which terms are
#' active and in the values of eta_w and mu_w. Verified against the printed
#' tables of the papers named below.
#'
#' "2021"    McCauley ME, McCauley P, Riedy SM, Banks S, Ecker AJ, Kalachev LV,
#'           Rangan S, Dinges DF, Van Dongen HPA (2021). Fatigue risk management
#'           based on self-reported fatigue: Expanding a biomathematical model
#'           of fatigue-related performance deficits to also predict subjective
#'           sleepiness. Transp Res F 79:94-106. doi:10.1016/j.trf.2021.04.006
#'           -- Table 2. No sleep inertia, no reciprocal feedback.
#'
#' "2024jtb" McCauley ME, McCauley P, Kalachev LV, Riedy SM, Banks S, Ecker AJ,
#'           Dinges DF, Van Dongen HPA (2024). Biomathematical modeling of
#'           fatigue due to sleep inertia. J Theor Biol 590:111851.
#'           doi:10.1016/j.jtbi.2024.111851 -- Tables 1 and 2. Adds sleep
#'           inertia h(t); nu and gamma newly estimated there; all other
#'           parameters held at their 2021 values.
#'
#' "2024fe"  McCauley ME, McCauley P, Kalachev LV, Banks S, Dinges DF,
#'           Van Dongen HPA (2024). The dynamics of neurobehavioral impairment
#'           and recovery sleep. Front Environ Health 3:1362755.
#'           doi:10.3389/fenvh.2024.1362755 -- Tables 1 and 2. Adds the
#'           reciprocal feedback parameter zeta and REFITS eta_w and mu_w.
#'
#' Lineage check performed against McCauley P, Kalachev LV, Smith AD, Belenky G,
#' Dinges DF, Van Dongen HPA (2009), J Theor Biol 256:227-239, Eq. (22) and
#' Eq. (27): alpha_22 = 0.00743 and sigma_22 = -0.0397 reproduce W_c = 20.2 h
#' and correspond to eta_w = 0.0074 and -eta_s, consistent with the later papers.
#'
#' @param metric  "PVT" or "KSS"
#' @param version "2021", "2024jtb" or "2024fe"
mccauley_params <- function(metric = c("PVT", "KSS"),
                            version = c("2024fe", "2024jtb", "2021")) {
  metric  <- match.arg(metric)
  version <- match.arg(version)
  pvt     <- identical(metric, "PVT")

  p <- list(
    metric = metric, version = version,
    # --- shared across all three versions (2021 Table 2 / JTB 2024 Table 2) ---
    p_f     = if (pvt) 0      else 1,      # floor of metric scale
    alpha_w = if (pvt) 0.028  else 0.22,   # h^-1  homeostatic build-up
    alpha_s = if (pvt) 0.26   else 0.037,  # h^-1  homeostatic dissipation
    W_c     = if (pvt) 20.2   else 22.02,  # h     bifurcation threshold
    xi      = if (pvt) 1.09   else 0.51,   #       metric scaling factor
    beta    = 0.26,                        #       impact of u onto p
    mu_s    = -1.5,                        #       circadian offset, sleep
    phi     = 21.2,                        # h     circadian phase vs midnight
    lambda  = 0.49,                        # h^-1  circadian amplitude rate
    tau     = 24.0, T_day = 24.0
  )

  # --- version-specific -------------------------------------------------
  if (version == "2024fe") {
    p$eta_w <- 0.0126                      # refitted in Frontiers 2024
    p$mu_w  <- if (pvt) 0.466 else 0.82    # refitted, metric-specific
    p$zeta  <- 1.31                        # reciprocal feedback, new in 2024fe
    p$nu    <- 1.37; p$gamma <- 0.71       # sleep inertia
    p$use_inertia <- TRUE
  } else if (version == "2024jtb") {
    p$eta_w <- 0.0074                      # 2021 value, NOT refitted
    p$mu_w  <- 0.33                        # 2021 value, same for both metrics
    p$zeta  <- 0                           # no reciprocal feedback
    p$nu    <- 1.37; p$gamma <- 0.71       # estimated in JTB 2024
    p$use_inertia <- TRUE
  } else {                                  # "2021"
    p$eta_w <- 0.0074
    p$mu_w  <- 0.33
    p$zeta  <- 0
    # gamma MUST be 0, not merely unused. gamma appears in the sleep equation
    # for p as the factor (1 - nu*gamma*xi/alpha_s); leaving gamma at 0.71 with
    # h held at zero would flip the sign of that coefficient and make p diverge
    # during sleep. Setting gamma = 0 reduces Table 1 exactly to the 2021 Eqs.
    p$nu    <- 1.37; p$gamma <- 0
    p$use_inertia <- FALSE                 # no sleep inertia in the 2021 model
  }

  # Equation 9: eta_s = eta_w * W_c / (T - W_c)
  p$eta_s <- p$eta_w * p$W_c / (p$T_day - p$W_c)
  p
}

# Back-compatible alias
mccauley2024_params <- function(metric = c("PVT", "KSS"))
  mccauley_params(metric, version = "2024fe")

model_versions <- c(
  "Frontiers 2024 -- inertia + reciprocal feedback (zeta)" = "2024fe",
  "J Theor Biol 2024 -- inertia, no reciprocal feedback"   = "2024jtb",
  "Transp Res F 2021 -- no sleep inertia"                  = "2021"
)

#' Circadian process, Equation (T1-6): c(t) = sin(2*pi*(t - phi)/tau)
#' @param t_clock clock time in hours
circadian_c <- function(t_clock, phi, tau) sin(2 * pi * (t_clock - phi) / tau)

#' One Euler step of the McCauley et al. (2024) system.
#' Implements Table 1 equations T1-1 through T1-5 exactly as printed.
#' @param st   list(p, u, h, kappa)
#' @param wake TRUE during wakefulness, FALSE during sleep
mc2024_step <- function(st, wake, t_clock, P, dt) {
  p <- st$p; u <- st$u; h <- st$h; k <- st$kappa
  cc <- circadian_c(t_clock, P$phi, P$tau)
  pf <- P$p_f
  xi <- P$xi

  # The 2021 model has no sleep-inertia process: hold h at zero throughout,
  # which reduces Table 1 exactly to the 2021 equations (verified algebraically
  # against that paper's Eqs. 1b and 5b).
  if (!isTRUE(P$use_inertia)) h <- 0

  if (wake) {
    # (T1-5a) g_w(t) = kappa(t) [ c(t) + mu_w ]
    g  <- k * (cc + P$mu_w)
    # (T1-1a) dp/dt = -alpha_w [ p - p_f - beta_w u ]
    #                 + xi_h (alpha_w - nu_w) h + xi_kappa g_w
    dp <- -P$alpha_w * (p - pf - P$beta * u) +
            xi * (P$alpha_w - P$nu) * h + xi * g
    # (T1-2a) du/dt = -eta_w [ zeta_w xi_u (p - p_f) - u ]
    #                 + eta_w zeta_w xi_u xi_h h
    du <- -P$eta_w * (P$zeta * xi * (p - pf) - u) +
            P$eta_w * P$zeta * xi * xi * h
    # (T1-3a) dh/dt = -nu_w h
    dh <- -P$nu * h
    # (T1-4a) dkappa/dt = lambda_w kappa (1 - kappa)
    dk <- P$lambda * k * (1 - k)
  } else {
    # (T1-5b) g_s(t) = kappa(t) [ c(t) + mu_s ]
    g  <- k * (cc + P$mu_s)
    # (T1-1b) dp/dt = -alpha_s [ (1 - nu_s gamma xi_h / alpha_s)(p - p_f)
    #                            - beta_s (u - 1/eta_s) ]
    #                 + xi_h (alpha_s - nu_s (1 + gamma xi_h)) h + xi_kappa g_s
    dp <- -P$alpha_s * ((1 - P$nu * P$gamma * xi / P$alpha_s) * (p - pf) -
                          P$beta * (u - 1 / P$eta_s)) +
            xi * (P$alpha_s - P$nu * (1 + P$gamma * xi)) * h + xi * g
    # (T1-2b) du/dt = eta_s [ zeta_s xi_u (p - p_f) - u ]
    #                 - eta_s zeta_s xi_u xi_h h + 1
    du <- P$eta_s * (P$zeta * xi * (p - pf) - u) -
            P$eta_s * P$zeta * xi * xi * h + 1
    # (T1-3b) dh/dt = -nu_s [ (1 + gamma xi_h) h - gamma (p - p_f) ]
    dh <- -P$nu * ((1 + P$gamma * xi) * h - P$gamma * (p - pf))
    # (T1-4b) dkappa/dt = -lambda_s kappa
    dk <- -P$lambda * k
  }

  if (!isTRUE(P$use_inertia)) dh <- 0

  list(p     = p + dt * dp,
       u     = u + dt * du,
       h     = h + dt * dh,
       kappa = min(max(k + dt * dk, 1e-9), 1))
}

#' Equilibrate initial conditions on a repeating baseline day.
#'
#' The paper states: "For determining the initial values of p, u, h and
#' kappa, we assumed an equilibrium state for the sleep/wake schedule of
#' the baseline day(s) in each study, and iterated across days using that
#' schedule until a steady state over days at the scheduled time of
#' awakening was reached (at a tolerance level of eps = 1e-9)."
#'
#' @param base_onset,base_offset baseline sleep window in clock hours
#' @return list(p, u, h, kappa) at the baseline awakening time
mc2024_equilibrate <- function(P, base_onset = 23, base_offset = 7,
                               dt = 0.05, tol = 1e-9, max_days = 400) {

  # Bifurcation guard. The model has a critical daily wake duration W_c:
  # below it impairment stabilises over days, above it impairment escalates
  # without bound, so no steady state exists and equilibration would diverge.
  sleep_dur <- (base_offset - base_onset) %% 24
  wake_dur  <- P$T_day - sleep_dur
  if (wake_dur >= P$W_c) {
    stop(sprintf(paste0(
      "Baseline day has %.2f h of wakefulness, at or above the critical ",
      "threshold W_c = %.2f h for the %s arm. Above W_c the model escalates ",
      "without a steady state, so initial conditions cannot be equilibrated ",
      "on this baseline. Use a baseline with more sleep, or switch off ",
      "equilibration."), wake_dur, P$W_c, P$metric))
  }

  st <- list(p = P$p_f, u = 0.5, h = 0, kappa = 0.5)

  # Build a 24 h wake indicator starting at the baseline awakening hour
  n      <- round(P$T_day / dt)
  offs   <- seq(0, by = dt, length.out = n)
  clock  <- (base_offset + offs) %% 24
  # asleep between onset and offset (handles wrap past midnight)
  if (base_onset < base_offset) {
    asleep <- clock >= base_onset & clock < base_offset
  } else {
    asleep <- clock >= base_onset | clock < base_offset
  }

  for (d in seq_len(max_days)) {
    prev <- st
    for (i in seq_len(n)) {
      st <- mc2024_step(st, wake = !asleep[i], t_clock = clock[i], P = P, dt = dt)
    }
    if (!all(is.finite(unlist(st))))
      stop("Equilibration diverged; check the baseline schedule and parameters.")
    delta <- max(abs(c(st$p - prev$p, st$u - prev$u,
                       st$h - prev$h, st$kappa - prev$kappa)))
    if (is.finite(delta) && delta < tol) break
  }
  st
}

#' Integrate the model over a sleep/wake schedule for one metric.
#' @param schedule_h data.frame(onset_h, offset_h) in hours from t = 0
#' @return data.frame(time_h, clock_h, wake, p, u, h, kappa, c)
simulate_mc2024 <- function(schedule_h, P, dt = 0.05, t0_clock = 8,
                            base_onset = 23, base_offset = 7,
                            equilibrate = TRUE) {
  if (nrow(schedule_h) == 0) stop("Schedule must contain at least one sleep period.")
  schedule_h <- schedule_h[schedule_h$offset_h > schedule_h$onset_h, , drop = FALSE]
  if (nrow(schedule_h) == 0) stop("No sleep period has offset later than onset.")

  t_end <- max(schedule_h$offset_h) + 8
  times <- seq(0, t_end, by = dt)
  n     <- length(times)
  clock <- (t0_clock + times) %% 24

  W <- rep(TRUE, n)                       # TRUE = awake
  for (j in seq_len(nrow(schedule_h))) {
    W[times >= schedule_h$onset_h[j] & times < schedule_h$offset_h[j]] <- FALSE
  }

  st <- if (equilibrate) {
    mc2024_equilibrate(P, base_onset, base_offset, dt = max(dt, 0.05))
  } else {
    list(p = P$p_f, u = 0.5, h = 0, kappa = 0.5)
  }

  pv <- uv <- hv <- kv <- numeric(n)
  pv[1] <- st$p; uv[1] <- st$u; hv[1] <- st$h; kv[1] <- st$kappa
  for (i in 2:n) {
    st <- mc2024_step(st, wake = W[i], t_clock = clock[i], P = P, dt = dt)
    pv[i] <- st$p; uv[i] <- st$u; hv[i] <- st$h; kv[i] <- st$kappa
  }

  data.frame(time_h = times, clock_h = clock, wake = as.integer(W),
             p = pv, u = uv, h = hv, kappa = kv,
             c = circadian_c(clock, P$phi, P$tau))
}

#' Run both metric arms and assemble the results table.
#'
#' The risk score appended here is NOT from the paper -- see the
#' Provenance tab. It is a weighted combination of the two published
#' predictions, with weights and cut-points chosen by the app author.
simulate_fatigue <- function(schedule_h, dt = 0.05, t0_clock = 8,
                             base_onset = 23, base_offset = 7,
                             equilibrate = TRUE,
                             P_pvt = mccauley_params("PVT"),
                             P_kss = mccauley_params("KSS"),
                             risk_w = c(0.5, 0.5), pvt_ceiling = 20) {

  rp <- simulate_mc2024(schedule_h, P_pvt, dt, t0_clock,
                        base_onset, base_offset, equilibrate)
  rk <- simulate_mc2024(schedule_h, P_kss, dt, t0_clock,
                        base_onset, base_offset, equilibrate)

  out <- rp[, c("time_h", "clock_h", "wake", "c")]
  out$PVT     <- pmax(rp$p, 0)                    # lapses cannot be negative
  out$KSS     <- pmin(pmax(rk$p, 1), 9)           # KSS is bounded 1-9
  out$PVT_raw <- rp$p
  out$KSS_raw <- rk$p
  out$u_PVT   <- rp$u;  out$u_KSS   <- rk$u
  out$h_PVT   <- rp$h;  out$h_KSS   <- rk$h
  out$kappa   <- rp$kappa

  w <- risk_w / sum(risk_w)
  out$Risk <- pmin(pmax(
    w[1] * (out$KSS - 1) / 8 +
    w[2] * pmin(out$PVT / pvt_ceiling, 1), 0), 1)
  out
}

# Helper: detect sleep-period rectangles for plotly shapes
sleep_shapes <- function(res, y0 = 0, y1 = 1) {
  runs   <- rle(res$wake)
  idx    <- 0
  shapes <- list()
  for (i in seq_along(runs$values)) {
    s_i <- idx + 1
    e_i <- idx + runs$lengths[i]
    if (runs$values[i] == 0L) {
      shapes[[length(shapes) + 1]] <- list(
        type = "rect", layer = "below",
        fillcolor = "rgba(174, 214, 241, 0.35)",
        line = list(width = 0),
        x0 = res$time_h[s_i], x1 = res$time_h[min(e_i, nrow(res))],
        y0 = y0, y1 = y1, yref = "paper")
    }
    idx <- e_i
  }
  shapes
}

# Helper: convert #RRGGBB to rgba(...) string for plotly
hex_to_rgba <- function(hex, alpha = 0.25) {
  r <- strtoi(substr(hex, 2L, 3L), 16L)
  g <- strtoi(substr(hex, 4L, 5L), 16L)
  b <- strtoi(substr(hex, 6L, 7L), 16L)
  sprintf("rgba(%d,%d,%d,%.2f)", r, g, b, alpha)
}

# Helper: build a distinct colour per group (Set2 pastel)
group_colours <- function(groups) {
  n <- length(groups)
  if (n <= 1) return(setNames("#e74c3c", groups))
  setNames(RColorBrewer::brewer.pal(max(min(n, 8), 3), "Set2")[seq_len(n)], groups)
}

# Published parameter tables for all three verified versions, for the
# Provenance tab. Transcribed from the printed tables of each paper.
published_reference_params <- local({
  mk <- function(v, lbl) {
    Pp <- mccauley_params("PVT", v); Pk <- mccauley_params("KSS", v)
    nm <- c("p_f","alpha_w","alpha_s","W_c","xi","beta","eta_w","eta_s",
            "zeta","nu","gamma","mu_w","mu_s","phi","lambda","tau")
    data.frame(
      Version   = lbl,
      Parameter = nm,
      PVT = vapply(nm, function(x) {
              if (!Pp$use_inertia && x %in% c("nu","gamma")) "n/a (no inertia)"
              else if (x == "zeta" && Pp$zeta == 0) "0 (absent)"
              else format(signif(Pp[[x]], 6))}, character(1)),
      KSS = vapply(nm, function(x) {
              if (!Pk$use_inertia && x %in% c("nu","gamma")) "n/a (no inertia)"
              else if (x == "zeta" && Pk$zeta == 0) "0 (absent)"
              else format(signif(Pk[[x]], 6))}, character(1)),
      row.names = NULL, stringsAsFactors = FALSE)
  }
  rbind(mk("2024fe",  "Frontiers 2024"),
        mk("2024jtb", "J Theor Biol 2024"),
        mk("2021",    "Transp Res F 2021"))
})

# ============================================================
# 2a. COLUMN-NAME TIME PARSING
# ============================================================
# Wide-format questionnaire exports usually encode the measurement moment in
# the column name itself. This parser recognises the common Dutch/English/
# numeric conventions and proposes a "hours from t=0" vector, so the user does
# not have to type the time grid by hand.
#
# Recognised clock forms : 9  9h  9u  9:00  9.30  9h00  9u30  0900  21:45
# Recognised day markers : dag1  dag_2  day3  d2  _t1 (treated as day index)
# Pure occasion indices  : kss1, vas3, pvt_2  -> no clock information

#' Parse one column name into a clock time and (optional) day index
#' @return list(kind = "clock"|"index"|"none", clock = numeric h, day = numeric)
parse_time_token <- function(nm) {
  s <- tolower(trimws(as.character(nm)))
  s <- gsub("[[:space:]]+", "", s)

  out <- list(kind = "none", clock = NA_real_, day = NA_real_)

  # --- 1. Day marker: dag1 / dag_2 / day3 / d2 (d only when followed by digits
  #        and not part of a longer word) ---
  day_pat <- "(dag|day)[_\\-\\.]?([0-9]+)"
  m <- regexec(day_pat, s, perl = TRUE)
  g <- regmatches(s, m)[[1]]
  if (length(g) >= 3) {
    out$day <- as.numeric(g[3])
    s <- sub(day_pat, "", s, perl = TRUE)
  } else {
    m2 <- regexec("(^|[_\\-\\.])d([0-9]+)([_\\-\\.]|$)", s, perl = TRUE)
    g2 <- regmatches(s, m2)[[1]]
    if (length(g2) >= 3) {
      out$day <- as.numeric(g2[3])
      s <- sub("(^|[_\\-\\.])d([0-9]+)([_\\-\\.]|$)", "\\1", s, perl = TRUE)
    }
  }

  # --- 2. Occasion-index form: letters followed by digits (kss1, vas_3) ---
  # Checked BEFORE bare numbers so "kss1" is not misread as 01:00.
  if (grepl("^[a-z]+[_\\-]?[0-9]+[_\\-]?$", s, perl = TRUE)) {
    idx <- as.numeric(gsub("[^0-9]", "", s))
    out$kind <- "index"; out$clock <- idx
    return(out)
  }

  # --- 3. hh:mm / hh.mm / hhhmm / hhumm  (explicit separator + 2-digit min) ---
  m3 <- regexec("([0-9]{1,2})[:hu\\.]([0-9]{2})", s, perl = TRUE)
  g3 <- regmatches(s, m3)[[1]]
  if (length(g3) >= 3) {
    hh <- as.numeric(g3[2]); mm <- as.numeric(g3[3])
    if (hh <= 23 && mm <= 59) {
      out$kind <- "clock"; out$clock <- hh + mm / 60
      return(out)
    }
  }

  # --- 4. hh followed by a bare h/u marker: 9h, 21u ---
  m4 <- regexec("([0-9]{1,2})[hu]([^0-9]|$)", s, perl = TRUE)
  g4 <- regmatches(s, m4)[[1]]
  if (length(g4) >= 2) {
    hh <- as.numeric(g4[2])
    if (hh <= 23) { out$kind <- "clock"; out$clock <- hh; return(out) }
  }

  # --- 5. Compact 4-digit military time: 0900, 2145 ---
  m5 <- regexec("(^|[^0-9])([0-9]{4})([^0-9]|$)", s, perl = TRUE)
  g5 <- regmatches(s, m5)[[1]]
  if (length(g5) >= 3) {
    v  <- g5[3]
    hh <- as.numeric(substr(v, 1, 2)); mm <- as.numeric(substr(v, 3, 4))
    if (hh <= 23 && mm <= 59) {
      out$kind <- "clock"; out$clock <- hh + mm / 60
      return(out)
    }
  }

  # --- 6. Bare 1-2 digit number: treat as whole clock hour ---
  m6 <- regexec("(^|[^0-9])([0-9]{1,2})([^0-9]|$)", s, perl = TRUE)
  g6 <- regmatches(s, m6)[[1]]
  if (length(g6) >= 3) {
    hh <- as.numeric(g6[3])
    if (hh <= 23) { out$kind <- "clock"; out$clock <- hh; return(out) }
  }

  out
}

#' Propose a "hours from t=0" vector for a set of column names
#'
#' Clock times are referenced to `t0_clock` (the study reference hour). Day
#' markers, when present, resolve multi-day designs directly. Without day
#' markers the series is "unwrapped" across midnight on the assumption that
#' measurement occasions are in chronological order -- the standard layout of
#' a wide-format export.
#'
#' @return list(times = numeric vector, method = character, detail = character)
propose_times_from_names <- function(col_names, t0_clock = 8, fallback_step = 8) {
  n <- length(col_names)
  if (n == 0) return(list(times = numeric(0), method = "none", detail = "No columns selected."))

  parsed <- lapply(col_names, parse_time_token)
  kinds  <- vapply(parsed, function(x) x$kind, character(1))
  clocks <- vapply(parsed, function(x) x$clock, numeric(1))
  days   <- vapply(parsed, function(x) x$day,   numeric(1))

  n_clock <- sum(kinds == "clock")

  # --- Case A: clock times recovered from the majority of column names ---
  if (n_clock >= max(2, ceiling(0.6 * n))) {

    # Interpolate any column whose clock could not be read
    if (anyNA(clocks[kinds == "clock"]) || n_clock < n) {
      idx_ok <- which(kinds == "clock")
      clocks <- approx(idx_ok, clocks[idx_ok], xout = seq_len(n), rule = 2)$y
    }

    if (sum(!is.na(days)) >= max(2, ceiling(0.6 * n))) {
      # Day markers present: fill gaps, then resolve absolute hours
      idx_d <- which(!is.na(days))
      days  <- approx(idx_d, days[idx_d], xout = seq_len(n), rule = 2)$y
      days  <- round(days)
      raw   <- (days - min(days, na.rm = TRUE)) * 24 + clocks - t0_clock
      if (min(raw, na.rm = TRUE) < 0) raw <- raw + 24
      return(list(times = round(raw, 3), method = "clock+day",
                  detail = sprintf("Clock times and day markers detected; referenced to %.2f h.", t0_clock)))
    }

    # No day markers: unwrap the clock series across midnight
    raw <- clocks - t0_clock
    if (!is.na(raw[1]) && raw[1] < 0) raw[1] <- raw[1] + 24
    if (n > 1) {
      for (i in 2:n) {
        while (!is.na(raw[i]) && !is.na(raw[i - 1]) && raw[i] < raw[i - 1]) raw[i] <- raw[i] + 24
      }
    }
    return(list(times = round(raw, 3), method = "clock",
                detail = sprintf("Clock times detected in column names; referenced to %.2f h and unwrapped across midnight.", t0_clock)))
  }

  # --- Case B: occasion indices only -> evenly spaced fallback grid ---
  if (all(kinds == "index") && n >= 1) {
    return(list(times = seq(0, by = fallback_step, length.out = n), method = "index",
                detail = sprintf("Column names carry occasion numbers only (no clock time); proposing an even %g h spacing -- edit if your protocol differs.", fallback_step)))
  }

  list(times = seq(0, by = fallback_step, length.out = n), method = "fallback",
       detail = sprintf("No time information recognised in column names; proposing an even %g h spacing -- please edit.", fallback_step))
}


# ============================================================
# 2b. SHARED HELPERS: empirical overlay drawing + comparison outputs
# ============================================================

#' Overlay one empirical outcome (individual traces and/or mean +/- SD band)
#' onto an existing plotly object. Reused for KSS, VAS-derived KSS-equivalent,
#' and PVT overlays so the drawing logic is written once.
#'
#' @param p            existing plotly object
#' @param emp_df       long-format data.frame with columns: subject, group,
#'                     time_h, and `value_col`
#' @param value_col    name of the column holding the numeric outcome
#' @param yaxis        "y" or "y2"
#' @param y_clip       optional c(lo, hi) to clip the +/-SD ribbon
add_outcome_overlay <- function(p, emp_df, value_col, yaxis = "y",
                                show_indiv = FALSE, show_mean = TRUE,
                                label_prefix = "Observed",
                                line_dash = "solid", marker_symbol = "circle",
                                y_clip = NULL) {
  if (is.null(emp_df) || nrow(emp_df) == 0 || !(show_indiv || show_mean)) return(p)

  groups   <- sort(unique(emp_df$group))
  n_groups <- length(groups)
  grp_cols <- group_colours(groups)

  if (show_indiv) {
    for (grp in groups) {
      d_g <- emp_df[emp_df$group == grp, ]
      x_v <- y_v <- c()
      for (sid in unique(d_g$subject)) {
        d_s <- d_g[d_g$subject == sid, ]
        d_s <- d_s[order(d_s$time_h), ]
        x_v <- c(x_v, d_s$time_h, NA_real_)
        y_v <- c(y_v, d_s[[value_col]], NA_real_)
      }
      lbl <- if (n_groups <= 1) paste0(label_prefix, " (individual)")
             else paste0(label_prefix, ": ", grp, " (individual)")
      p <- add_trace(p, x = x_v, y = y_v, type = "scatter", mode = "lines",
                     yaxis = yaxis,
                     line = list(color = hex_to_rgba(unname(grp_cols[[grp]]), 0.30),
                                 width = 0.9, dash = line_dash),
                     name = lbl, legendgroup = paste0("indiv_", label_prefix, "_", grp),
                     showlegend = TRUE, hoverinfo = "none")
    }
  }

  if (show_mean) {
    for (grp in groups) {
      d_grp <- emp_df[emp_df$group == grp, ] %>%
        dplyr::group_by(time_h) %>%
        dplyr::summarise(m = mean(.data[[value_col]], na.rm = TRUE),
                         s = sd(.data[[value_col]],   na.rm = TRUE),
                         n = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(time_h)

      col <- unname(grp_cols[[grp]])
      lbl <- if (n_groups <= 1) paste0(label_prefix, " (mean +/- SD)")
             else paste0(label_prefix, ": ", grp, " (mean +/- SD)")

      upper <- d_grp$m + d_grp$s
      lower <- d_grp$m - d_grp$s
      if (!is.null(y_clip)) {
        upper <- pmin(upper, y_clip[2])
        lower <- pmax(lower, y_clip[1])
      }

      x_rib <- c(d_grp$time_h, rev(d_grp$time_h))
      y_rib <- c(upper, rev(lower))
      p <- add_trace(p, x = x_rib, y = y_rib, type = "scatter", mode = "none",
                     fill = "toself", fillcolor = hex_to_rgba(col, 0.15),
                     line = list(color = "transparent"), yaxis = yaxis,
                     legendgroup = paste0("mean_", label_prefix, "_", grp),
                     name = paste0(lbl, " [ribbon]"), showlegend = FALSE,
                     hoverinfo = "none")

      p <- add_trace(p, data = d_grp, x = ~time_h, y = ~m,
                     type = "scatter", mode = "lines+markers", yaxis = yaxis,
                     legendgroup = paste0("mean_", label_prefix, "_", grp),
                     name = lbl,
                     line   = list(color = col, width = 2, dash = line_dash),
                     marker = list(color = col, size = 7, symbol = marker_symbol,
                                   line = list(color = col, width = 2)),
                     customdata  = ~n,
                     hovertemplate = paste0("<b>t = %{x:.1f} h</b><br>",
                                            label_prefix, " = %{y:.2f}<br>",
                                            "n = %{customdata}<extra></extra>"))
    }
  }
  p
}

#' Generic wide-format empirical outcome importer (House Rules 1, 3, 6, 8)
#'
#' Wires up: file upload -> column auto-detect -> column/ID/group pickers
#' -> time-mapping text box -> long-format reshape. Call once per outcome
#' (KSS / VAS / PVT) with a unique `prefix`; every input/output ID used is
#' `prefix`-qualified so multiple instances can coexist without collision.
#'
#' @return reactive() yielding a long-format data.frame with columns
#'         subject, group, occasion, time_h, value  (value = RAW value,
#'         no unit/scale transform applied -- that is done by the caller)
make_wide_outcome_importer <- function(input, output, session, rec, rv,
                                       prefix, label, col_regex,
                                       default_step_h = 8) {

  file_id   <- paste0(prefix, "_file")
  raw_slot  <- paste0(prefix, "_raw")
  auto_slot <- paste0(prefix, "_cols_auto")

  # ---- Parse uploaded file -----------------------------------
  observeEvent(input[[file_id]], {
    f <- input[[file_id]]
    req(f)
    ext <- tolower(tools::file_ext(f$name))

    df <- tryCatch({
      if (ext %in% c("xlsx", "xls")) {
        readxl::read_excel(f$datapath)
      } else {
        raw1  <- readLines(f$datapath, n = 1L, encoding = "UTF-8", warn = FALSE)
        delim <- if (grepl(";", raw1, fixed = TRUE)) ";" else ","
        readr::read_delim(f$datapath, delim = delim,
                          locale         = readr::locale(encoding = "UTF-8"),
                          show_col_types = FALSE, trim_ws = TRUE)
      }
    }, error = function(e) {
      showNotification(paste0(label, " file error: ", e$message), type = "error")
      NULL
    })

    if (!is.null(df)) {
      names(df) <- gsub("^\uFEFF", "", names(df))  # strip UTF-8 BOM
      auto_cols <- grep(col_regex, names(df), ignore.case = TRUE, value = TRUE)
      rv[[raw_slot]]  <- df
      rv[[auto_slot]] <- auto_cols
      showNotification(
        sprintf("Loaded %s file: %d subjects, %d column(s) auto-detected.",
                label, nrow(df), length(auto_cols)),
        type = "message", duration = 5L
      )
    }
  })

  output[[paste0(prefix, "_loaded")]] <- reactive({ !is.null(rv[[raw_slot]]) })
  outputOptions(output, paste0(prefix, "_loaded"), suspendWhenHidden = FALSE)

  # --- Dynamic UI: column picker (multi-select with search) ---
  output[[paste0(prefix, "_ui_cols")]] <- renderUI({
    req(!is.null(rv[[raw_slot]]))
    pickerInput(paste0(prefix, "_cols_selected"),
                label    = paste0(label, " time-point columns:"),
                choices  = names(rv[[raw_slot]]),
                selected = rv[[auto_slot]],
                multiple = TRUE,
                options  = list(`actions-box`          = TRUE,
                                `live-search`          = TRUE,
                                `selected-text-format` = "count > 5",
                                `count-selected-text`  = "{0} columns selected"))
  })

  # --- Dynamic UI: subject ID column ---
  output[[paste0(prefix, "_ui_id")]] <- renderUI({
    req(!is.null(rv[[raw_slot]]))
    selectInput(paste0(prefix, "_id_col"), "Subject ID column:",
                choices  = c("(row number)" = "__rownum__", names(rv[[raw_slot]])),
                selected = names(rv[[raw_slot]])[1])
  })

  # --- Dynamic UI: optional grouping variable ---
  output[[paste0(prefix, "_ui_group")]] <- renderUI({
    req(!is.null(rv[[raw_slot]]))
    df <- rv[[raw_slot]]
    cat_cols <- Filter(function(cn) {
      v <- df[[cn]]
      n_uniq <- length(unique(v[!is.na(v)]))
      n_uniq >= 2L && n_uniq <= 10L && !cn %in% rv[[auto_slot]]
    }, names(df))
    selectInput(paste0(prefix, "_group_col"), "Group variable (optional):",
                choices  = c("None" = "__none__", cat_cols),
                selected = "__none__")
  })

  # --- Time proposal derived from the column names themselves ---------
  time_proposal <- reactive({
    sel  <- input[[paste0(prefix, "_cols_selected")]]
    cols <- if (!is.null(sel) && length(sel) > 0) sel else rv[[auto_slot]]
    t0   <- suppressWarnings(as.numeric(input$t0_clock))
    if (!is.finite(t0)) t0 <- 8
    propose_times_from_names(cols, t0_clock = t0, fallback_step = default_step_h)
  })

  # --- Dynamic UI: measurement time mapper ---
  # Structure re-renders when the column selection changes; the *values* are
  # additionally refreshed by the observer below whenever the study reference
  # hour changes, so the grid always tracks t0_clock.
  output[[paste0(prefix, "_ui_times")]] <- renderUI({
    req(!is.null(rv[[raw_slot]]))
    prop   <- time_proposal()
    n_cols <- length(prop$times)
    tagList(
      textInput(paste0(prefix, "_times_txt"),
                label = sprintf("Measurement times (h from t=0), %d values, same order as columns:", n_cols),
                value = paste(prop$times, collapse = ", "),
                placeholder = "e.g.  0, 8, 16, 24, 32"),
      helpText(HTML(sprintf("<b>Auto-detected (%s):</b> %s", prop$method, prop$detail))),
      actionButton(paste0(prefix, "_btn_redetect"), "Re-detect from column names",
                   class = "btn-xs btn-default", icon = icon("sync"))
    )
  })

  # Refresh proposed times when the study reference hour changes, or on demand.
  observeEvent(list(input$t0_clock, input[[paste0(prefix, "_btn_redetect")]]), {
    req(!is.null(rv[[raw_slot]]))
    prop <- time_proposal()
    req(length(prop$times) > 0)
    updateTextInput(session, paste0(prefix, "_times_txt"),
                    value = paste(prop$times, collapse = ", "))
  }, ignoreInit = TRUE)

  # --- Summary text ---
  output[[paste0(prefix, "_summary_txt")]] <- renderText({
    req(!is.null(rv[[raw_slot]]))
    df   <- rv[[raw_slot]]
    sel  <- input[[paste0(prefix, "_cols_selected")]]
    cols <- if (!is.null(sel) && length(sel) > 0) sel else rv[[auto_slot]]
    cols <- intersect(cols, names(df))
    req(length(cols) > 0)
    vals <- suppressWarnings(as.numeric(unlist(df[, cols])))
    ok   <- vals[!is.na(vals)]
    sprintf(paste0("Subjects (rows) : %d\n%s columns    : %d\n",
                   "Observations    : %d  (non-NA)\nMean            : %.2f\n",
                   "SD              : %.2f\nRange           : [%g, %g]"),
            nrow(df), label, length(cols), length(ok),
            mean(ok), sd(ok), min(ok), max(ok))
  })

  # --- Preview table (first 10 rows) ---
  output[[paste0(prefix, "_tbl_preview")]] <- renderDT({
    req(!is.null(rv[[raw_slot]]))
    df   <- rv[[raw_slot]]
    sel  <- input[[paste0(prefix, "_cols_selected")]]
    cols <- if (!is.null(sel) && length(sel) > 0) sel else rv[[auto_slot]]
    id_c  <- input[[paste0(prefix, "_id_col")]]
    id_c  <- if (!is.null(id_c) && id_c != "__rownum__") id_c else NULL
    grp_c <- input[[paste0(prefix, "_group_col")]]
    grp_c <- if (!is.null(grp_c) && grp_c != "__none__") grp_c else NULL
    show  <- unique(c(id_c, grp_c, intersect(cols, names(df))))
    datatable(head(df[, show, drop = FALSE], 10L),
              options  = list(pageLength = 5L, dom = "t", scrollX = TRUE),
              rownames = FALSE, class = "compact stripe")
  })

  # ---- Long-format reactive (returned to caller) ---------------
  reactive({
    req(!is.null(rv[[raw_slot]]))
    sel_cols <- input[[paste0(prefix, "_cols_selected")]]
    req(!is.null(sel_cols) && length(sel_cols) > 0)
    times_txt <- input[[paste0(prefix, "_times_txt")]]
    req(!is.null(times_txt) && nchar(trimws(times_txt)) > 0)

    val_cols <- intersect(sel_cols, names(rv[[raw_slot]]))
    if (length(val_cols) == 0L) return(NULL)

    t_vec <- tryCatch(
      as.numeric(unlist(strsplit(trimws(times_txt), "[,;\\s]+"))),
      error = function(e) numeric(0L))
    t_vec <- t_vec[!is.na(t_vec)]

    if (length(t_vec) != length(val_cols)) {
      showNotification(
        sprintf("%s: times entered (%d) must match columns selected (%d).",
                label, length(t_vec), length(val_cols)),
        type = "warning", duration = 5L)
      return(NULL)
    }
    time_map <- setNames(t_vec, val_cols)

    df <- rv[[raw_slot]]
    id_col <- input[[paste0(prefix, "_id_col")]]
    if (is.null(id_col) || id_col == "__rownum__") {
      df$.__subj__. <- seq_len(nrow(df))
      id_col <- ".__subj__."
    }
    grp_col   <- input[[paste0(prefix, "_group_col")]]
    use_group <- !is.null(grp_col) && grp_col != "__none__" && grp_col %in% names(df)

    keep <- unique(c(id_col, val_cols, if (use_group) grp_col))
    df_s <- df[, intersect(keep, names(df)), drop = FALSE]
    for (col in val_cols) df_s[[col]] <- suppressWarnings(as.numeric(df_s[[col]]))

    df_long <- tidyr::pivot_longer(df_s, cols = dplyr::all_of(val_cols),
                                   names_to = "occasion", values_to = "value") %>%
      dplyr::mutate(time_h = time_map[occasion]) %>%
      dplyr::filter(!is.na(value), !is.na(time_h))

    names(df_long)[names(df_long) == id_col] <- "subject"
    if (use_group) {
      names(df_long)[names(df_long) == grp_col] <- "group"
      df_long$group <- as.character(df_long$group)
    } else {
      df_long$group <- "All"
    }

    rec$record(paste0(prefix, "_import"),
      sprintf("Empirical %s data loaded: %d subjects, %d occasions, %d observations",
              label, length(unique(df_long$subject)), length(val_cols), nrow(df_long)),
      sprintf("# %s columns: %s\n# Time mapping: %s",
              label, paste(val_cols, collapse = ", "),
              paste(names(time_map), time_map, sep = "=", collapse = ", ")))

    df_long
  })
}

#' Generic model-vs-observed comparison outputs (House Rules 1, 6, 8)
#'
#' Builds fit-statistic value boxes, a time-series overlay, a predicted-vs-
#' observed scatter, a residuals bar chart, a comparison table, and a CSV
#' download -- all namespaced by `prefix`. Called once per outcome that has
#' a directly comparable model arm (KSS, VAS-as-KSS-equivalent, PVT lapses).
#'
#' @param emp_reactive     reactive() returning long-format empirical data
#'                         with columns subject, group, time_h, `value_col`
#' @param sim_results_fn   reactive() returning the model simulation data.frame
#' @param value_col        name of the empirical outcome column (e.g. "KSS")
#' @param model_col        name of the matching model column (e.g. "KSS")
#' @param y_range          c(lo, hi) axis range for this outcome
#' @return invisible reactive() of the merged comparison data.frame
make_comparison_outputs <- function(input, output, session, prefix, label,
                                    emp_reactive, sim_results_fn,
                                    value_col, model_col, y_range,
                                    err_type_input = "err_type") {

  # ---- Group selector -------------------------------------------------
  # Fit statistics and every plot below are computed on the SELECTED groups
  # only, so deselecting a group genuinely re-estimates RMSE / MAE / r rather
  # than merely hiding a trace.
  output[[paste0(prefix, "_ui_groups")]] <- renderUI({
    emp <- tryCatch(emp_reactive(), error = function(e) NULL)
    req(!is.null(emp) && nrow(emp) > 0)
    grps <- sort(unique(as.character(emp$group)))
    if (length(grps) <= 1) return(NULL)   # nothing to choose from
    pickerInput(paste0(prefix, "_groups_selected"),
                label    = "Groups included in plots and fit statistics:",
                choices  = grps, selected = grps, multiple = TRUE,
                options  = list(`actions-box` = TRUE,
                                `selected-text-format` = "count > 3",
                                `count-selected-text`  = "{0} of {1} groups"))
  })

  comp_data <- reactive({
    res <- sim_results_fn()
    req(!is.null(res))
    emp <- tryCatch(emp_reactive(), error = function(e) NULL)
    req(!is.null(emp) && nrow(emp) > 0)
    wake_res <- res[res$wake == 1L, ]
    req(nrow(wake_res) > 0)

    # Restrict to the user-selected groups (all groups when the picker has
    # not been rendered, i.e. a single-group dataset).
    sel_grp <- input[[paste0(prefix, "_groups_selected")]]
    if (!is.null(sel_grp) && length(sel_grp) > 0) {
      emp <- emp[as.character(emp$group) %in% sel_grp, , drop = FALSE]
    }
    req(nrow(emp) > 0)

    emp_summary <- emp %>%
      dplyr::group_by(time_h, group) %>%
      dplyr::summarise(
        obs_mean = mean(.data[[value_col]], na.rm = TRUE),
        obs_sd   = sd(.data[[value_col]],   na.rm = TRUE),
        obs_sem  = obs_sd / sqrt(dplyr::n()),
        obs_n    = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(group, time_h)

    emp_summary %>%
      dplyr::rowwise() %>%
      dplyr::mutate(pred = approx(wake_res$time_h, wake_res[[model_col]],
                                  xout = time_h, method = "linear", rule = 2)$y) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(residual = pred - obs_mean, abs_error = abs(residual))
  })

  fit_stat <- function(fn) {
    tryCatch({
      cd <- comp_data(); req(!is.null(cd) && nrow(cd) >= 2)
      d <- cd[!is.na(cd$pred) & !is.na(cd$obs_mean), ]
      if (nrow(d) < 2) return(NA_real_)
      fn(d$pred, d$obs_mean)
    }, error = function(e) NA_real_)
  }

  span <- diff(y_range)

  # Subtitle suffix naming the groups the statistic was computed on
  grp_suffix <- reactive({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    if (is.null(cd) || nrow(cd) == 0) return("")
    g <- sort(unique(as.character(cd$group)))
    if (length(g) == 1 && identical(g, "All")) return("")
    if (length(g) <= 2) paste0(" | ", paste(g, collapse = " + "))
    else paste0(" | ", length(g), " groups")
  })

  output[[paste0(prefix, "_vb_rmse")]] <- renderValueBox({
    val <- fit_stat(function(p, o) sqrt(mean((p - o)^2)))
    col <- if (is.na(val)) "blue" else if (val < 0.10 * span) "green"
            else if (val < 0.20 * span) "yellow" else "red"
    valueBox(if (is.na(val)) "--" else sprintf("%.2f", val),
             paste0("RMSE (", label, ")", grp_suffix()),
             icon = icon("bullseye"), color = col)
  })
  output[[paste0(prefix, "_vb_mae")]] <- renderValueBox({
    val <- fit_stat(function(p, o) mean(abs(p - o)))
    col <- if (is.na(val)) "blue" else if (val < 0.10 * span) "green"
            else if (val < 0.20 * span) "yellow" else "red"
    valueBox(if (is.na(val)) "--" else sprintf("%.2f", val),
             paste0("MAE (", label, ")", grp_suffix()),
             icon = icon("ruler-combined"), color = col)
  })
  output[[paste0(prefix, "_vb_r")]] <- renderValueBox({
    val <- fit_stat(function(p, o) if (sd(p) < 1e-9 || sd(o) < 1e-9) NA_real_ else cor(p, o))
    col <- if (is.na(val)) "blue" else if (val > 0.70) "green"
            else if (val > 0.40) "yellow" else "red"
    valueBox(if (is.na(val)) "--" else sprintf("%.3f", val),
             paste0("Pearson r (pred ~ obs)", grp_suffix()),
             icon = icon("chart-line"), color = col)
  })

  output[[paste0(prefix, "_plt_ts")]] <- renderPlotly({
    res <- sim_results_fn(); req(!is.null(res))
    cd  <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) > 0)
    wake <- res[res$wake == 1L, ]

    groups   <- sort(unique(cd$group))
    grp_cols <- group_colours(groups)
    err_col  <- if (isTRUE(input[[err_type_input]] == "sem")) "obs_sem" else "obs_sd"

    p <- plot_ly(wake, x = ~time_h, y = wake[[model_col]],
                type = "scatter", mode = "lines",
                name = paste0(label, " model"),
                line = list(color = "#2980b9", width = 2.5))

    for (grp in groups) {
      d_g <- cd[cd$group == grp, ]
      col <- unname(grp_cols[[grp]])
      err <- d_g[[err_col]]
      lbl <- if (length(groups) <= 1) "Observed mean" else paste0("Observed: ", grp)
      p <- add_trace(p, data = d_g, x = ~time_h, y = ~obs_mean,
                     type = "scatter", mode = "markers+lines",
                     error_y = list(type = "data", array = err, color = col,
                                    thickness = 2, width = 8),
                     name = lbl, marker = list(color = col, size = 10),
                     line = list(color = col, width = 1.5, dash = "dot"),
                     text = ~paste0("t = ", round(time_h, 1), " h<br>Obs = ", round(obs_mean, 2),
                                    " +/- ", round(err, 2), " (n = ", obs_n, ")<br>",
                                    "Model = ", round(pred, 2), "<br>Residual = ", round(residual, 2)),
                     hoverinfo = "text")
    }
    p %>% layout(xaxis = list(title = "Time (h from study start)", zeroline = FALSE),
                yaxis = list(title = label, range = y_range, gridcolor = "#eee"),
                shapes = sleep_shapes(res),
                legend = list(orientation = "h", y = -0.28),
                margin = list(t = 5, r = 10),
                plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe")
  })

  output[[paste0(prefix, "_plt_scatter")]] <- renderPlotly({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) >= 2)
    groups   <- sort(unique(cd$group))
    grp_cols <- group_colours(groups)
    err_col  <- if (isTRUE(input[[err_type_input]] == "sem")) "obs_sem" else "obs_sd"

    vals  <- c(cd$pred, cd$obs_mean)
    ax_lo <- max(floor(min(vals, na.rm = TRUE) - 0.05 * span), y_range[1])
    ax_hi <- min(ceiling(max(vals, na.rm = TRUE) + 0.05 * span), y_range[2])
    ax_range <- c(ax_lo, ax_hi)

    p <- plot_ly() %>%
      add_lines(x = ax_range, y = ax_range, name = "y = x (perfect fit)",
               line = list(color = "#999", dash = "dash", width = 1.5))

    for (grp in groups) {
      d_g <- cd[cd$group == grp & !is.na(cd$pred), ]
      col <- unname(grp_cols[[grp]])
      err <- d_g[[err_col]]
      lbl <- if (length(groups) <= 1) "Occasions" else grp
      p <- add_trace(p, data = d_g, x = ~pred, y = ~obs_mean,
                     type = "scatter", mode = "markers",
                     error_y = list(type = "data", array = err, color = col, thickness = 2, width = 5),
                     name = lbl, marker = list(color = col, size = 10,
                                               line = list(color = "white", width = 1)),
                     text = ~paste0("t = ", round(time_h, 1), " h<br>Model = ", round(pred, 2),
                                    "<br>Obs = ", round(obs_mean, 2), " +/- ", round(err, 2),
                                    "<br>Residual = ", round(residual, 2)),
                     hoverinfo = "text")
    }
    p %>% layout(xaxis  = list(title = paste0("Model predicted ", label), range = ax_range, gridcolor = "#eee"),
                yaxis  = list(title = paste0("Observed mean ", label),   range = ax_range, gridcolor = "#eee"),
                legend = list(orientation = "h", y = -0.28), margin = list(t = 5),
                plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe")
  })

  output[[paste0(prefix, "_plt_residuals")]] <- renderPlotly({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) > 0)
    groups <- sort(unique(cd$group))
    p <- plot_ly()
    for (grp in groups) {
      d_g <- cd[cd$group == grp & !is.na(cd$residual), ]
      d_g <- d_g[order(d_g$time_h), ]
      col_vec <- ifelse(d_g$residual >= 0, "#2980b9", "#c0392b")
      lbl <- if (length(groups) <= 1) "Residual" else paste0("Residual: ", grp)
      p <- add_bars(p, data = d_g, x = ~time_h, y = ~residual, name = lbl,
                    marker = list(color = col_vec, line = list(color = "white", width = 0.5)),
                    text = ~paste0("t = ", round(time_h, 1), " h<br>Residual = ", sprintf("%+.2f", residual)),
                    hoverinfo = "text")
    }
    p %>% layout(xaxis = list(title = "Time (h from study start)", zeroline = FALSE),
                yaxis = list(title = "Residual (pred - obs)", zeroline = TRUE,
                             zerolinecolor = "#333", zerolinewidth = 1.5, gridcolor = "#eee"),
                barmode = "group", legend = list(orientation = "h", y = -0.35), margin = list(t = 5),
                plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe")
  })

  output[[paste0(prefix, "_tbl")]] <- renderDT({
    cd <- tryCatch(comp_data(), error = function(e) NULL)
    req(!is.null(cd) && nrow(cd) > 0)
    err_col <- if (isTRUE(input[[err_type_input]] == "sem")) "obs_sem" else "obs_sd"
    err_lbl <- if (isTRUE(input[[err_type_input]] == "sem")) "Obs SEM" else "Obs SD"
    disp <- cd %>% dplyr::arrange(time_h, group) %>%
      dplyr::transmute(
        `Time (h)` = round(time_h, 2),
        Group      = group,
        `Obs mean` = round(obs_mean, 2),
        !!err_lbl  := round(.data[[err_col]], 2),
        N          = obs_n,
        Model      = round(pred, 2),
        Residual   = round(residual, 2),
        `Abs error`= round(abs_error, 2)
      )
    datatable(disp, options = list(pageLength = 20, dom = "tp", scrollX = TRUE),
              rownames = FALSE, class = "compact stripe") %>%
      formatStyle("Residual", backgroundColor = styleInterval(
        c(-0.15 * span, -0.05 * span, 0.05 * span, 0.15 * span),
        c("#f5b7b1", "#fdebd0", "transparent", "#fdebd0", "#f5b7b1")))
  })

  output[[paste0(prefix, "_dl")]] <- downloadHandler(
    filename = function() paste0("farima_comparison_", prefix, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      cd <- tryCatch(comp_data(), error = function(e) NULL)
      req(!is.null(cd))
      write.csv(cd, file, row.names = FALSE)
    }
  )

  invisible(comp_data)
}


# ============================================================
# 3. UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "FaRiMaTo 4.0"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("Schedule",   tabName = "schedule", icon = icon("bed")),
      menuItem("Parameters", tabName = "params",   icon = icon("sliders-h")),
      menuItem("Results",    tabName = "results",  icon = icon("chart-line")),
      menuItem("Comparison", tabName = "compare",  icon = icon("not-equal")),
      menuItem("Provenance", tabName = "prov",     icon = icon("triangle-exclamation")),
      menuItem("Export",     tabName = "export",   icon = icon("download"))
    ),
    hr(),
    div(style = "padding: 10px;",
      p(style = "color:#aaa; font-size:11px;",
        "FaRiMaTo v4.0",
        br(), "McCauley / Van Dongen model",
        br(), "2021 / 2024jtb / 2024fe versions",
        br(), tags$b(style="color:#f1c40f;", "Risk score is NOT from the paper."))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-header .box-title { font-weight: 600; }
      pre { font-size: 11px; }
    "))),

    tabItems(

      # ========== SCHEDULE TAB ==========
      tabItem(tabName = "schedule",
        fluidRow(
          box(title = "Study reference time", width = 4, status = "primary", solidHeader = TRUE,
            numericInput("t0_clock", "Clock hour at t = 0 (from midnight):",
                        value = 8, min = 0, max = 23.5, step = 0.5),
            helpText("Sets the circadian phase reference for all predictions.")
          ),
          box(title = "Input mode", width = 4, status = "primary", solidHeader = TRUE,
            radioButtons("input_mode", NULL,
              choices = c("Manual schedule entry" = "manual",
                          "Upload file (onset / offset)"  = "file",
                          "Upload actigraphy (datetime + state)" = "actig"),
              selected = "manual"),
            conditionalPanel("input.input_mode == 'file'",
              fileInput("schedule_file", "File (CSV / xlsx):",
                       accept = c(".csv", ".txt", ".xlsx", ".xls")),
              helpText("Columns required: onset, offset (ISO datetime or decimal hours from t=0)")
            ),
            conditionalPanel("input.input_mode == 'actig'",
              fileInput("actig_file", "Actigraphy file (CSV):",
                       accept = c(".csv", ".txt")),
              helpText("Columns required: datetime (ISO), sleep_wake (0 = sleep, 1 = wake)")
            )
          ),
          box(title = "Run", width = 4, status = "success", solidHeader = TRUE,
            actionButton("btn_example", "Load 7-day example",
                        class = "btn-default btn-block", icon = icon("magic")),
            br(),
            checkboxInput("auto_recalc",
                         "Live recalculation (auto-update on any change)",
                         value = TRUE),
            conditionalPanel("!input.auto_recalc",
              actionButton("btn_run", "Run Simulation",
                          class = "btn-success btn-lg btn-block")
            ),
            br(),
            verbatimTextOutput("sim_status_txt", placeholder = TRUE)
          )
        ),

        # Manual entry panel
        conditionalPanel("input.input_mode == 'manual'",
          fluidRow(
            box(title = "Sleep schedule (hours from t = 0)", width = 12,
                status = "info", solidHeader = TRUE,
              fluidRow(
                column(3, numericInput("new_onset", "Sleep onset (h from t=0):",
                                      value = 15, min = 0, max = 500, step = 0.5)),
                column(3, numericInput("new_offset", "Wake time (h from t=0):",
                                      value = 23, min = 0, max = 500, step = 0.5)),
                column(3, br(),
                  actionButton("btn_add", "Add period", class = "btn-info"),
                  actionButton("btn_del", "Delete selected", class = "btn-warning")
                ),
                column(3, br(),
                  helpText("Onset 15 h -> 23:00 if study starts at 08:00.",
                          br(), "Select row(s) to delete.")
                )
              ),
              DTOutput("tbl_schedule_manual")
            )
          )
        ),

        # File preview
        conditionalPanel("input.input_mode != 'manual'",
          fluidRow(
            box(title = "Parsed schedule preview (hours from t=0)", width = 12,
                status = "info", solidHeader = TRUE,
              DTOutput("tbl_file_preview")
            )
          )
        ),

        # ---- Empirical KSS data import (optional) -------------------
        fluidRow(
          box(title = "Empirical KSS overlay (optional)",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsible = TRUE, collapsed = FALSE,
            fluidRow(
              column(4,
                fileInput("kss_file", "Upload KSS file (CSV / xlsx):",
                         accept = c(".csv", ".txt", ".xlsx", ".xls")),
                helpText(
                  "Wide format: one row per subject, one column per KSS measurement.",
                  br(), "Column names matching kss1, kss2, ... are auto-detected.",
                  br(), "UTF-8 with BOM, comma- and semicolon-delimited files are supported."
                )
              ),
              column(4,
                conditionalPanel("output.kss_loaded",
                  uiOutput("kss_ui_cols"),
                  uiOutput("kss_ui_id"),
                  uiOutput("kss_ui_group")
                )
              ),
              column(4,
                conditionalPanel("output.kss_loaded",
                  uiOutput("kss_ui_times"),
                  br(),
                  checkboxInput("kss_show_indiv", "Show individual subject traces", value = FALSE),
                  checkboxInput("kss_show_mean",  "Show mean +/- SD band", value = TRUE)
                )
              )
            ),
            conditionalPanel("output.kss_loaded",
              hr(),
              fluidRow(
                column(5, verbatimTextOutput("kss_summary_txt")),
                column(7, DTOutput("kss_tbl_preview"))
              )
            )
          )
        ),

        # ---- Empirical VAS data import (optional) --------------------
        fluidRow(
          box(title = "Empirical VAS sleepiness overlay (optional, 0-100 scale)",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
            fluidRow(
              column(4,
                fileInput("vas_file", "Upload VAS file (CSV / xlsx):",
                         accept = c(".csv", ".txt", ".xlsx", ".xls")),
                helpText(
                  "Wide format, one column per VAS measurement.",
                  br(), "Columns matching vas1, vas2, ... are auto-detected;",
                  " otherwise select manually below (e.g. custom-named VAS columns)."
                )
              ),
              column(4,
                conditionalPanel("output.vas_loaded",
                  uiOutput("vas_ui_cols"),
                  uiOutput("vas_ui_id"),
                  uiOutput("vas_ui_group")
                )
              ),
              column(4,
                conditionalPanel("output.vas_loaded",
                  uiOutput("vas_ui_times"),
                  br(),
                  checkboxInput("vas_show_indiv", "Show individual subject traces", value = FALSE),
                  checkboxInput("vas_show_mean",  "Show mean +/- SD band", value = TRUE)
                )
              )
            ),
            conditionalPanel("output.vas_loaded",
              hr(),
              fluidRow(
                column(4, numericInput("vas_min", "VAS raw minimum:", value = 0,   step = 1)),
                column(4, numericInput("vas_max", "VAS raw maximum:", value = 100, step = 1)),
                column(4, radioButtons("vas_direction", "Scale direction:",
                                       choices = c("0 = alert, high = sleepy" = "direct",
                                                   "0 = sleepy, high = alert" = "reverse"),
                                       selected = "direct"))
              ),
              helpText(
                "Rescaled linearly onto the KSS range (1-9) using the endpoints above,",
                " so it can be plotted and compared against the model's KSS arm.",
                br(), "This is a transparent linear endpoint transform, not a validated",
                " regression equation, and no published conversion is implemented here.",
                " Treat the KSS-equivalent as a rescaling of your own VAS, nothing more."
              ),
              fluidRow(
                column(5, verbatimTextOutput("vas_summary_txt")),
                column(7, DTOutput("vas_tbl_preview"))
              )
            )
          )
        ),

        # ---- Empirical PVT data import (optional) --------------------
        fluidRow(
          box(title = "Empirical PVT overlay (optional)",
              width = 12, status = "warning", solidHeader = TRUE,
              collapsible = TRUE, collapsed = TRUE,
            fluidRow(
              column(4,
                fileInput("pvt_file", "Upload PVT lapse file (CSV / xlsx):",
                         accept = c(".csv", ".txt", ".xlsx", ".xls")),
                helpText(
                  "Expected metric: ", tags$b("lapses (count per 10-min bout)"),
                  " -- the same unit as the model's PVT arm.",
                  br(), "Columns matching pvt1, pvt2, ... are auto-detected.",
                  br(), tags$em("Reaction time is not accepted here: the model predicts lapse",
                                " counts, not latencies, so RT has no comparable model arm.")
                )
              ),
              column(4,
                conditionalPanel("output.pvt_loaded",
                  uiOutput("pvt_ui_cols"),
                  uiOutput("pvt_ui_id"),
                  uiOutput("pvt_ui_group")
                )
              ),
              column(4,
                conditionalPanel("output.pvt_loaded",
                  uiOutput("pvt_ui_times"),
                  br(),
                  checkboxInput("pvt_show_indiv", "Show individual subject traces", value = FALSE),
                  checkboxInput("pvt_show_mean",  "Show mean +/- SD band", value = TRUE)
                )
              )
            ),
            conditionalPanel("output.pvt_loaded",
              hr(),
              fluidRow(
                column(5, verbatimTextOutput("pvt_summary_txt")),
                column(7, DTOutput("pvt_tbl_preview"))
              )
            )
          )
        )

      ),  # end schedule tabItem


      # ========== PARAMETERS TAB ==========
      tabItem(tabName = "params",
        fluidRow(
          box(width = 12, status = "success", solidHeader = TRUE,
              title = "Published parameter values -- McCauley et al. (2024), Table 2",
            p("The sliders below are initialised to the ",
              tags$b("published calibrated values"),
              " and the equations are those printed in Table 1 of that paper. ",
              "Moving a slider departs from the published calibration; the ",
              tags$b("Restore published values"), " button puts it back."),
            p(tags$b("Note on metric-specific parameters: "),
              "p_f, alpha_w, alpha_s, W_c, mu_w and xi take different published ",
              "values for PVT and for KSS. The model is run twice, once per metric. ",
              "There is no link function converting between them."),
            fluidRow(
              column(7,
                selectInput("model_version", "Published model version:",
                            choices = model_versions, selected = "2024fe",
                            width = "100%")),
              column(5, br(),
                actionButton("btn_restore_pub", "Restore published values",
                             class = "btn-success", icon = icon("rotate-left")))
            ),
            uiOutput("version_note")
          )
        ),

        fluidRow(
          box(title = "Simulation settings", width = 6, status = "warning", solidHeader = TRUE,
            selectInput("dt_choice", "Integration step dt (Euler):",
                       choices = c("3 min (0.05 h, default)" = 0.05,
                                   "1 min (0.0167 h)"        = 0.0167,
                                   "30 s (0.0083 h)"         = 0.0083),
                       selected = 0.05),
            helpText("The published model was not integrated with a stated step size; ",
                     "3 min is used here as a default. Check that results are stable ",
                     "when you reduce it.")
          ),
          box(title = "Initial conditions (Section 2.4 of the paper)",
              width = 6, status = "warning", solidHeader = TRUE,
            checkboxInput("equilibrate",
                          "Equilibrate on a repeating baseline day (as in the paper)",
                          value = TRUE),
            fluidRow(
              column(6, numericInput("base_onset",  "Baseline sleep onset (clock h):",
                                     value = 23, min = 0, max = 23.75, step = 0.25)),
              column(6, numericInput("base_offset", "Baseline wake time (clock h):",
                                     value = 7,  min = 0, max = 23.75, step = 0.25))
            ),
            helpText("The paper iterates the baseline schedule until a steady state at ",
                     "the scheduled awakening time (tolerance 1e-9). Set the baseline to ",
                     "your study's own pre-study schedule."),
            helpText(tags$b("Bifurcation limit: "),
                     "the baseline must leave less than W_c of daily wakefulness ",
                     "(20.2 h for PVT, 22.02 h for KSS at published values), i.e. at ",
                     "least ~3.8 h of baseline sleep. Above W_c the model escalates ",
                     "without a steady state and equilibration is refused.")
          )
        ),

        fluidRow(
          box(title = "Shared parameters (same value for both metrics)",
              width = 6, status = "primary", solidHeader = TRUE,
            sliderInput("beta",   "beta -- impact of u(t) onto p(t):",        0, 1,    0.26, 0.01),
            sliderInput("zeta",   "zeta -- impact of p(t) onto u(t) [2024]:", 0, 3,    1.31, 0.01),
            sliderInput("eta_w",  "eta_w -- build-up rate for u(t) (h-1):",   0, 0.05, 0.0126, 0.0001),
            sliderInput("nu",     "nu -- decay rate for h(t):",               0, 4,    1.37, 0.01),
            sliderInput("gamma",  "gamma -- asymptote scaling for h(t):",     0, 2,    0.71, 0.01),
            sliderInput("lambda", "lambda -- circadian amplitude rate (h-1):",0, 2,    0.49, 0.01)
          ),
          box(title = "Circadian process (shared)", width = 6, status = "info", solidHeader = TRUE,
            sliderInput("phi",   "phi -- circadian phase vs midnight (h):", 0, 24, 21.2, 0.1),
            sliderInput("tau",   "tau -- circadian period (h):",         23, 25, 24.0, 0.1),
            sliderInput("mu_s",  "mu_s -- circadian offset during sleep:", -3, 1, -1.5, 0.05),
            helpText("The published model has a single sine circadian term (Equation T1-6) ",
                     "with a fixed phase. It has no chronotype parameter; the earlier ",
                     "chronotype presets in this app were not from any source and have ",
                     "been removed.")
          )
        ),

        fluidRow(
          box(title = "PVT-specific parameters (lapses / 10 min)",
              width = 6, status = "danger", solidHeader = TRUE,
            sliderInput("pf_pvt",      "p_f -- floor of metric scale:",        0, 2,   0,     0.1),
            sliderInput("alpha_w_pvt", "alpha_w -- build-up rate (h-1):",      0, 0.5, 0.028, 0.001),
            sliderInput("alpha_s_pvt", "alpha_s -- dissipation rate (h-1):",   0, 0.5, 0.26,  0.001),
            sliderInput("Wc_pvt",      "W_c -- bifurcation threshold (h):",   12, 24,  20.2,  0.01),
            sliderInput("mu_w_pvt",    "mu_w -- circadian offset, wake:",      0, 2,   0.466, 0.001),
            sliderInput("xi_pvt",      "xi -- metric scaling factor:",         0, 3,   1.09,  0.01)
          ),
          box(title = "KSS-specific parameters (scale 1-9)",
              width = 6, status = "danger", solidHeader = TRUE,
            sliderInput("pf_kss",      "p_f -- floor of metric scale:",        0, 2,   1,     0.1),
            sliderInput("alpha_w_kss", "alpha_w -- build-up rate (h-1):",      0, 0.5, 0.22,  0.001),
            sliderInput("alpha_s_kss", "alpha_s -- dissipation rate (h-1):",   0, 0.5, 0.037, 0.001),
            sliderInput("Wc_kss",      "W_c -- bifurcation threshold (h):",   12, 24,  22.02, 0.01),
            sliderInput("mu_w_kss",    "mu_w -- circadian offset, wake:",      0, 2,   0.82,  0.001),
            sliderInput("xi_kss",      "xi -- metric scaling factor:",         0, 3,   0.51,  0.01)
          )
        ),

        fluidRow(
          box(title = "Composite risk score -- NOT from the paper",
              width = 12, status = "danger", solidHeader = TRUE,
            p(tags$b(style = "color:#c0392b;",
                     "This score is the app author's construction. "),
              "The published model predicts PVT lapses and KSS sleepiness; it does ",
              "not define a combined risk index, and the weights and cut-points below ",
              "have no empirical basis. Report the two predictions, not this score."),
            fluidRow(
              column(4, sliderInput("w_kss", "Weight on KSS:", 0, 1, 0.5, 0.05)),
              column(4, sliderInput("w_pvt", "Weight on PVT:", 0, 1, 0.5, 0.05)),
              column(4, numericInput("pvt_ceiling",
                                     "PVT lapses treated as maximum:", value = 20,
                                     min = 1, max = 100, step = 1))
            ),
            textOutput("weight_check")
          )
        )
      ),  # end params tabItem


      # ========== RESULTS TAB ==========
      tabItem(tabName = "results",
        conditionalPanel("!output.sim_done",
          fluidRow(box(width = 12, status = "warning",
            h4("No valid schedule yet -- add at least one sleep period on the Schedule tab.")
          ))
        ),
        conditionalPanel("output.sim_done",
          fluidRow(
            valueBoxOutput("vb_kss",      width = 4),
            valueBoxOutput("vb_pvt",      width = 4),
            valueBoxOutput("vb_highrisk", width = 4)
          ),
          fluidRow(
            box(title = "Subjective Sleepiness (KSS 1-9) & PVT Lapses / 10 min -- with empirical overlays",
                width = 12, status = "primary", solidHeader = TRUE,
              plotlyOutput("plt_kss_pvt", height = "340px"),
              helpText("Circles = observed KSS, triangles = VAS->KSS-equivalent, squares = observed PVT lapses (when loaded).")
            )
          ),
          fluidRow(
            box(title = "Composite Risk Score  (green < 0.33 | amber 0.33-0.65 | red >= 0.66)",
                width = 12, status = "danger", solidHeader = TRUE,
              fluidRow(
                column(5,
                  checkboxInput("risk_overlay_sleep",
                                "Overlay observed sleepiness (KSS / VAS) on right axis",
                                value = FALSE)),
                column(7,
                  conditionalPanel("input.risk_overlay_sleep",
                    radioButtons("risk_overlay_scale", NULL, inline = TRUE,
                      choices = c("Right axis in KSS units (1-9)" = "kss",
                                  "Mapped onto the risk scale [0-1]" = "risk"),
                      selected = "kss")))
              ),
              plotlyOutput("plt_risk", height = "280px"),
              conditionalPanel("input.risk_overlay_sleep && input.risk_overlay_scale == 'risk'",
                helpText("Mapped values use only the KSS term of the composite, (KSS-1)/8, so they are NOT an observed risk score -- the PVT term cannot be derived from a sleepiness rating."))
            )
          ),
          fluidRow(
            box(title = "State variables: u(t), h(t), kappa(t), c(t)",
                width = 12, status = "info", solidHeader = TRUE,
              plotlyOutput("plt_states", height = "300px")
            )
          ),
          fluidRow(
            box(title = "Hourly prediction table (wake periods only)",
                width = 12, status = "success", solidHeader = TRUE,
              DTOutput("tbl_preds"),
              br(),
              downloadButton("dl_csv", "Download full predictions CSV",
                            class = "btn-sm btn-success")
            )
          )
        )
      ),  # end results tabItem


      # ========== COMPARISON TAB ==========
      tabItem(tabName = "compare",

        conditionalPanel("!output.sim_done || (!output.kss_loaded && !output.vas_loaded && !output.pvt_loaded)",
          fluidRow(
            box(width = 12, status = "warning", solidHeader = FALSE,
              h4("Load at least one empirical dataset to activate this tab:"),
              tags$ol(
                tags$li("A valid schedule/simulation (Schedule tab)"),
                tags$li("KSS, VAS, and/or PVT data (Schedule tab, collapsible boxes below the schedule)")
              )
            )
          )
        ),

        conditionalPanel("output.sim_done && (output.kss_loaded || output.vas_loaded || output.pvt_loaded)",
          fluidRow(
            box(width = 12, status = "primary", solidHeader = FALSE,
              radioButtons("err_type", "Error bars show:",
                           choices  = c("+/- SD" = "sd", "+/- SEM" = "sem"),
                           selected = "sd", inline = TRUE),
              helpText("SD = individual variability; SEM = precision of the mean. Applies to all sub-tabs below.")
            )
          ),

          tabsetPanel(id = "compare_subtabs",

            # ---- KSS sub-tab ----
            tabPanel("KSS",
              conditionalPanel("!output.kss_loaded",
                fluidRow(box(width = 12, status = "warning", h5("No empirical KSS data loaded.")))
              ),
              conditionalPanel("output.kss_loaded",
                fluidRow(
                  box(width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("kss_ui_groups"))
                ),
                fluidRow(
                  valueBoxOutput("kss_vb_rmse", width = 4),
                  valueBoxOutput("kss_vb_mae",  width = 4),
                  valueBoxOutput("kss_vb_r",    width = 4)
                ),
                fluidRow(
                  box(title = "Model vs Observed KSS -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("kss_plt_ts", height = "360px")),
                  box(title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("kss_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("kss_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("kss_tbl"), br(),
                    downloadButton("kss_dl", "Download CSV", class = "btn-sm btn-success"))
                )
              )
            ),

            # ---- VAS sub-tab ----
            tabPanel("VAS (-> KSS-equiv)",
              conditionalPanel("!output.vas_loaded",
                fluidRow(box(width = 12, status = "warning", h5("No empirical VAS data loaded.")))
              ),
              conditionalPanel("output.vas_loaded",
                fluidRow(
                  box(width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("vas_ui_groups"))
                ),
                fluidRow(
                  valueBoxOutput("vas_vb_rmse", width = 4),
                  valueBoxOutput("vas_vb_mae",  width = 4),
                  valueBoxOutput("vas_vb_r",    width = 4)
                ),
                fluidRow(
                  box(title = "Model KSS vs VAS-derived KSS-equivalent -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("vas_plt_ts", height = "360px")),
                  box(title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("vas_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("vas_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("vas_tbl"), br(),
                    downloadButton("vas_dl", "Download CSV", class = "btn-sm btn-success"))
                ),
                helpText("KSS-equivalent obtained via linear endpoint rescaling of the raw VAS (see Schedule tab note); treat as approximate.")
              )
            ),

            # ---- PVT sub-tab ----
            tabPanel("PVT",
              conditionalPanel("!output.pvt_loaded",
                fluidRow(box(width = 12, status = "warning", h5("No empirical PVT data loaded.")))
              ),
              conditionalPanel("output.pvt_loaded",
                fluidRow(
                  box(width = 12, status = "primary", solidHeader = FALSE,
                      uiOutput("pvtlap_ui_groups"))
                ),
                fluidRow(
                  valueBoxOutput("pvtlap_vb_rmse", width = 4),
                  valueBoxOutput("pvtlap_vb_mae",  width = 4),
                  valueBoxOutput("pvtlap_vb_r",    width = 4)
                ),
                fluidRow(
                  box(title = "Model vs Observed PVT lapses -- time series", width = 8,
                      status = "primary", solidHeader = TRUE,
                    plotlyOutput("pvtlap_plt_ts", height = "360px")),
                  box(title = "Predicted vs Observed", width = 4,
                      status = "warning", solidHeader = TRUE,
                    plotlyOutput("pvtlap_plt_scatter", height = "360px"))
                ),
                fluidRow(
                  box(title = "Residuals (Model - Observed)", width = 12,
                      status = "danger", solidHeader = TRUE,
                    plotlyOutput("pvtlap_plt_residuals", height = "200px"))
                ),
                fluidRow(
                  box(title = "Comparison table", width = 12,
                      status = "success", solidHeader = TRUE,
                    DTOutput("pvtlap_tbl"), br(),
                    downloadButton("pvtlap_dl", "Download CSV", class = "btn-sm btn-success"))
                )
              )
            )
          )
        )
      ),  # end compare tabItem


      # ========== PROVENANCE TAB ==========
      tabItem(tabName = "prov",
        fluidRow(
          box(width = 12, status = "danger", solidHeader = TRUE,
              title = "Parameter provenance -- read before reporting any result",
            h4("None of the parameter values used in this app come from a publication."),
            p("The functional forms (saturating-exponential homeostat, ",
              "two-harmonic circadian term, allostatic modulation, exponential ",
              "decay of sleep inertia) are conventional in the fatigue-modelling ",
              "literature. The ", tags$b("numbers"), " are not. They were chosen ",
              "to produce output in plausible ranges and have not been estimated ",
              "from data."),
            p(tags$b("Correction of record: "), "earlier builds of this tool labelled ",
              "these values with specific citations. Those attributions were wrong ",
              "and have been removed. In particular, the sleep-inertia constant 0.71 ",
              "was described as a published decay time in hours; in the one source ",
              "verified here it is a dimensionless scaling factor with a different role."),
            p(tags$b("What to do: "), "use the Comparison tab to estimate the output-arm ",
              "coefficients on your own KSS / VAS / PVT data, and report them as your ",
              "estimates in your sample. Do not present any default below as a ",
              "literature value.")
          )
        ),
        fluidRow(
          box(width = 12, status = "warning", solidHeader = TRUE,
              title = "Status of each parameter in this app",
            DTOutput("tbl_prov_status")
          )
        ),
        fluidRow(
          box(width = 12, status = "info", solidHeader = TRUE,
              title = "Verified published reference values (NOT used by this app)",
            p("Transcribed from Table 2 of: McCauley ME, McCauley P, Kalachev LV, ",
              "Banks S, Dinges DF, Van Dongen HPA (2024). The dynamics of ",
              "neurobehavioral impairment and recovery sleep: improved ",
              "biomathematical modeling for fatigue risk management in operational ",
              "settings. ", tags$em("Frontiers in Environmental Health"), " 3:1362755. ",
              tags$a(href = "https://doi.org/10.3389/fenvh.2024.1362755",
                     target = "_blank", "doi:10.3389/fenvh.2024.1362755")),
            p(tags$b("These are not wired into the simulation."), " They belong to a ",
              "different ODE system, in which the PVT/KSS distinction is handled ",
              "inside the differential equations via metric-specific rate constants, ",
              "a metric floor, a bifurcation threshold and a scaling factor -- not by ",
              "a link function applied to a shared impairment term, which is what this ",
              "app currently does. That paper explicitly argues against the ",
              "link-function approach. The table is kept here as the target for a ",
              "future re-implementation."),
            DTOutput("tbl_prov_published")
          )
        )
      ),  # end provenance tabItem


      # ========== EXPORT TAB ==========
      tabItem(tabName = "export",
        fluidRow(
          box(title = "Reproducible R Script (House Rule 6)", width = 12,
              status = "primary", solidHeader = TRUE,
            downloadButton("dl_rscript", "Download R script",
                          class = "btn-primary"),
            br(), br(),
            verbatimTextOutput("rscript_preview")
          )
        ),
        fluidRow(
          box(title = "Pipeline Log (House Rule 8)", width = 12,
              status = "info", solidHeader = TRUE,
            verbatimTextOutput("pipeline_log")
          )
        )
      )  # end export tabItem

    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage


# ============================================================
# 4. SERVER
# ============================================================
server <- function(input, output, session) {

  rec <- new_recorder()   # one step recorder per session

  # ---- Reactive storage -------------------------------------------
  rv <- reactiveValues(
    # Empty by default -- the user defines their own protocol.
    schedule_manual = data.frame(
      onset_h  = numeric(0),
      offset_h = numeric(0)
    ),
    kss_raw       = NULL, kss_cols_auto = character(0),
    vas_raw       = NULL, vas_cols_auto = character(0),
    pvt_raw       = NULL, pvt_cols_auto = character(0)
  )

  # ---- Manual schedule UI -----------------------------------------
  output$tbl_schedule_manual <- renderDT({
    datatable(rv$schedule_manual,
              selection = "multiple",
              editable  = list(target = "cell", disable = list(columns = 0L)),
              options   = list(pageLength = 15, dom = "t", ordering = FALSE),
              rownames  = FALSE,
              colnames  = c("Onset (h from t=0)", "Offset (h from t=0)"),
              class     = "compact stripe")
  })

  observeEvent(input$tbl_schedule_manual_cell_edit, {
    info <- input$tbl_schedule_manual_cell_edit
    rv$schedule_manual[info$row, info$col + 1] <- as.numeric(info$value)
  })

  observeEvent(input$btn_add, {
    new_row <- data.frame(onset_h  = as.numeric(input$new_onset),
                          offset_h = as.numeric(input$new_offset))
    rv$schedule_manual <- rbind(rv$schedule_manual, new_row)
  })

  observeEvent(input$btn_del, {
    sel <- input$tbl_schedule_manual_rows_selected
    if (!is.null(sel) && length(sel) > 0)
      rv$schedule_manual <- rv$schedule_manual[-sel, , drop = FALSE]
  })

  # ---- Load 7-day example -----------------------------------------
  observeEvent(input$btn_example, {
    rv$schedule_manual <- data.frame(
      onset_h  = 15 + (0:6) * 24,
      offset_h = 23 + (0:6) * 24
    )
    updateRadioButtons(session, "input_mode", selected = "manual")
  })

  # ---- Parse uploaded schedule file ---------------------------------
  parsed_file_schedule <- reactive({
    req(input$input_mode %in% c("file", "actig"))

    if (input$input_mode == "file") {
      req(input$schedule_file)
      f   <- input$schedule_file
      ext <- tolower(tools::file_ext(f$name))
      df  <- tryCatch({
        if (ext %in% c("xlsx", "xls")) readxl::read_excel(f$datapath)
        else readr::read_delim(f$datapath, show_col_types = FALSE)
      }, error = function(e) { showNotification(paste("Parse error:", e$message), type = "error"); NULL })
      req(!is.null(df))

      if (is.numeric(df[[1]])) {
        colnames(df)[1:2] <- c("onset_h", "offset_h")
        return(df[, c("onset_h", "offset_h")])
      } else {
        parse_orders <- c("YmdHM", "YmdHMS", "Ymd HM", "Ymd HMS", "dmy HM", "dmy HMS")
        onset_dt  <- lubridate::parse_date_time(as.character(df[[1]]), orders = parse_orders)
        offset_dt <- lubridate::parse_date_time(as.character(df[[2]]), orders = parse_orders)
        t0 <- min(onset_dt, na.rm = TRUE)
        return(data.frame(
          onset_h  = as.numeric(difftime(onset_dt,  t0, units = "hours")),
          offset_h = as.numeric(difftime(offset_dt, t0, units = "hours"))
        ))
      }

    } else {  # actigraphy
      req(input$actig_file)
      f  <- input$actig_file
      df <- tryCatch(readr::read_csv(f$datapath, show_col_types = FALSE),
                     error = function(e) { showNotification(paste("Parse error:", e$message), type = "error"); NULL })
      req(!is.null(df))

      names(df) <- tolower(names(df))
      if ("sleep_wake" %in% names(df)) df$state_col <- df$sleep_wake
      else if ("state"      %in% names(df)) df$state_col <- df$state
      else { showNotification("Column 'sleep_wake' or 'state' not found.", type = "error"); return(NULL) }

      df$dt <- lubridate::parse_date_time(as.character(df$datetime),
                                          orders = c("YmdHMS","YmdHM","Ymd"))
      df <- df[order(df$dt), ]
      sleep_flag <- as.integer(df$state_col == 0)
      d <- diff(sleep_flag)
      onset_idx  <- which(d ==  1) + 1
      offset_idx <- which(d == -1) + 1
      if (sleep_flag[1] == 1)           onset_idx  <- c(1L, onset_idx)
      if (sleep_flag[nrow(df)] == 1)    offset_idx <- c(offset_idx, nrow(df))
      k  <- min(length(onset_idx), length(offset_idx))
      t0 <- df$dt[onset_idx[1]]
      data.frame(
        onset_h  = as.numeric(difftime(df$dt[onset_idx[1:k]],  t0, units = "hours")),
        offset_h = as.numeric(difftime(df$dt[offset_idx[1:k]], t0, units = "hours"))
      )
    }
  })

  output$tbl_file_preview <- renderDT({
    df <- parsed_file_schedule(); req(!is.null(df))
    datatable(round(df, 3),
              options = list(pageLength = 10, dom = "tp"),
              rownames = FALSE,
              colnames = c("Onset (h from t=0)", "Offset (h from t=0)"),
              class = "compact stripe")
  })

  # ---- Active schedule (from whichever input mode is selected) ----
  active_schedule <- reactive({
    if (input$input_mode == "manual") rv$schedule_manual
    else parsed_file_schedule()
  })

  # ---- Collect parameters -----------------------------------------
  # Build the two published parameter sets, overridden by slider values.
  current_params <- reactive({
    v  <- if (is.null(input$model_version)) "2024fe" else input$model_version
    Pp <- mccauley_params("PVT", v)
    Pk <- mccauley_params("KSS", v)

    # Shared parameters
    for (nm in c("beta", "zeta", "eta_w", "nu", "gamma", "lambda",
                 "phi", "tau", "mu_s")) {
      val <- input[[nm]]
      if (!is.null(val) && is.finite(val)) { Pp[[nm]] <- val; Pk[[nm]] <- val }
    }
    # Metric-specific
    if (!is.null(input$pf_pvt))      Pp$p_f     <- input$pf_pvt
    if (!is.null(input$alpha_w_pvt)) Pp$alpha_w <- input$alpha_w_pvt
    if (!is.null(input$alpha_s_pvt)) Pp$alpha_s <- input$alpha_s_pvt
    if (!is.null(input$Wc_pvt))      Pp$W_c     <- input$Wc_pvt
    if (!is.null(input$mu_w_pvt))    Pp$mu_w    <- input$mu_w_pvt
    if (!is.null(input$xi_pvt))      Pp$xi      <- input$xi_pvt

    if (!is.null(input$pf_kss))      Pk$p_f     <- input$pf_kss
    if (!is.null(input$alpha_w_kss)) Pk$alpha_w <- input$alpha_w_kss
    if (!is.null(input$alpha_s_kss)) Pk$alpha_s <- input$alpha_s_kss
    if (!is.null(input$Wc_kss))      Pk$W_c     <- input$Wc_kss
    if (!is.null(input$mu_w_kss))    Pk$mu_w    <- input$mu_w_kss
    if (!is.null(input$xi_kss))      Pk$xi      <- input$xi_kss

    # Equation 9: eta_s = eta_w * W_c / (T - W_c)  -- recomputed, never a free slider
    Pp$eta_s <- Pp$eta_w * Pp$W_c / (Pp$T_day - Pp$W_c)
    Pk$eta_s <- Pk$eta_w * Pk$W_c / (Pk$T_day - Pk$W_c)

    list(PVT = Pp, KSS = Pk)
  })

  # Restore every slider to the published Table 2 value
  restore_published <- function() {
    v  <- if (is.null(input$model_version)) "2024fe" else input$model_version
    Pp <- mccauley_params("PVT", v); Pk <- mccauley_params("KSS", v)
    for (nm in c("beta", "zeta", "eta_w", "nu", "gamma", "lambda",
                 "phi", "tau", "mu_s"))
      updateSliderInput(session, nm, value = Pp[[nm]])
    updateSliderInput(session, "pf_pvt",      value = Pp$p_f)
    updateSliderInput(session, "alpha_w_pvt", value = Pp$alpha_w)
    updateSliderInput(session, "alpha_s_pvt", value = Pp$alpha_s)
    updateSliderInput(session, "Wc_pvt",      value = Pp$W_c)
    updateSliderInput(session, "mu_w_pvt",    value = Pp$mu_w)
    updateSliderInput(session, "xi_pvt",      value = Pp$xi)
    updateSliderInput(session, "pf_kss",      value = Pk$p_f)
    updateSliderInput(session, "alpha_w_kss", value = Pk$alpha_w)
    updateSliderInput(session, "alpha_s_kss", value = Pk$alpha_s)
    updateSliderInput(session, "Wc_kss",      value = Pk$W_c)
    updateSliderInput(session, "mu_w_kss",    value = Pk$mu_w)
    updateSliderInput(session, "xi_kss",      value = Pk$xi)
    invisible(TRUE)
  }

  observeEvent(input$btn_restore_pub, {
    restore_published()
    showNotification("Published values restored.", type = "message")
  })

  # Switching version resets every slider to that version's published values
  observeEvent(input$model_version, {
    restore_published()
    showNotification(
      sprintf("Loaded published values for: %s",
              names(model_versions)[model_versions == input$model_version]),
      type = "message", duration = 5)
  }, ignoreInit = TRUE)

  output$version_note <- renderUI({
    v <- if (is.null(input$model_version)) "2024fe" else input$model_version
    txt <- switch(v,
      "2024fe" = paste0("Frontiers in Environmental Health 3:1362755 (2024). ",
        "Sleep inertia h(t) active; reciprocal feedback zeta = 1.31 active; ",
        "eta_w and mu_w refitted relative to 2021 (eta_w = 0.0126, ",
        "mu_w = 0.466 PVT / 0.82 KSS)."),
      "2024jtb" = paste0("Journal of Theoretical Biology 590:111851 (2024). ",
        "Sleep inertia h(t) active (nu and gamma estimated there); NO reciprocal ",
        "feedback; eta_w = 0.0074 and mu_w = 0.33 held at their 2021 values."),
      "2021" = paste0("Transportation Research Part F 79:94-106 (2021). ",
        "No sleep inertia and no reciprocal feedback: h(t) is held at zero, ",
        "which reduces the equations exactly to that paper's Eqs. 1-5."))
    helpText(tags$b("Selected version: "), txt)
  })

  # Flag whether the current parameters still equal the published ones
  output$using_published <- reactive({
    cp  <- tryCatch(current_params(), error = function(e) NULL)
    if (is.null(cp)) return(TRUE)
    v   <- if (is.null(input$model_version)) "2024fe" else input$model_version
    ref <- list(PVT = mccauley_params("PVT", v), KSS = mccauley_params("KSS", v))
    same <- TRUE
    for (m in c("PVT", "KSS"))
      for (nm in names(ref[[m]]))
        if (is.numeric(ref[[m]][[nm]]) &&
            !isTRUE(all.equal(ref[[m]][[nm]], cp[[m]][[nm]], tolerance = 1e-8)))
          same <- FALSE
    same
  })
  outputOptions(output, "using_published", suspendWhenHidden = FALSE)

  output$weight_check <- renderText({
    s <- input$w_kss + input$w_pvt
    if (s <= 0) "Weights must not both be zero."
    else sprintf("Weights are renormalised to sum to 1 (currently %.2f + %.2f).",
                 input$w_kss / s, input$w_pvt / s)
  })

  # ============================================================
  # LIVE SIMULATION  -- recalculates automatically as inputs change
  # ============================================================
  sim_inputs_live <- reactive({
    list(
      sched   = tryCatch(active_schedule(), error = function(e) NULL),
      params  = current_params(),
      dt      = as.numeric(input$dt_choice),
      t0c     = as.numeric(input$t0_clock),
      b_on    = as.numeric(input$base_onset),
      b_off   = as.numeric(input$base_offset),
      equil   = isTRUE(input$equilibrate),
      risk_w  = c(input$w_kss, input$w_pvt),
      pvt_ceil = as.numeric(input$pvt_ceiling)
    )
  })
  sim_inputs_debounced <- debounce(sim_inputs_live, millis = 500)

  sim_results <- reactive({
    if (isTRUE(input$auto_recalc)) {
      inp <- sim_inputs_debounced()
    } else {
      input$btn_run                       # take a dependency on the button
      inp <- isolate(sim_inputs_live())   # but read current inputs only on click
    }
    sched <- inp$sched
    if (is.null(sched) || nrow(sched) == 0) return(NULL)
    tryCatch(
      simulate_fatigue(sched, dt = inp$dt, t0_clock = inp$t0c,
                       base_onset = inp$b_on, base_offset = inp$b_off,
                       equilibrate = inp$equil,
                       P_pvt = inp$params$PVT, P_kss = inp$params$KSS,
                       risk_w = inp$risk_w, pvt_ceiling = inp$pvt_ceil),
      error = function(e) {
        showNotification(paste("Simulation error:", e$message), type = "error")
        NULL
      }
    )
  })

  # Keep the step recorder in sync with whatever is currently simulated
  # (fires silently on every recalculation; only the exported script content changes).
  observeEvent(sim_results(), {
    res <- sim_results()
    req(!is.null(res))
    sched <- isolate(active_schedule())
    p     <- isolate(current_params())
    dt    <- isolate(as.numeric(input$dt_choice))
    t0c   <- isolate(as.numeric(input$t0_clock))

    rec$record("schedule",
      paste0("Sleep schedule: ", nrow(sched), " period(s)"),
      sprintf("schedule_h <- data.frame(\n  onset_h  = c(%s),\n  offset_h = c(%s)\n)",
              paste(round(sched$onset_h,  2), collapse = ", "),
              paste(round(sched$offset_h, 2), collapse = ", "))
    )
    rec$record("params",
      "Model parameters (McCauley et al. 2024 parameterisation)",
      paste0("P_pvt <- list(\n",
             paste(mapply(function(k, v)
               sprintf("  %s = %s", k, if (is.numeric(v)) signif(v, 6) else shQuote(v)),
               names(p$PVT), p$PVT), collapse = ",\n"),
             "\n)\n\nP_kss <- list(\n",
             paste(mapply(function(k, v)
               sprintf("  %s = %s", k, if (is.numeric(v)) signif(v, 6) else shQuote(v)),
               names(p$KSS), p$KSS), collapse = ",\n"),
             "\n)")
    )
    rec$record("simulation",
      sprintf("Simulation run (dt = %.4f h, %d time points, %.1f h total)",
              dt, nrow(res), max(res$time_h)),
      sprintf(paste0(
        "# Paste mc2024_step(), mc2024_equilibrate(), simulate_mc2024()\n",
        "# and simulate_fatigue() from FaRiMaTo.R here\n",
        "results <- simulate_fatigue(schedule_h,\n",
        "                            dt = %.4f, t0_clock = %.1f,\n",
        "                            P_pvt = P_pvt, P_kss = P_kss)"),
        dt, t0c)
    )
  }, ignoreNULL = TRUE)

  output$sim_status_txt <- renderText({
    res <- sim_results()
    if (is.null(res)) {
      "No valid schedule yet -- add at least one sleep period."
    } else {
      mode_lbl <- if (isTRUE(input$auto_recalc)) "Live" else "Manual"
      sprintf("%s -- %d time points, %.1f h simulated.",
              mode_lbl, nrow(res), max(res$time_h))
    }
  })

  output$sim_done <- reactive(!is.null(sim_results()))
  outputOptions(output, "sim_done", suspendWhenHidden = FALSE)

  # ---- Value boxes ------------------------------------------------
  output$vb_kss <- renderValueBox({
    res <- sim_results(); req(!is.null(res))
    val <- round(mean(res$KSS[res$wake == 1], na.rm = TRUE), 2)
    col <- if (val < 4) "green" else if (val < 6) "yellow" else "red"
    valueBox(val, "Mean KSS (during wake)", icon = icon("moon"), color = col)
  })
  output$vb_pvt <- renderValueBox({
    res <- sim_results(); req(!is.null(res))
    val <- round(max(res$PVT[res$wake == 1], na.rm = TRUE), 1)
    col <- if (val < 5) "green" else if (val < 10) "yellow" else "red"
    valueBox(val, "Peak PVT lapses / 10 min", icon = icon("chart-bar"), color = col)
  })
  output$vb_highrisk <- renderValueBox({
    res <- sim_results(); req(!is.null(res))
    dt  <- diff(res$time_h[1:2])
    val <- round(sum(res$wake == 1 & res$Risk >= 0.66) * dt, 1)
    col <- if (val < 2) "green" else if (val < 6) "yellow" else "red"
    valueBox(paste0(val, " h"), "Waking hours in High Risk zone",
             icon = icon("exclamation-triangle"), color = col)
  })

  # ---- Plot helpers -----------------------------------------------
  make_xaxis <- function() list(title = "Time (h from study start)", zeroline = FALSE)

  # ============================================================
  # EMPIRICAL DATA IMPORTERS  (House Rules 1, 3, 6, 8)
  # ============================================================
  kss_long_data_raw <- make_wide_outcome_importer(
    input, output, session, rec, rv,
    prefix = "kss", label = "KSS", col_regex = "^kss[0-9]+_?$", default_step_h = 8)

  kss_long_data <- reactive({
    d <- kss_long_data_raw()
    req(!is.null(d))
    d$KSS <- d$value
    d
  })

  vas_long_data_raw <- make_wide_outcome_importer(
    input, output, session, rec, rv,
    prefix = "vas", label = "VAS", col_regex = "^vas[0-9]+_?$", default_step_h = 8)

  # Rescale raw VAS (0-100 by default) onto the KSS range (1-9)
  vas_kss_equiv <- reactive({
    d <- vas_long_data_raw()
    req(!is.null(d) && nrow(d) > 0)
    lo <- input$vas_min; hi <- input$vas_max
    req(is.finite(lo), is.finite(hi), hi > lo)
    frac <- (d$value - lo) / (hi - lo)
    frac <- pmin(pmax(frac, 0), 1)
    if (identical(input$vas_direction, "reverse")) frac <- 1 - frac
    d$VAS_raw <- d$value
    d$KSS     <- 1 + frac * 8
    d
  })

  pvt_long_data_raw <- make_wide_outcome_importer(
    input, output, session, rec, rv,
    prefix = "pvt", label = "PVT", col_regex = "^pvt[0-9]+_?$", default_step_h = 8)

  # PVT is lapses-only: the model's PVT arm predicts lapse counts, so this is
  # the single metric with a directly comparable model quantity.
  pvt_emp <- reactive({
    d <- pvt_long_data_raw()
    req(!is.null(d) && nrow(d) > 0)
    d$PVT <- d$value
    d
  })

  # ============================================================
  # KSS + PVT MAIN PLOT  (model predictions + all empirical overlays)
  # ============================================================
  output$plt_kss_pvt <- renderPlotly({
    res <- sim_results(); req(!is.null(res))
    wake <- res[res$wake == 1L, ]

    p <- plot_ly(wake, x = ~time_h) %>%
      add_lines(y = ~KSS, name = "KSS model prediction", yaxis = "y",
                line = list(color = "#2980b9", width = 2.5)) %>%
      add_lines(y = ~PVT, name = "PVT model prediction", yaxis = "y2",
                line = list(color = "#c0392b", width = 2, dash = "dash"))

    kss_emp <- tryCatch(kss_long_data(), error = function(e) NULL)
    p <- add_outcome_overlay(p, kss_emp, "KSS", yaxis = "y",
                             show_indiv = isTRUE(input$kss_show_indiv),
                             show_mean  = isTRUE(input$kss_show_mean),
                             label_prefix = "KSS observed",
                             marker_symbol = "circle", y_clip = c(1, 9))

    vas_emp <- tryCatch(vas_kss_equiv(), error = function(e) NULL)
    p <- add_outcome_overlay(p, vas_emp, "KSS", yaxis = "y",
                             show_indiv = isTRUE(input$vas_show_indiv),
                             show_mean  = isTRUE(input$vas_show_mean),
                             label_prefix = "VAS->KSS-equiv", line_dash = "dot",
                             marker_symbol = "triangle-up", y_clip = c(1, 9))

    pvt_emp_d <- tryCatch(pvt_emp(), error = function(e) NULL)
    p <- add_outcome_overlay(p, pvt_emp_d, "PVT", yaxis = "y2",
                             show_indiv = isTRUE(input$pvt_show_indiv),
                             show_mean  = isTRUE(input$pvt_show_mean),
                             label_prefix = "PVT observed",
                             marker_symbol = "square", y_clip = c(0, 30))

    p %>% layout(
      xaxis  = make_xaxis(),
      yaxis  = list(title = "KSS (1-9)", range = c(1, 9), tickvals = 1:9, gridcolor = "#eee"),
      yaxis2 = list(title = "PVT lapses", overlaying = "y", side = "right",
                    range = c(0, 22), zeroline = FALSE),
      shapes = sleep_shapes(res),
      legend = list(orientation = "h", y = -0.30),
      margin = list(r = 60, t = 10),
      plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe"
    )
  })

  # ---- Risk score plot (optional observed-sleepiness overlay) -------
  output$plt_risk <- renderPlotly({
    res <- sim_results(); req(!is.null(res))

    overlay   <- isTRUE(input$risk_overlay_sleep)
    on_risk   <- overlay && identical(input$risk_overlay_scale, "risk")
    # When mapping onto the risk scale the observed data share the primary
    # axis; in KSS units they get their own right-hand axis.
    ov_axis   <- if (on_risk) "y" else "y2"

    p <- plot_ly(res, x = ~time_h, y = ~Risk, name = "Composite risk (model)",
                 type = "scatter", mode = "lines",
                 line = list(color = "#8e44ad", width = 2)) %>%
      add_lines(y = rep(0.33, nrow(res)), name = "Low / Moderate",
                line = list(color = "#f39c12", dash = "dot", width = 1)) %>%
      add_lines(y = rep(0.66, nrow(res)), name = "Moderate / High",
                line = list(color = "#e74c3c", dash = "dot", width = 1))

    if (overlay) {
      # Convert observed sleepiness to whichever scale is displayed
      to_scale <- function(d) {
        if (on_risk) d$plotval <- (d$KSS - 1) / 8 else d$plotval <- d$KSS
        d
      }
      clip <- if (on_risk) c(0, 1) else c(1, 9)

      kss_emp <- tryCatch(kss_long_data(), error = function(e) NULL)
      if (!is.null(kss_emp) && nrow(kss_emp) > 0) {
        p <- add_outcome_overlay(p, to_scale(kss_emp), "plotval", yaxis = ov_axis,
                                 show_indiv = isTRUE(input$kss_show_indiv),
                                 show_mean  = isTRUE(input$kss_show_mean),
                                 label_prefix = "KSS observed",
                                 marker_symbol = "circle", y_clip = clip)
      }

      vas_emp <- tryCatch(vas_kss_equiv(), error = function(e) NULL)
      if (!is.null(vas_emp) && nrow(vas_emp) > 0) {
        p <- add_outcome_overlay(p, to_scale(vas_emp), "plotval", yaxis = ov_axis,
                                 show_indiv = isTRUE(input$vas_show_indiv),
                                 show_mean  = isTRUE(input$vas_show_mean),
                                 label_prefix = "VAS->KSS-equiv", line_dash = "dot",
                                 marker_symbol = "triangle-up", y_clip = clip)
      }
    }

    lay <- list(
      xaxis = make_xaxis(),
      yaxis = list(title = "Risk [0-1]", range = c(0, 1.02), gridcolor = "#eee"),
      shapes = c(
        list(
          list(type = "rect", layer = "below",
               fillcolor = "rgba(39,174,96,0.08)", line = list(width = 0),
               x0 = 0, x1 = max(res$time_h), y0 = 0,    y1 = 0.33),
          list(type = "rect", layer = "below",
               fillcolor = "rgba(243,156,18,0.10)", line = list(width = 0),
               x0 = 0, x1 = max(res$time_h), y0 = 0.33, y1 = 0.66),
          list(type = "rect", layer = "below",
               fillcolor = "rgba(231,76,60,0.10)", line = list(width = 0),
               x0 = 0, x1 = max(res$time_h), y0 = 0.66, y1 = 1.0)
        ),
        sleep_shapes(res)
      ),
      legend = list(orientation = "h", y = -0.35),
      margin = list(t = 10, r = if (overlay && !on_risk) 60 else 10),
      plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe"
    )
    if (overlay && !on_risk) {
      lay$yaxis2 <- list(title = "Observed sleepiness (KSS 1-9)",
                         overlaying = "y", side = "right",
                         range = c(1, 9), tickvals = 1:9, zeroline = FALSE)
    }

    do.call(plotly::layout, c(list(p), lay))
  })

  # ---- State variables plot ---------------------------------------
  output$plt_states <- renderPlotly({
    res <- sim_results(); req(!is.null(res))

    plot_ly(res, x = ~time_h) %>%
      add_lines(y = ~u_PVT, name = "u(t) allostatic [PVT arm]",
                line = list(color = "#8e44ad", width = 2)) %>%
      add_lines(y = ~u_KSS, name = "u(t) allostatic [KSS arm]",
                line = list(color = "#8e44ad", width = 2, dash = "dot")) %>%
      add_lines(y = ~h_PVT, name = "h(t) sleep inertia [PVT arm]",
                line = list(color = "#2980b9", width = 1.6)) %>%
      add_lines(y = ~h_KSS, name = "h(t) sleep inertia [KSS arm]",
                line = list(color = "#2980b9", width = 1.6, dash = "dot")) %>%
      add_lines(y = ~kappa, name = "kappa(t) circadian amplitude",
                line = list(color = "#27ae60", width = 2)) %>%
      add_lines(y = ~c, name = "c(t) circadian process",
                line = list(color = "#e67e22", width = 1.6, dash = "dash")) %>%
      layout(
        xaxis = make_xaxis(),
        yaxis = list(title = "State value", gridcolor = "#eee"),
        shapes = sleep_shapes(res),
        legend = list(orientation = "h", y = -0.35),
        margin = list(t = 10),
        plot_bgcolor = "#fefefe", paper_bgcolor = "#fefefe"
      )
  })

  # ---- Prediction table -------------------------------------------
  output$tbl_preds <- renderDT({
    res <- sim_results(); req(!is.null(res))
    dt_val <- diff(res$time_h[1:2])
    every  <- max(1L, round(1 / dt_val))  # ~hourly rows
    wake   <- res[res$wake == 1L, ]
    disp   <- wake[seq(1, nrow(wake), by = every), ] %>%
      transmute(
        `Time (h)`   = round(time_h,  2),
        `Clock (h)`  = round(clock_h, 2),
        `KSS`        = round(KSS,     2),
        `PVT lapses` = round(PVT,     2),
        `u (PVT)`    = round(u_PVT,   3),
        `u (KSS)`    = round(u_KSS,   3),
        `kappa`      = round(kappa,   3),
        `c(t)`       = round(c,       3),
        `Risk`       = round(Risk,    3)
      )
    datatable(disp,
              options = list(pageLength = 20, scrollX = TRUE, dom = "tp"),
              rownames = FALSE, class = "compact stripe") %>%
      formatStyle("Risk",
                  backgroundColor = styleInterval(
                    c(0.33, 0.66),
                    c("#a9dfbf", "#fdebd0", "#f5b7b1")))
  })

  # ---- CSV download -----------------------------------------------
  output$dl_csv <- downloadHandler(
    filename = function() paste0("farima_predictions_", Sys.Date(), ".csv"),
    content  = function(file) {
      res <- sim_results(); req(!is.null(res))
      write.csv(res, file, row.names = FALSE)
      rec$record("export_csv", "Full predictions exported to CSV",
                 'write.csv(results, "farima_predictions.csv", row.names = FALSE)')
    }
  )

  # ============================================================
  # MODEL vs EMPIRICAL COMPARISON  -- one factory call per outcome
  # ============================================================
  make_comparison_outputs(input, output, session, prefix = "kss", label = "KSS",
                          emp_reactive = kss_long_data, sim_results_fn = sim_results,
                          value_col = "KSS", model_col = "KSS", y_range = c(1, 9))

  make_comparison_outputs(input, output, session, prefix = "vas", label = "VAS (KSS-equiv)",
                          emp_reactive = vas_kss_equiv, sim_results_fn = sim_results,
                          value_col = "KSS", model_col = "KSS", y_range = c(1, 9))

  make_comparison_outputs(input, output, session, prefix = "pvtlap", label = "PVT lapses",
                          emp_reactive = pvt_emp, sim_results_fn = sim_results,
                          value_col = "PVT", model_col = "PVT", y_range = c(0, 30))

  # ---- Provenance tables -------------------------------------------
  output$tbl_prov_status <- renderDT({
    cp  <- current_params()
    v   <- if (is.null(input$model_version)) "2024fe" else input$model_version
    ref <- list(PVT = mccauley_params("PVT", v), KSS = mccauley_params("KSS", v))
    nms <- setdiff(names(ref$PVT), c("metric", "version", "use_inertia"))
    rows <- lapply(nms, function(nm) {
      pub_p <- ref$PVT[[nm]]; pub_k <- ref$KSS[[nm]]
      cur_p <- cp$PVT[[nm]];  cur_k <- cp$KSS[[nm]]
      if (!is.numeric(pub_p) || !is.numeric(cur_p)) return(NULL)
      changed <- !isTRUE(all.equal(pub_p, cur_p, tolerance = 1e-8)) ||
                 !isTRUE(all.equal(pub_k, cur_k, tolerance = 1e-8))
      data.frame(
        Parameter = nm,
        `Published PVT` = signif(pub_p, 6),
        `Published KSS` = signif(pub_k, 6),
        `In use PVT`    = signif(cur_p, 6),
        `In use KSS`    = signif(cur_k, 6),
        Status = if (changed) "CHANGED from published" else "published value",
        check.names = FALSE, stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    datatable(do.call(rbind, rows), rownames = FALSE, class = "compact stripe",
              options = list(pageLength = 25, dom = "t", scrollX = TRUE)) %>%
      formatStyle("Status",
                  color = styleEqual("CHANGED from published", "#c0392b"),
                  fontWeight = styleEqual("CHANGED from published", "bold"))
  })

  output$tbl_prov_published <- renderDT({
    datatable(published_reference_params, rownames = FALSE,
              colnames = c("Version", "Parameter", "PVT value", "KSS value"),
              filter = "top",
              class = "compact stripe",
              options = list(pageLength = 20, dom = "t", scrollX = TRUE))
  })

  # ---- R script / pipeline log export ------------------------------
  output$dl_rscript <- downloadHandler(
    filename = function() paste0("farima_script_", Sys.Date(), ".R"),
    content  = function(file) writeLines(rec$build_script(), file)
  )

  output$rscript_preview <- renderText({ rec$build_script() })
  output$pipeline_log    <- renderText({ rec$build_log()    })

}  # end server


# ============================================================
# 5. LAUNCH
# ============================================================
shinyApp(ui = ui, server = server)
