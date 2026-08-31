# ==============================================================================
# tests/smoke_test.R — structural smoke test for the merged app
#
# Checks the things a merge can break that a human would only notice by
# clicking every tab:
#   1. every ui/*.R and server/*.R file parses;
#   2. the whole UI assembles and renders to HTML (so every tabItem, every
#      conditionalPanel and every output placeholder is well formed);
#   3. every sidebar menuItem points at a tabItem that exists, and no two
#      tabItems claim the same tabName;
#   4. the entire server body registers under a mock session — i.e. all 16
#      server files source into one environment without an id clash, a missing
#      helper, or an outputOptions() call on an output that was never defined;
#   5. no output id is assigned twice (in a single server environment the
#      second assignment silently wins).
#
# It does NOT run any statistics: it needs only the UI packages (shiny,
# shinydashboard, DT, plotly), not
# fda/refund/rmfanova, so it can run anywhere.
#
# Run with:   Rscript tests/smoke_test.R      (from the FCK directory)
# ==============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(DT)       # DTOutput / dataTableOutput placeholders in the UI
  library(plotly)   # plotlyOutput placeholders in the UI
})

app_dir <- if (dir.exists("ui")) "." else "FCK"
stopifnot(dir.exists(file.path(app_dir, "ui")),
          dir.exists(file.path(app_dir, "server")))

ui_files     <- sort(list.files(file.path(app_dir, "ui"), full.names = TRUE, pattern = "[.]R$"))
server_files <- sort(list.files(file.path(app_dir, "server"), full.names = TRUE, pattern = "[.]R$"))
all_files    <- c(file.path(app_dir, "app.R"), ui_files, server_files)

fail <- function(...) { cat("FAIL:", ..., "\n"); quit(status = 1) }
ok   <- function(...) cat("ok  :", ..., "\n")

# -- 1. everything parses ------------------------------------------------------
for (f in all_files) {
  e <- tryCatch({ parse(f, encoding = "UTF-8"); NULL },
                error = function(e) conditionMessage(e))
  if (!is.null(e)) fail(f, "does not parse:", e)
}
ok(length(all_files), "files parse")

# -- 2. UI assembles and renders ----------------------------------------------
fck_source <- function(file, envir = parent.frame()) {
  eval(parse(file, encoding = "UTF-8"), envir = envir)
  invisible(NULL)
}
for (f in ui_files) fck_source(f, envir = globalenv())

app_src   <- paste(readLines(file.path(app_dir, "app.R"), warn = FALSE), collapse = "\n")
ui_names  <- unique(regmatches(app_src, gregexpr("ui_tab_[A-Za-z0-9_]+", app_src))[[1]])
missing   <- ui_names[!vapply(ui_names, exists, logical(1))]
if (length(missing)) fail("app.R references undefined UI objects:", paste(missing, collapse = ", "))

ui_obj <- do.call(tabItems, lapply(ui_names, get))
html   <- htmltools::renderTags(ui_obj)$html
if (!nzchar(html)) fail("UI rendered to empty HTML")
ok("UI assembles and renders (", length(ui_names), "tabs,", nchar(html), "chars of HTML )")

# -- 3. sidebar <-> tabItem wiring --------------------------------------------
# Read the tab ids back out of the RENDERED html: shinydashboard turns
# tabItem(tabName = "x") into <div id="shiny-tab-x">, so this checks what the
# browser will actually receive rather than what the source says.
tab_names <- regmatches(
  html, gregexpr('(?<=id="shiny-tab-)[A-Za-z0-9_]+', html, perl = TRUE))[[1]]
menu_names <- unique(regmatches(
  app_src, gregexpr('(?<=tabName = ")[A-Za-z0-9_]+', app_src, perl = TRUE))[[1]])

dupes <- tab_names[duplicated(tab_names)]
if (length(dupes)) fail("two tabItems share a tabName (only one is reachable):",
                        paste(unique(dupes), collapse = ", "))
orphan_menu <- setdiff(menu_names, tab_names)
if (length(orphan_menu)) fail("sidebar entries with no tabItem:",
                              paste(orphan_menu, collapse = ", "))
orphan_tab <- setdiff(tab_names, menu_names)
if (length(orphan_tab)) fail("tabItems no sidebar entry can reach:",
                             paste(orphan_tab, collapse = ", "))
ok("all", length(tab_names), "tabs are uniquely named and reachable from the sidebar")

# -- 4 + 5. the server registers, and no output is defined twice ---------------
assigned <- character(0)
for (f in server_files) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
  txt <- gsub("(?m)#.*$", "", txt, perl = TRUE)
  assigned <- c(assigned, regmatches(
    txt, gregexpr("(?<=output\\$)[A-Za-z0-9_.]+(?=\\s*<-)", txt, perl = TRUE))[[1]])
}
dup_out <- unique(assigned[duplicated(assigned)])
if (length(dup_out)) fail("output assigned more than once:", paste(dup_out, collapse = ", "))
ok(length(assigned), "outputs, each assigned exactly once")

server <- function(input, output, session) {
  for (f in server_files) fck_source(f, envir = environment())
}

res <- tryCatch({
  testServer(shinyApp(ui_obj, server), { NULL })
  NULL
}, error = function(e) conditionMessage(e))
if (!is.null(res)) fail("server body failed to register:", res)
ok("all", length(server_files), "server files register in one environment")

cat("\nSmoke test passed.\n")
