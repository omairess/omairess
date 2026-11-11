# ==============================================================================
# server/00a_help.R -- explanations behind an (i), status messages in plain view
# ==============================================================================
# Every explanatory note in the app is shown behind a small (i): clicking it
# opens the text in a pop-up, as the equation (i)s do, so the controls and the
# results are not buried under paragraphs. The browser moves the (i) next to
# the label of the control just above it (ui/00_theme.R: dance-help script).
#
#   dance_help(...)    an explanation: hidden, behind an (i). `title` names the
#                      pop-up (default: the label it sits next to)
#   dance_status(...)  a STATUS, warning or instruction that must be seen
#                      without a click (nothing selected yet, a test's method, a
#                      missing package): shown in plain view, as before
#   helpText(...)      every helpText() in the app is an explanation, so it is
#                      routed to dance_help(); this file is sourced into the UI's
#                      environment as well (ui/00_theme.R)
dance_help <- function(..., title = NULL) {
  shiny::tags$span(class = "dance-help", `data-title` = title,
    shiny::tags$a(href = "#", class = "dance-help-icon", title = "Explanation", role = "button",
                  onclick = "return danceHelp(this, event);",
                  shiny::icon("info-circle"), shiny::tags$span(class = "dance-help-more", "About this")),
    shiny::tags$div(class = "dance-help-body", style = "display:none", ...))
}
dance_status <- function(...) shiny::helpText(...)
helpText <- function(...) dance_help(...)
