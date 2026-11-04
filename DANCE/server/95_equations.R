# ==============================================================================
# server/95_equations.R -- the (i) buttons: one observer for the whole app
#
# Every (i) link (dance_eq_link(), server/08n_helpers_equations.R) sets
# input$dance_eq to the key of its quantity; this opens the formula, written for
# the model the current cosinor run fitted.
# ==============================================================================
observeEvent(input$dance_eq, {
  key <- input$dance_eq$key %||% ""
  if (!nzchar(key)) return()
  eq <- dance_equation_fixed(key) %||% dance_equation(key, dance_eq_context(values$harmonic_model))
  showModal(modalDialog(title = eq$title, HTML(eq$html), easyClose = TRUE, size = "m",
                        footer = modalButton("Close")))
})
