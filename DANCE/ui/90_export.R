# ==============================================================================
# ui/90_export.R — merged Export tab
#
# Union of the two source apps' export tabs:
#   WaPaa1_3.R lines 1068-1113  (fPCA, fANOVA, post-hoc, clustering, smoothed
#                                data, plot bundles, and the R-code generator)
#   CIRCAREG.R lines  657-673   (FoSR coefficients, cosinor parameters and
#                                the cosinor summary report)
#
# Two ids had to move:
#   * CIRCAREG's "export_scores_csv" (FoSR beta coefficients) collided with
#     WaPaa's "export_scores_csv" (fPCA scores) -> export_fosr_coefs_csv.
#   * CIRCAREG's "Download Reproduction R Code" button is NOT here: its handler
#     wrote a single placeholder comment ("R Code Generation logic here
#     (omitted for brevity)"). WaPaa's generator below is real, so that one
#     button is the app's code export rather than shipping a stub beside it.
# ==============================================================================

ui_tab_export <- tabItem(
  tabName = "export",
  fluidRow(
    box(
      title = "Export Options",
      collapsible = TRUE, collapsed = FALSE,
      status = "info",
      solidHeader = TRUE,
      width = 12,

      h4("Functional PCA / warping"),
      downloadButton("export_scores_csv", "Download PCA Scores (CSV)", class = "btn-primary"),
      downloadButton("export_loadings_csv", "Download PCA Loadings (CSV)", class = "btn-primary"),
      conditionalPanel(
        condition = "input.pca_type == 'twpca'",
        downloadButton("export_warping_csv", "Download Warping Scores (CSV)", class = "btn-primary")
      ),
      hr(),

      h4("Functional ANOVA"),
      downloadButton("export_fanova_results_csv", "Download FANOVA Results (CSV)", class = "btn-primary"),
      downloadButton("export_pairwise_results_csv", "Download Post-hoc Results (CSV)", class = "btn-primary"),
      hr(),

      h4("Circadian / regression"),
      downloadButton("export_fosr_coefs_csv", "Download FoSR Coefficients (CSV)", class = "btn-primary"),
      downloadButton("export_harmonic_params", "Download Cosinor Individual Parameters (CSV)", class = "btn-primary"),
      downloadButton("export_harmonic_summary", "Download Cosinor Summary Report (TXT)", class = "btn-primary"),
      helpText("Cosinor pairwise results and their plot download from the",
               "'Cosinor: pairwise tests' tab itself."),
      hr(),

      h4("Clustering Results"),
      downloadButton("export_cluster_membership_csv", "Download Cluster Membership (CSV)", class = "btn-primary"),
      downloadButton("export_cluster_means_csv", "Download Cluster Mean Functions (CSV)", class = "btn-primary"),
      downloadButton("export_cluster_stats_csv", "Download Cluster Statistics (CSV)", class = "btn-primary"),
      downloadButton("export_silhouette_csv", "Download Silhouette Data (CSV)", class = "btn-primary"),
      conditionalPanel(
        condition = "output.has_group_labels",
        downloadButton("export_cluster_group_csv", "Download Cluster-Group Composition (CSV)", class = "btn-primary")
      ),
      hr(),

      h4("Export Smoothed Data"),
      helpText("The shared smoothed curves at original time resolution with original time labels —",
               "these are the curves every analysis in this app was run on."),
      downloadButton("export_smoothed_csv", "Download Smoothed Curves (CSV)", class = "btn-primary"),
      downloadButton("export_smoothed_long_csv", "Download Smoothed Curves - Long Format (CSV)", class = "btn-primary"),
      helpText("The matching observed/interpolated/extrapolated map downloads from",
               "the 'Missing data & filled points' panel on the smoothing tab — export",
               "both together, or the smoothed curves travel without any record of",
               "which of their values were measured."),
      hr(),

      h4("Export Plots"),
      downloadButton("export_plots", "Download All Plots (PDF)", class = "btn-success"),
      downloadButton("export_fanova_plots", "Download FANOVA Plots (PDF)", class = "btn-success"),
      downloadButton("export_pairwise_plots", "Download Post-hoc Plots (PDF)", class = "btn-success"),
      hr(),

      # ---- the APA report -----------------------------------------------
      # Every other export on this tab writes NUMBERS. This one writes the
      # sentences a paper needs: which analytic choices were made, and what
      # the results were, in APA 7th-edition style. It reports only analyses
      # that were actually run, and each Results subsection closes with what
      # its numbers do not establish -- the corrections this app has had to
      # make are worth nothing if the document produced for publication
      # reintroduces the overstatement.
      h4("Publication report (APA 7)"),
      textInput("apa_report_title", "Report title (optional):",
                placeholder = "Functional data analysis: methods and results",
                width = "100%"),
      downloadButton("download_apa_html", "Download Report (HTML, print to PDF)",
                     class = "btn-success"),
      downloadButton("download_apa_md", "Download Report (Markdown)",
                     class = "btn-primary"),
      helpText(HTML(paste(
        "A Methods and Results document written from the analyses in this",
        "session. The <b>Statistical analysis</b> section states the choices",
        "&mdash; basis and smoothing parameter and how it was selected,",
        "registration method and its boundary rule, estimator, multiplicity",
        "correction, permutation count and the smallest <i>p</i> it can",
        "produce. The <b>Results</b> section reports the numbers in APA style",
        "(exact <i>df</i>, <i>p</i> to three decimals with no leading zero,",
        "effect sizes named for what they are), and each subsection ends with",
        "what those numbers do <i>not</i> establish. Open the HTML in a browser",
        "and print to PDF for a formatted copy; the Markdown is for pasting",
        "into a manuscript.<br><br>",
        "<b>Check every number against the app before submission.</b> This is",
        "an aid to writing, not a substitute for reading your own results."))),
      hr(),

      h4("Report preview"),
      div(style = "max-height:26em;overflow:auto;background:#f8f9fa;padding:.6em;border:1px solid #dee2e6;border-radius:4px;",
          verbatimTextOutput("apa_report_preview")),
      hr(),

      h4("Export R Code"),
      downloadButton("export_code", "Download Analysis Code (R)", class = "btn-warning"),
      helpText("Reproduces the whole pipeline in plain R: import, smoothing,",
               "functional PCA/ANOVA/clustering, and the cosinor and FoSR",
               "fits. The cosinor section carries the app's own fitting function",
               "verbatim rather than a re-implementation of it."),
      hr(),

      h4("Analysis Code Preview"),
      verbatimTextOutput("code_preview")
    )
  ),

  # Neither source app could be reopened: close the browser and the import,
  # the smoothing and every fitted model were gone.
  fluidRow(
    box(
      title = "Session", status = "success", solidHeader = TRUE, width = 12,
      collapsible = TRUE, collapsed = FALSE,
      h4("Save this session"),
      helpText("Writes one .rds holding the data, the variable selection, the",
               "smoothed curves and every result you have run, plus the exact",
               "package versions they were produced with."),
      downloadButton("save_session", "Save session (.rds)", class = "btn-success"),
      hr(),
      h4("Restore a session"),
      fileInput("load_session", "Open a saved .rds", accept = ".rds"),
      verbatimTextOutput("session_status")
    )
  )
)
