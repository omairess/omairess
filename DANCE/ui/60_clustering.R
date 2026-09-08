# ==========================================================================
# ui/60_clustering.R
#
# PORTED VERBATIM by tools/port_fck.py — do not hand-edit the ranges
# below without updating that script's manifest.  Provenance:
#   WaPaa1_3.R lines 792-1065  (Functional clustering)
# ==========================================================================
ui_tab_clustering <- tabItem(
        tabName = "kmeans",
        fluidRow(
          box(
            title = "Optimal Number of Clusters",
            collapsible = TRUE, collapsed = FALSE,
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            h4("Find Optimal k"),
            p("Run this analysis first to determine the optimal number of clusters before performing the final clustering."),
            fluidRow(
              column(3,
                numericInput("max_clusters_test", "Maximum k to test:",
                             value = 10, min = 5, max = 20, step = 1),
                helpText("Tests k from 2 to this value")
              ),
              column(3,
                radioButtons("opt_data_type", "Data to use:",
                             choices = list("Raw Data" = "raw",
                                            "Smoothed Data" = "smoothed"),
                             selected = "smoothed")
              ),
              column(3,
                radioButtons("opt_kmeans_method", "Method:",
                             choices = list("Standard" = "standard",
                                            "Functional" = "functional",
                                            "Hierarchical" = "hierarchical"),
                             selected = "standard"),
                conditionalPanel(
                  condition = "input.opt_kmeans_method == 'functional' && input.opt_data_type == 'raw'",
                  div(style = "color: orange; font-weight: bold; font-size: 11px;",
                      "⚠ Functional requires smoothed data")
                ),
                conditionalPanel(
                  condition = "input.opt_kmeans_method == 'hierarchical'",
                  selectInput("opt_hclust_linkage", "Linkage:",
                              choices = list("Ward.D2" = "ward.D2",
                                             "Complete" = "complete",
                                             "Average" = "average",
                                             "Single" = "single"),
                              selected = "ward.D2")
                )
              ),
              column(3,
                br(),
                actionButton("run_cluster_optimization", "Analyze Optimal k",
                             class = "btn-warning", icon = icon("chart-line"))
              ),
              column(3,
                br(),
                verbatimTextOutput("optimization_status")
              )
            ),
            hr(),
            fluidRow(
              column(6,
                plotlyOutput("elbow_plot", height = "350px")
              ),
              column(6,
                plotlyOutput("silhouette_plot", height = "350px")
              )
            ),
            hr(),
            verbatimTextOutput("optimization_recommendation"),
            conditionalPanel(
              condition = "input.opt_kmeans_method == 'hierarchical'",
              hr(),
              h4("Dendrogram"),
              p("Shows the full agglomerative tree. Use it together with the elbow and silhouette plots to choose the optimal cut point (k)."),
              plotOutput("opt_dendrogram_plot", height = "400px")
            )
          )
        ),
        fluidRow(
          box(
            title = "Clustering Settings",
            collapsible = TRUE, collapsed = FALSE,
            status = "primary",
            solidHeader = TRUE,
            width = 4,
            numericInput("n_clusters", "Number of Clusters (k):",
                         value = 3, min = 2, max = 20, step = 1),
            radioButtons("cluster_data_type", "Data to Cluster:",
                         choices = list("Raw Data" = "raw",
                                        "Smoothed Data" = "smoothed"),
                         selected = "smoothed"),
            helpText("Choose whether to cluster raw or smoothed functional data."),
            hr(),
            radioButtons("cluster_standardize", "Standardize Before Clustering:",
                         choices = list("None" = "none",
                                        "Within participants (person-mean centering)" = "within",
                                        "Between participants (z-score per time point)" = "between"),
                         selected = "none"),
            helpText(HTML("<b>None:</b> Use data as-is.<br>
                          <b>Within participants:</b> Center and scale each participant's own time series (removes between-person differences in level and variability).<br>
                          <b>Between participants:</b> Z-score each time point across participants (standardizes each occasion to mean 0, SD 1).")),
            hr(),
            radioButtons("clustering_method", "Clustering Method:",
                         choices = list("Standard (Point-wise) K-Means" = "standard",
                                        "Functional K-Means" = "functional",
                                        "DCF (Density Core Finding)" = "dcf",
                                        "Hierarchical Clustering" = "hierarchical"),
                         selected = "standard"),
            helpText(HTML("<b>Standard K-Means:</b> Clusters based on values at each time point.<br>
                          <b>Functional K-Means:</b> Uses functional distances (fda.usc, requires smoothed data).<br>
                          <b>DCF:</b> Density-based clustering - automatically finds number of clusters (requires Python).<br>
                          <b>Hierarchical:</b> Agglomerative clustering with chosen linkage; cuts the tree at k clusters.")),
            conditionalPanel(
              condition = "input.clustering_method == 'functional' && input.cluster_data_type == 'raw'",
              div(style = "color: orange; font-weight: bold;",
                  "⚠ Functional K-Means requires smoothed data. Please select 'Smoothed Data' above.")
            ),
            # Hierarchical-specific parameters
            conditionalPanel(
              condition = "input.clustering_method == 'hierarchical'",
              wellPanel(
                style = "background-color: #f5fff5;",
                h5("Hierarchical Parameters"),
                selectInput("hclust_linkage", "Linkage Method:",
                            choices = list("Ward.D2" = "ward.D2",
                                           "Complete" = "complete",
                                           "Average (UPGMA)" = "average",
                                           "Single" = "single"),
                            selected = "ward.D2"),
                helpText("Ward.D2 minimises total within-cluster variance and is recommended for compact clusters.")
              )
            ),
            # DCF-specific parameters
            conditionalPanel(
              condition = "input.clustering_method == 'dcf'",
              wellPanel(
                style = "background-color: #f0f8ff;",
                h5("DCF Parameters"),
                numericInput("dcf_k", "k (neighborhood size):",
                             value = 10, min = 2, max = 100, step = 1),
                helpText("Number of neighbors for density estimation. Larger k = smoother density."),
                numericInput("dcf_beta", "β (fluctuation parameter):",
                             value = 0.1, min = 0.01, max = 1, step = 0.01),
                helpText("Controls density variation within cluster cores. Smaller β = tighter cores."),
                hr(),
                p("DCF automatically determines the number of clusters.", style = "font-style: italic; color: #666;"),
                actionButton("check_dcf_setup", "Check Python/DCF Setup", class = "btn-sm btn-info")
              )
            ),
            # Partitioning method parameters (Standard and Functional K-Means only)
            conditionalPanel(
              condition = "input.clustering_method == 'standard' || input.clustering_method == 'functional'",
              hr(),
              numericInput("kmeans_nstart", "Number of Random Restarts:",
                           value = 10, min = 1, max = 100, step = 1),
              helpText("The algorithm is run this many times with different random initialisations; the best solution (lowest within-cluster SS) is kept."),
              numericInput("kmeans_iter", "Maximum Iterations per Restart:",
                           value = 100, min = 10, max = 1000, step = 10)
            ),
            hr(),
            actionButton("run_clustering", "Run Clustering",
                         class = "btn-primary btn-lg", icon = icon("play")),
            br(), br(),
            verbatimTextOutput("clustering_status")
          ),
          box(
            title = "Cluster Statistics",
            collapsible = TRUE, collapsed = FALSE,
            status = "success",
            solidHeader = TRUE,
            width = 8,
            h4("Cluster Summary"),
            DTOutput("cluster_summary_table"),
            hr(),
            h4("Fit Statistics"),
            verbatimTextOutput("clustering_fit_stats")
          )
        ),
        fluidRow(
          box(
            title = "Cluster Membership",
            collapsible = TRUE, collapsed = FALSE,
            status = "info",
            solidHeader = TRUE,
            width = 12,
            DTOutput("cluster_membership_table"),
            hr(),
            h4("Save Cluster Membership to Dataset"),
            fluidRow(
              column(3,
                textInput("cluster_var_name", "New Variable Name:",
                          value = "Cluster", placeholder = "e.g., Cluster_K3")
              ),
              column(3,
                br(),
                actionButton("save_cluster_membership", "Add to Dataset",
                             class = "btn-success", icon = icon("plus"))
              ),
              column(3,
                br(),
                downloadButton("download_data_with_clusters", "Download Dataset",
                               class = "btn-primary")
              ),
              column(3,
                br(),
                verbatimTextOutput("cluster_save_status")
              )
            ),
            helpText("Click 'Add to Dataset' to add cluster membership to your in-memory data (usable in other modules). Click 'Download Dataset' to export the full dataset including cluster membership to a CSV file.")
          )
        ),
        fluidRow(
          box(
            title = "Group Composition by Cluster",
            collapsible = TRUE, collapsed = FALSE,
            status = "success",
            solidHeader = TRUE,
            width = 12,
            conditionalPanel(
              condition = "output.has_group_labels",
              uiOutput("cluster_group_var_selector"),
              h4("Group Counts by Cluster"),
              DTOutput("cluster_group_counts_table"),
              hr(),
              h4("Group Percentages by Cluster (%)"),
              DTOutput("cluster_group_pct_table"),
              hr(),
              verbatimTextOutput("cluster_group_test")
            ),
            conditionalPanel(
              condition = "!output.has_group_labels",
              p("No group variable available. Group composition analysis requires a group variable in your data.")
            )
          )
        ),
        fluidRow(
          box(
            title = "Cluster Mean Functions",
            collapsible = TRUE, collapsed = FALSE,
            status = "warning",
            solidHeader = TRUE,
            width = 12,
            selectInput("tick_freq_kmeans", "X-axis tick frequency:",
                        choices = list("Every 15 min" = 0.25, "Every 30 min" = 0.5,
                                       "Every hour" = 1, "Every 2 hours" = 2,
                                       "Every 4 hours" = 4, "All" = 0),
                        selected = 1, width = "200px"),
            plotlyOutput("cluster_means_plot", height = "500px")
          )
        ),
        fluidRow(
          box(
            title = "Individual Curves by Cluster",
            collapsible = TRUE, collapsed = FALSE,
            status = "info",
            solidHeader = TRUE,
            width = 12,
            plotlyOutput("cluster_individuals_plot", height = "600px")
          )
        ),
        fluidRow(
          box(
            title = "Detailed Silhouette Analysis",
            collapsible = TRUE, collapsed = FALSE,
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            p("Silhouette plot showing individual sample coefficients. Each bar represents one subject."),
            p("Values close to 1 indicate well-clustered samples, values near 0 indicate borderline cases, negative values indicate potential misclassifications."),
            plotlyOutput("detailed_silhouette_plot", height = "600px")
          )
        ),
        conditionalPanel(
          condition = "input.clustering_method == 'hierarchical'",
          fluidRow(
            box(
              title = "Dendrogram",
              collapsible = TRUE, collapsed = FALSE,
              status = "success",
              solidHeader = TRUE,
              width = 12,
              p("Agglomerative dendrogram. The red dashed line shows the cut point for the selected k."),
              plotOutput("dendrogram_plot", height = "450px")
            )
          )
        )
      )
