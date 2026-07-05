# Inventory & Merge-Readiness Report: BootSON.R, Dagger_zero.R, PsychoNetrix.R

Repo: `/home/user/omairess` — three standalone R/Shiny apps. All three read completely (BootSON.R 7075 lines, Dagger_zero.R 6224 lines, PsychoNetrix.R 4163 lines).

---

## 1. Per-app structure, packages, and analysis pipeline

### 1.1 BootSON.R — bootnet GGM/RIN estimation + NSON (Nested-Specificity Network) viz

**Packages** (`BootSON.R:11-28`): shiny, shinyWidgets, DT, shinydashboard, colourpicker, dplyr, purrr, readr, vroom, readxl, tools, **bootnet, qgraph, mgm, relaimpo, igraph, Matrix, EGAnet, IsingFit, networktools, huge**, ggplot2, gridExtra, ggraph, plus GitHub-installed `NetworkComparisonTest` (development branch, `:21-28`). Package install is unconditional `install_if_missing()` (`:7-10`) — no version pinning/recording, no "no-auto-update" guard (violates house rule 0).

**UI** — single `tabsetPanel(id="tabs")` (`:1187-2228`) with tabs: **Data Input** (`:1190`), **Variables** (`:1242`), **Estimation & Visualization** (`:1367`, containing a nested `tabsetPanel` with Network/Centrality/Bridge Symptoms/Network Comparison/Data sub-tabs, `:1633-1863`), **NSON** (`:1868`, its own nested tabset: GNSS Table/Classic Network/NSON Plot/Edge Table/Conditional Probability Graph/Syndromic Depth & Centrality (k-way), `:2037-2224`).

**Server key reactives/observers**:
- File load: `raw_dat` reactive (`:2275-2348`) — csv/tsv/txt via vroom/readr, xlsx/xls via readxl; auto letter-suffix fix for BDI-style items (`:2322-2345`). Missing-data handling: `dat_unbalanced`/`dat` reactives implementing listwise/pairwise/keep policy and split-group equal-sampling with `set.seed(seed_value)` (`:2410-2412`).
- Variable typing: continuous/discrete/count buckets via `observeEvent` pairs (`:2555-2602`), `auto_detect` (`:2478-2519`).
- Estimation: `estimate_single_network()` (`:3294-3769`) dispatches to `bootnet::estimateNetwork()` for EBICglasso/pcor/cor/huge (`:3458,3475,3489,3516`), `IsingFit` (`:3614`), `mgm` (`:3674,3682`), `relimp` (`:3731,3741`). `est_result` eventReactive (`:2763-2833`) handles split-group looping.
- Bootstrap edge-significance: `run_boot_sig` observer calls `bootnet::bootnet(net, nBoots, type="nonparametric", statistics="edge")` (`:2712-2761`), then `apply_boot_threshold()` (`:2903-2945`) zeroes non-significant edges via empirical two-sided p-values + `p.adjust`.
- Bridge/case-dropping bootstrap: `bridge_stability_result` eventReactive calls `bootnet::bootnet(net, boots=..., type="case", statistics="bridgeExpectedInfluence", communities=...)` (`:6028-6056`) and prints `bootnet::corStability()` (CS-coefficient) at `:6083,6089` — **but this exists only for the Bridge-Symptoms tab**, not for the main Centrality tab.
- NCT: `run_nct` observer (`:5333-5593`) calls `NetworkComparisonTest::NCT()` (`:5532-5547`) with explicit `paired = is_paired` (`:5406-5471` determines paired vs independent via `input$nct_test_type`, `:1785-1788` UI).
- Plotting: `plot_single_network()` (`:4278-4780`) — one giant function driving `qgraph::qgraph()`/`plot.default` for both single and split-group cases; NSON directed/classic plots via `ggraph`/`igraph` (`:6357-6506, 6517-6598, 6600-6687`).
- Export: `downloadRDS` (`:5181-5271`, full bundle: raw data, classification, estimate, viz/centrality/bridge settings, session info) and `downloadRscript` (`:5273-5300`) driven by `generate_r_code()` (`:2974-3292`, a single monolithic string-builder invoked fresh on every download click).

### 1.2 Dagger_zero.R — bnlearn DAG structure learning ("The DAGger")

**Packages** (`:5-25`): shiny, bnlearn, qgraph, RColorBrewer, sortable, DT, shinycssloaders, shinyWidgets, plotly, igraph, visNetwork, shinyjs, jsonlite, colourpicker, readxl, haven. Same unconditional `install.packages()` pattern (`:9-10`), no version pinning.

**UI** — custom fixed sidebar/main-panel layout (not `sidebarLayout`) with collapsible sections (`:2780-3562`): Data Input panel (`create_data_panel()`, `:2436-2512`), Variables (drag-and-drop via `sortable::rank_list`, `:4005-4079`), Structure Learning (`create_algorithm_panel()`, `:2515-2599`), Network Constraints (blacklist/whitelist multi-select, `:2312-2379`), Split Analysis (`:2382-2429`), Bootstrap Analysis (`:2602-2622`), Visualization (`:2625-2774`) — then a results `tabsetPanel` (`:3328-3558`): Network Visualization, Network Statistics, Bootstrap Results, Model Comparison, Folded Temporal Graph, R Code Export, Save/Export.

**Server key reactives/observers**:
- File load: `read_data_file()` (`:161-178`) supports csv/txt/tsv/xlsx/xls/sav (haven), with haven-labelled → factor conversion. Long→wide via custom `long_to_wide()` (`:182-197`, merge-based, columns renamed `var_<time>`).
- Analysis core: `run_analysis()` (`:1713-1962`) — hill-climbing/tabu/mmhc/pc.stable/gs/iamb/rsmax2 initial fit (`:1770-1798`), then **`boot.strength(analysis_data, R=params$boot_r, algorithm=..., algorithm.args=...)`** (`:1851-1853`, default `bootR=500`, `:2606`), **`averaged.network(boot, threshold=params$threshold)`** (`:1861`), optional `cextend()` for full directionality (`:1864-1871`), `arc.strength()` (`:1877`) and `bnlearn::score()` (`:1903,1911`) on both original and averaged nets.
- Split (multi-group) analysis: loops `run_analysis()` per group level (`:4194-4272`), computes a shared `calculate_common_layout()` (`:288-406`, spring/circle/cascade with Katz centrality via `igraph::alpha_centrality`).
- Folded temporal graph: `detect_temporal_structure()` (`:1369-1440`), `infer_temporal_order()` (`:1466-1500`), `build_folded_edges()`/`plot_folded_temporal_graph()` (`:1503-1710`) — a unique, sophisticated feature specific to this app's panel-data use case.
- Plotting: `create_network_plot()` (qgraph static, `:972-1239`), `create_interactive_network_plot()` (visNetwork, `:547-920`), `create_dag_diff_plot()` (visNetwork group-vs-group diff, `:409-542`).
- Export: `generate_r_code()` (`:1965-2309`, monolithic, regenerated on demand), `downloadPlot` (PNG/PDF, `:5638-5699`), `exportAllResults` (ZIP of many CSVs, `:5464-5631`), `saveAnalysis`/`saveSettings` (RDS/JSON, `:5948-6028`).
- **No CS-coefficient / Markov-equivalence caveat anywhere** — UI copy repeatedly uses unqualified causal language: "gold standard for causal discovery" (`:2543`), "Excellent for identifying causal relationships" (`:2543`), "ensure known causal relationships are included" (`:2373`).

### 1.3 PsychoNetrix.R — psychonetrics GGM/Ising/CFA/LNM/RNM/LRNM/panel models

**Packages** (`:9-14`): shiny, shinydashboard, shinyWidgets, colourpicker, DT, **psychonetrics**, qgraph, dplyr, readxl, haven, ggplot2, viridis, Matrix, corpcor. Per-package `requireNamespace` install loop (`:16-21`) — again no version pinning/recording.

**UI** — `shinydashboard` with `sidebarMenu` (`:29-39`) and 7 `tabItem`s: 1. Data Import (`:58-117`), 2. Variable Setup (`:120-234`), 3. Model (`:237-400`, family radio: ggm/ising/lvm/lnm/rnm/lrnm/dlvm1/panelgvar/ri_clpm at `:242-254`), 4. Results (`:403-459`, `tabBox` with Fit Indices/Parameters/Modification Indices/Matrices/Model Summary/Compare Models/Factor Scores & Loadings), 5. Plots (`:462-587`), 6. Advanced (`:590-706`, model comparison + bootstrapping + R-code export), About (`:709-748`).

**Server key reactives/observers**:
- File load `observeEvent(list(input$file_upload,...))` (`:775-816`) — csv/txt/tsv/xlsx/xls/sav/por; haven-labelled→numeric conversion (`:797-800`). Long→wide via base R `reshape()` inside `analysis_data()` (`:882-897`), triggered by `radioButtons("data_format")` (`:71-74`) + `id_var`/`time_var` text inputs (`:78-79`).
- Data transform pipeline (nonparanormal/dichotomize-mean/median/log1p/sqrt/zscore) inside `analysis_data()` (`:914-980`) — unique to this app, none of the others transform data pre-estimation this way.
- Model fit: `observeEvent(input$run_model)` (`:1261-1486`) dispatches to `ggm()`, `Ising()`, `lvm()/lnm()/rnm()/lrnm()` (builder switch `:1387-1392`), `dlvm1()`, `panelgvar()`, `ri_clpm()` — all via `do_runmodel()` = `setoptimizer() %>% runmodel()` (`:1965-1967`), then optional `prune()`/`stepup`/`modelsearch` (`:1351-1353` etc.).
- Bootstrap: `run_boot` observer (`:3440-3498`) uses `loop_psychonetrics({...bootstrap="nonparametric"...}, reps=input$boot_reps, nCores=...)` + `aggregate_bootstraps()` — a **nonparametric parameter bootstrap**, not a case-dropping/CS-coefficient stability check; feeds `CIplot()`/manual CI plot (`ci_plot` type, `:2996-3076`).
- Model comparison: second model (`model2`) fit via near-duplicate code path (`:3220-3302`), compared via `psychonetrics:::compare` (`:3304,3336-3354`).
- Plotting: `make_plot()` (`:2177-3122`) — one function handling factor_network / factor_ri_network / latent_network / ri_latent_network / residual_network / ggm_network / ega_ggm / ega_ising / beta_network / within_network / between_network / cor_heat / ci_plot / semplot, each a `qgraph()` call with shared helper args (`qgraph_label_args()`, `qgraph_color_args()`, `pie_args()`).
- Export: `build_r_code()` (`:3510-4146`, monolithic) is generated **once**, only inside `input$run_model`'s success branch (`rv$r_code <- build_r_code()` at `:1481`) — i.e. it goes stale if the user later changes plot/appearance inputs, unlike BootSON/Dagger which regenerate on every download click.
- **No dedicated "export everything" bundle** — `save_model1` (`:3307-3314`) only `saveRDS()`s the bare fitted model object; no raw data / settings / plot-parameter bundle equivalent to BootSON's `downloadRDS` or Dagger's `saveAnalysis`/`exportAllResults`.

---

## 2. Candidates for a shared module (and where implementations diverge even though intent is identical)

| Concern | BootSON | Dagger_zero | PsychoNetrix | Divergence |
|---|---|---|---|---|
| **File dispatch** | csv/tsv/txt (vroom/readr) + xlsx/xls (readxl) only, `:1197-1198` | csv/txt/tsv/xlsx/xls/sav (haven), `:2455-2456,161-178` | csv/txt/tsv/xlsx/xls/sav/por, `:64-65,780-794` | Three independent `switch`/`tryCatch` implementations; **none support .rds/.RData raw-data import** (house rule 1 requires it); accepted-extension sets differ (BootSON lacks SPSS entirely) |
| **Wide↔long reshape** | Not implemented at all — always assumes wide | Custom `long_to_wide()` merge-based, output cols `var_<time>` (underscore), `:182-197`; reshaped data **is** shown in the existing preview table (rule-2 "inspectable" satisfied) | Base R `reshape(idvar=,timevar=,direction="wide")`, `:888-891`; output columns use R's default `.`-separated naming (`var.1`); **the pivoted `analysis_data()` is never rendered back to the user** — `data_table` shows only `rv$raw_data` (pre-pivot) — rule-2 "inspectable" violated | Column-naming convention differs (`var_1` vs `var.1`) — if unified, Dagger's own `detect_temporal_structure()` regex (`_t\d+$`/`_w\d+$`/`_\d+$`) would **not** recognize PsychoNetrix-style dot-separated wide columns |
| **Variable-selection UI** | Four parallel base `selectInput(..., multiple=TRUE, selectize=FALSE)` list-boxes with Add/Remove buttons (`:1259-1301`) — **not** a searchable pickerInput | `sortable::rank_list` drag-and-drop across three containers (`:4013-4078`) — **not** pickerInput either | `shinyWidgets::pickerInput(multiple=TRUE, liveSearch=TRUE, actionsBox=TRUE)` (`:130-136`) — the only one actually compliant with house rule 3 | Three structurally incompatible paradigms (dual-listbox / drag-drop / picker) — a real UI rewrite, not a trivial merge |
| **Missing-data policy** | `radioButtons("missing_policy", choices=listwise/pairwise/keep)`, `:1338-1342` | `radioButtons("missingDataMethod", choices=listwise/pairwise/none)`, `:2498-2505` | `selectInput("missing_method", choices=listwise/pairwise)`, `:168-171` | Same concept, 3 different input IDs/choice vocab/semantics (`keep` vs `none`) |
| **Group/split variable** | `selectInput("split_var")` (`:1305-1306`) | `selectInput("splitVar")` (`:2393`) | `selectizeInput("group_var")` (`:218-220`) | Same concept, 3 IDs; BootSON/Dagger cap at 10 groups (`:2774,4135`), PsychoNetrix has no explicit cap |
| **Node/border colour pickers** | `nodeColor` (#FFFFFF default, `:1519-1523`), `nodeBorderColor` (#8F8F8F, `:1524-1528`) | `nodeContinuousColor`(#e8f4f8)/`nodeDiscreteColor`(#fff3e0)/`nodeBorderColor`(#000000), `:2678-2692` | `node_color` (#72AFD3, `:521-522`)/`node_border_color` (#FFFFFF, `:523-524`) | **Exact ID collision** on `nodeBorderColor` between BootSON and Dagger with conflicting defaults (see §4). PsychoNetrix uses a snake_case variant with a third default. Three unrelated default colour choices |
| **Positive/negative edge colour** | Not exposed as a colourpicker anywhere — edge sign colour is baked into the `theme` parameter passed to `qgraph()` (`classic`/`colorblind`/`gray`/`Borkulo`, `:1552-1554`); NSON module hardcodes `theme_cols` per theme (`:6408-6414`) | Not exposed — arcs aren't signed (bnlearn structure, no partial-correlation sign) | Not exposed — relies on `qgraph`'s theme default pos/neg colouring | **Rule 4's pos/neg edge colourpicker is unimplemented in all three apps.** A shared module must also reconcile that Dagger's arcs are unsigned while BootSON/PsychoNetrix edges are signed |
| **Size sliders / plot export** | `vsize`/`esize`/`label_cex` numeric sliders exist (`:1577,1586,1562`) but **there is no PNG/PDF/SVG export of the network plot at all** — only `downloadRDS`/`downloadRscript` | `nodeSize`/`edgeWidth`/`labelSize` (`:2674,2734,2675`) + `downloadPlot` exporting **PNG (fixed 1200×800) or PDF** (`:5638-5699`) | `qgraph_vsize`/`qgraph_esize`/`qgraph_label_cex` (`:503,505,513`) + `download_plot` exporting **PDF only**, user-controlled width/height (`:3204-3217`) | Export format coverage: BootSON=none, Dagger=PNG/PDF, PsychoNetrix=PDF-only. **None support SVG.** Rule 5 is violated to varying degrees by all three |
| **R-script export mechanism** | `generate_r_code()` (`:2974-3292`), regenerated fresh **every time** `downloadRscript` fires (`:5285-5292`) | `generate_r_code()` (`:1965-2309`), regenerated fresh every time `downloadRCodeTab`/`rCodeDisplay` render (`:5260,5359`) | `build_r_code()` (`:3510-4146`), computed **once** at successful `run_model` and cached in `rv$r_code` (`:1481`) — can go stale if plot/appearance inputs change afterward | All three are monolithic string-builders, not an append-only step recorder. **None implement rule 6's id+description+code-fragment upsert pattern.** Freshness semantics differ (live vs cached) |
| **Pipeline log (plain-language, ordered)** | Absent — only ad hoc `message()`/console debug output (`debug_params()`, `:2955-2970`) and transient `showNotification()`s | Absent — only transient `showNotification()`s and `values$error_msg/success_msg` | Absent — only transient `showNotification()`s and a "transform banner" (`:1489-1500`) | **Rule 8 is unimplemented in all three.** Must be built from scratch as the shared recorder, sharing its underlying mechanism with rule 6 per the brief |
| **Export analysis data/inputs/results** | `downloadRDS` — full bundle (raw+prepped data, classification, estimate, viz/centrality/bridge/relimp settings, `session_info` w/ bootnet version), `:5181-5271` | `saveAnalysis` (RDS bundle, `:5948-5976`) + `exportAllResults` (ZIP of many CSVs incl. per-group, `:5464-5631`) — the most complete of the three | **No equivalent** — only bare-model `saveRDS()` (`save_model1`, `:3307-3314`) and separate loadings/factor-score CSVs; no combined data+settings+results export | PsychoNetrix is materially behind on rule 7; must be brought up to par, not just "merged" |
| **Package install/version recording (rule 0)** | Unconditional `install.packages()`, no version capture except inside `downloadRDS`'s `session_info$bootnet_version` (`:5252`) | Unconditional `install.packages()`, no version capture at all except raw `sessionInfo()` dump (`:5633-5635`, Save/Export tab) | Per-package `requireNamespace`+install, no version capture at all | None implement "record versions, no auto-update" as a first-class, visible requirement — needs a genuinely new shared installer/version-ledger module |

---

## 3. Genuinely unique per-app logic (must stay per-tab)

- **BootSON**: `estimate_single_network()` estimator dispatch (EBICglasso/pcor/cor/huge/IsingFit/mgm/relimp), bootstrap edge-significance thresholding (`apply_boot_threshold`), Bridge-Symptoms module (`networktools::bridge`, case-dropping bootstrap + `corStability`), the entire **NSON** sub-app (GNSS computation, conditional-probability graph, k-way syndromic depth — `compute_gnss`, `compute_kway_gnss`, `compute_cond_prob_matrix`, `:493-733`), NCT integration.
- **Dagger_zero**: bnlearn algorithm dispatch (hc/tabu/mmhc/pc.stable/gs/iamb/rsmax2) and 26 bnlearn score types, `boot.strength`/`averaged.network`/`cextend`, blacklist/whitelist constraint UI with multi-select cross-products (`:4140-4189`), DAG-difference plot between two split groups (`create_dag_diff_plot`), and the **Folded Temporal Graph** feature (`detect_temporal_structure`, `infer_temporal_order`, `build_folded_edges`, `plot_folded_temporal_graph`) — has no analogue in the other two apps.
- **PsychoNetrix**: the entire psychonetrics model-family dispatch (ggm/Ising/lvm/lnm/rnm/lrnm/dlvm1/panelgvar/ri_clpm), Lambda-matrix editor (`gen_lambda`/`lambda_editor_ui`/`simple_structure`), wave/beta-matrix panel-model UI (`detect_waves`, `beta_matrix_ui`), fit-index interpretive guide (`.fit_thresholds`/`.print_fit_guide`), Relative-Importance latent network via Johnson/LMG decomposition (`johnson_rw_from_cor`), factor-score computation (`fscores_data`, `fscores_merged_data`), EGA community detection with caching (`.ega_cache`), semPlot path diagrams, data-transformation pipeline (npn/dichotomize/log/sqrt/zscore), multi-group model comparison + invariance constraints (`groupEqual`).

---

## 4. Concrete conflicts (file:line cited)

### 4.1 Duplicate Shiny input IDs if naively concatenated

1. **`input$run`** — BootSON `actionButton("run", "Estimate Network", ...)` UI at `BootSON.R:1619`, consumed by `eventReactive(input$run, ...)` at `BootSON.R:2763`, and also reset by `observeEvent(input$run, {boot_sig_result(NULL)})` at `BootSON.R:2710`. Dagger_zero `actionButton("run", "Run Analysis", ...)` UI at `Dagger_zero.R:3314`, consumed by `observeEvent(input$run, {...})` at `Dagger_zero.R:4105`. **Exact string collision** — in a single merged `server()`, both `observeEvent`s would fire on the same button, running bootnet estimation and bnlearn structure learning simultaneously off one click.
2. **`input$nodeBorderColor`** — BootSON `colourpicker::colourInput("nodeBorderColor", ..., value = "#8F8F8F")` at `BootSON.R:1524-1528`; Dagger_zero `colourpicker::colourInput("nodeBorderColor", ..., value = "#000000")` at `Dagger_zero.R:2688-2692`. **Exact string collision with conflicting defaults** — whichever `colourInput` call is evaluated last in a merged UI silently overrides the other's default value/label.
3. **`input$layout`** — BootSON `shinyWidgets::pickerInput("layout", "Type", choices = c(spring, circle, groups, ega, pca))` at `BootSON.R:1537-1543`; Dagger_zero `selectInput("layout", "Layout:", choices = LAYOUTS)` where `LAYOUTS <- list("Spring"="spring","Circle"="circle","Cascade"="cascade")` at `Dagger_zero.R:2640` (`LAYOUTS` defined `Dagger_zero.R:127`). **Exact string collision, incompatible choice vocabularies** — `"cascade"` would crash BootSON's `plot_single_network()`; conversely `"ega"`/`"groups"`/`"pca"` selected in a merged app would break Dagger's `create_network_plot()`.

Naming-convention near-misses that are not exact collisions but would confuse a shared module and should be normalized during merge: `show_labels` (BootSON, `:1558`) vs `showLabels` (Dagger, `:2701`) vs `qgraph_labels_est` (PsychoNetrix, `:507`); `label_cex` (BootSON, `:1562`) vs `labelSize` (Dagger, `:2675`) vs `qgraph_label_cex` (PsychoNetrix, `:513`); `label_bold` (BootSON, `:1563`) vs `boldLabels` (Dagger, `:2671`); `nodeColor` (BootSON, `:1519`) vs `node_color` (PsychoNetrix, `:521`) vs `nodeContinuousColor`/`nodeDiscreteColor` (Dagger, `:2678,2683`); `downloadPlot` (Dagger, camelCase, `:3348`) vs `download_plot` (PsychoNetrix, snake_case, `:577`).

### 4.2 Wide/long format assumptions

- BootSON has **no** long-format support at all — always treats uploaded data as wide.
- Dagger_zero's `long_to_wide()` (`:182-197`) renames value columns to `<var>_<time>` (underscore-suffixed) and shows the pivoted result in the data preview.
- PsychoNetrix's long→wide pivot uses base `reshape(idvar=input$id_var, timevar=input$time_var, direction="wide")` (`:888-891`), which produces `<var>.<time>` (dot-suffixed) column names, and **never re-renders the pivoted result** back into `data_table` — violates house rule 2's "explicit inspectable" requirement.
- Consequence for a merge: Dagger's own `detect_temporal_structure()` (`:1369-1440`) pattern-matches `_t\d+$`/`_w\d+$`/`_\d+$` suffixes — it would **not** recognize PsychoNetrix-style `.`-suffixed wide columns, so a shared reshape module must standardize on one separator or the Folded-Temporal-Graph auto-detection silently fails post-merge.

### 4.3 Colour palette / pastel-default compliance (house rule 9)

None of the three apps uses `RColorBrewer::Pastel1`/`Pastel2` (confirmed via grep). Dagger_zero is the only one that even imports `RColorBrewer` (`:15`), but its edge-colour palette selector draws from the **sequential** category (`Blues` default, `:2739-2741`), not qualitative/pastel. BootSON's community-detection colouring uses `rainbow(n, s=0.4, v=0.95, alpha=0.85)` (`:4511`). PsychoNetrix and BootSON's NSON module both hardcode a fixed 10-colour hex array. **Net finding: rule 9 is violated by all three apps**, each with a different ad hoc palette.

### 4.4 Bootstrap / CS-coefficient / paired-NCT / boot.strength gating — explicit correctness-skill audit

- **bootnet (BootSON) — nonparametric edge-CI bootstrap**: implemented (`:2739-2744`), gated behind an opt-in checkbox `boot_sig_enable` (`:1594-1616`) — **off by default**.
- **bootnet (BootSON) — case-dropping bootstrap + CS-coefficient**: implemented **only for Bridge-Symptoms centrality** (`:6035-6038`, `corStability()` at `:6083,6089`). **The main "Centrality" tab has no case-dropping bootstrap and no CS-coefficient at all** — no CS<0.25 gate anywhere. Hard violation of the correctness skill.
- **Expected influence vs strength**: `ExpectedInfluence` is available (`:262-264`) but **default selection excludes it** (`c("Strength","Closeness","Betweenness")`, `:1699`) — contradicts the rule that EI must be shown whenever negative edges exist, with no auto-detection to force it.
- **corMethod pinning for EBICglasso**: always explicit (never silently defaulted at the API level, `:3458`), but its **UI default is `"cor"` (plain Pearson), not `"cor_auto"`** (`:1411-1417`), contradicting the ordinal-data rule given the app types Likert items as "Continuous."
- **EBIC gamma**: correctly, explicitly pinned (`tuning=0.5` default, always passed, `:3103,3459`) — compliant.
- **NCT paired handling**: compliant — explicit `paired=` radio + subject-matching logic (`:1785-1788,5409-5466,5538`).
- **NCT multiple-comparison correction**: **not implemented** — raw p<0.05 flagging with no `p.adjust()` anywhere in the NCT pathway (`:5541-5542,5721`). Direct violation.
- **bnlearn (Dagger_zero) — boot.strength + averaged.network + threshold**: correctly implemented, default R=500 (`:1851,1861,2606`). Compliant on mechanics.
- **bnlearn — Markov-equivalence / CPDAG caveat**: **absent** — unqualified causal language throughout UI copy (`:45,2543,2373`). Direct violation.
- **psychonetrics (PsychoNetrix) — estimator pinning**: defaults to `"default"` (psychonetrics' own ML→FIML auto-resolution) — borderline non-compliant in spirit.
- **psychonetrics — ordinal/DWLS mismatch**: only warns about DWLS instability if DWLS is already selected with ordered data; does **not** warn if estimator stays at "default"/ML while data_type="ordered" — the exact failure mode the correctness skill targets.
- **psychonetrics — bootstrap**: a nonparametric *parameter* bootstrap (CI-focused), not case-dropping/CS-coefficient — a different (and for SEM-family models, appropriate) concept; doesn't map onto the bootnet rule directly.

### 4.5 Step-recorder / R-script-export mechanism differences (critical for Stage 2)

All three implement the same anti-pattern: one big string-building function invoked at export time, not an append-only `{id, description, code}` recorder.
- **Regeneration timing**: BootSON/Dagger regenerate fresh on every download/view. PsychoNetrix caches once at `run_model` (`:1481`) — goes stale if plot/appearance/model2/bootstrap inputs change afterward.
- **Granularity**: PsychoNetrix's internal decomposition (`fit_code`/`prune_code`/`extract_code`/`fscores_code`/`plots_code`/`compare_code`/`boot_code`, `:3653-4145`) is the closest existing precedent for a step-recorder design.
- **No pipeline-log UI (rule 8) exists in any app** — must be built from scratch as shared infrastructure.

---

## 5. Proposed shared-module vs. per-tab architecture

**Shared modules:**
1. Data ingestion module — unified file dispatch for .txt/.csv/.xls/.xlsx/.sav/.rds/.RData (adding the two everyone is missing), uniform head+dims preview showing both raw and post-reshape frames.
2. Wide/long reshape module — one column-naming convention (recommend underscore, matching Dagger's existing temporal-detection regex), explicit previewed step.
3. Variable-selection module — standardize on `pickerInput(multiple=TRUE, liveSearch=TRUE, actionsBox=TRUE)` (PsychoNetrix's existing pattern, the only rule-3-compliant one).
4. Appearance module — node/border/**new pos-edge/neg-edge** colour pickers, unified size-slider naming, shared PNG/PDF/**SVG** export handler.
5. Step-recorder module (rules 6+8, unified) — ordered `{id, description, code}` records, upserted per stage; one renderer for the pipeline log, one for the R-script concatenation.
6. Results/export bundle module (rule 7) — one "export everything" handler extended to cover PsychoNetrix, which currently lacks one.
7. Package-install/version-ledger module (rule 0) — guarded installer recording `packageVersion()` for every dependency, never silently upgrading.

**Must stay per-tab:** BootSON's estimator dispatch, bootstrap-CI thresholding, Bridge-Symptoms module, NSON sub-app, NCT wiring. Dagger_zero's bnlearn algorithm/score dispatch, blacklist/whitelist builder, DAG-diff plot, Folded-Temporal-Graph. PsychoNetrix's model-family dispatch, Lambda/beta-matrix editors, fit-index guide, RI/Johnson-LMG decomposition, factor scores, EGA caching, semPlot, data-transform pipeline. The three stability concepts (bootnet edge-CI+CS-coefficient; bnlearn boot.strength/averaged.network/CPDAG; psychonetrics parameter bootstrap) are not interchangeable and must remain separate, model-family-specific gates.

**Decisions requiring careful reasoning to merge correctly (Stage 2 focus):**

1. Resolve the three confirmed exact ID collisions (`run`, `nodeBorderColor`, `layout`) via proper Shiny module namespacing (`ns()`), not just renaming — `layout`'s choice-vocabulary incompatibility is semantic, not just cosmetic.
2. Pick one wide/long naming convention before wiring Dagger's temporal auto-detection to a shared reshape module.
3. Add the missing CS-coefficient hard gate to BootSON's main Centrality tab (currently only Bridge has it) — new behavior, not something to "merge in."
4. Fix BootSON's EBICglasso `corMethod` default (`"cor"` → `"cor_auto"`) and add multiple-comparison correction to the NCT edge-level test — correctness bugs to fix during the merge.
5. Add an explicit CPDAG/Markov-equivalence caveat to Dagger_zero's results and scrub unqualified "causal discovery" UI copy before folding into shared results display.
6. Reconcile "estimator/corMethod/EBIC-gamma must never be left at package default" across bootnet (mostly compliant, one bug), bnlearn (compliant), and psychonetrics (`"default"` estimator, borderline) via a shared pinned-parameters audit panel.
7. Decide the palette-unification strategy for rule 9 given Dagger's edges are unsigned (strength/direction) vs BootSON/PsychoNetrix's signed GGM edges — "positive/negative edge colour" can't apply uniformly to bnlearn arcs.
8. Design the step-recorder's stage granularity to represent all three apps' differing branch structures as independently-upsertable records, fixing PsychoNetrix's stale-cache bug as part of that redesign.
9. Bring PsychoNetrix's export bundle up to parity with BootSON/Dagger before exposing one shared "Export" button across all three tabs.
