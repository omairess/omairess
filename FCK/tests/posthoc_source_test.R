# ==============================================================================
# tests/posthoc_source_test.R — the post-hoc tests must compare what they claim
#
# The bug this pins: with two or more scalar variables selected, the omnibus
# functional ANOVA ran on input$fanova_group_var while the post-hoc tests ran on
# values$group_labels — the FIRST scalar variable selected at import. Two
# different variables, no warning. The failure mode is invisible in the output,
# which is exactly why it needs a test rather than a look.
#
# Run with:   Rscript tests/posthoc_source_test.R      (from the FCK directory)
# Needs no packages: fck_posthoc_spec() is pure logic over lists.
# ==============================================================================

app_dir <- if (dir.exists("server")) "." else "FCK"
eval(parse(file.path(app_dir, "server/06_helpers_posthoc.R"), encoding = "UTF-8"))
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a)) b else a

failures <- 0
check <- function(label, cond, detail = "") {
  if (isTRUE(cond)) cat(sprintf("ok  : %s\n", label))
  else { failures <<- failures + 1; cat(sprintf("FAIL: %s  %s\n", label, detail)) }
}

# a stand-in fd object: only ncol(coefs) is read, as the curve count
fake_fd <- function(n) list(coefs = matrix(0, nrow = 5, ncol = n))

n <- 12
# Sex is selected FIRST, so it is values$group_labels — the old code's choice.
# The fANOVA was run on Condition, the SECOND variable.
sex       <- factor(rep(c("F", "M"), length.out = n))
condition <- factor(rep(c("Ctrl", "Sleep", "Recov"), each = 4))
subject   <- factor(rep(paste0("P", 1:4), times = 3))

values <- list(
  data = matrix(0, n, 20),
  fd_obj = fake_fd(n),
  group_labels = sex,                       # first scalar variable
  selected_group_vars = c("Sex", "Condition", "Subject"),
  group_variables = data.frame(Sex = sex, Condition = condition, Subject = subject),
  covariates = data.frame(Sex = sex, Condition = condition, Subject = subject),
  fanova_results = list(design = "between", group_var = "Condition",
                        group_labels = condition, fd_used = fake_fd(n))
)

# --- 1. following the omnibus picks the omnibus's variable, not the first one -
spec <- fck_posthoc_spec(list(posthoc_source = "fanova"), values)
check("follows the fANOVA variable", isTRUE(spec$ok) &&
      identical(levels(spec$group_labels), levels(condition)),
      paste("got", paste(levels(spec$group_labels), collapse = "/")))
check("does NOT fall back to the first scalar variable",
      !identical(levels(spec$group_labels), levels(sex)))
check("says so in the description", grepl("Condition", spec$description, fixed = TRUE))
check("marked as matching the omnibus", isTRUE(spec$matches_omnibus))

# --- 2. a level subset in the omnibus carries over ---------------------------
sub_idx <- which(condition != "Recov")
values2 <- values
values2$fanova_results$group_labels <- droplevels(condition[sub_idx])
values2$fanova_results$fd_used <- fake_fd(length(sub_idx))
spec2 <- fck_posthoc_spec(list(posthoc_source = "fanova"), values2)
check("inherits the omnibus level subset", isTRUE(spec2$ok) &&
      identical(levels(spec2$group_labels), c("Ctrl", "Sleep")))

# --- 3. choosing a variable here overrides, and is flagged as a new family ---
spec3 <- fck_posthoc_spec(list(posthoc_source = "custom", posthoc_design = "between",
                               posthoc_group_var = "Sex"), values)
check("custom choice is honoured", isTRUE(spec3$ok) &&
      identical(levels(spec3$group_labels), levels(sex)))
check("custom choice is NOT called post-hoc", isFALSE(spec3$matches_omnibus))

# --- 4. within-subjects, chosen here, without an RM omnibus ------------------
spec4 <- fck_posthoc_spec(list(posthoc_source = "custom", posthoc_design = "within",
                               posthoc_subject_var = "Subject",
                               posthoc_rm_var = "Condition"), values)
check("paired design available after a between-subjects omnibus",
      isTRUE(spec4$ok) && identical(spec4$design, "within"))
check("paired design reports its pairing",
      isTRUE(spec4$ok) && grepl("4 of 4 subjects", spec4$description))

# --- 5. refusals, rather than silently comparing the wrong thing -------------
one_level <- values; one_level$group_variables$Flat <- factor(rep("only", n))
one_level$selected_group_vars <- c(one_level$selected_group_vars, "Flat")
spec5 <- fck_posthoc_spec(list(posthoc_source = "custom", posthoc_design = "between",
                               posthoc_group_var = "Flat"), one_level)
check("refuses a single-level variable", isFALSE(spec5$ok) &&
      grepl("only 1 level", spec5$message))

mismatch <- values; mismatch$fd_obj <- fake_fd(n + 3)
spec6 <- fck_posthoc_spec(list(posthoc_source = "custom", posthoc_design = "between",
                               posthoc_group_var = "Sex"), mismatch)
check("refuses labels that do not line up with the curves",
      isFALSE(spec6$ok) && grepl("do not line up|values but there are", spec6$message))

no_pairs <- values
no_pairs$covariates$Subject <- factor(paste0("P", 1:n))   # each subject once
spec7 <- fck_posthoc_spec(list(posthoc_source = "custom", posthoc_design = "within",
                               posthoc_subject_var = "Subject",
                               posthoc_rm_var = "Condition"), no_pairs)
check("refuses a paired test with nothing paired",
      isFALSE(spec7$ok) && grepl("more than one level", spec7$message))

spec8 <- fck_posthoc_spec(list(posthoc_source = "fanova"), list(fanova_results = NULL))
check("refuses before the omnibus has run", isFALSE(spec8$ok))

# --- 6. RM columns come from the row-aligned frame --------------------------
v <- list(data = matrix(0, 4, 10),
          covariates = data.frame(ID = c("a", "b", "c", "d")),
          uploaded_data = data.frame(ID = c("a", "b", "c", "d", "e", "f")))  # unfiltered
check("prefers the row-aligned covariates",
      identical(as.character(fck_rm_column(v, "ID")), c("a", "b", "c", "d")))

v2 <- list(data = matrix(0, 4, 10), covariates = NULL,
           uploaded_data = data.frame(Other = c("a", "b", "c", "d", "e", "f")))
err <- tryCatch({ fck_rm_column(v2, "Other"); NA_character_ },
                error = function(e) conditionMessage(e))
check("errors instead of mispairing when only the raw frame has it",
      !is.na(err) && grepl("no longer lines up", err))

if (failures) { cat("\n", failures, " failure(s).\n", sep = ""); quit(status = 1) }
cat("\nPost-hoc source tests passed.\n")
