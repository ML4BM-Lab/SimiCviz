# Integration tests using bundled example data in inst/extdata/
# These tests exercise public API functions end-to-end with real data.

# ── Helpers ──────────────────────────────────────────────────────────────────

extdata <- function(...) system.file("extdata", ..., package = "SimiCviz")

# ── I/O ──────────────────────────────────────────────────────────────────────

test_that("read_weights_csv loads bundled weights correctly", {
  w <- SimiCviz::read_weights_csv(extdata("example_weights.csv"))
  expect_s3_class(w, "data.frame")
  expect_true(all(c("tf", "target", "weight", "label") %in% colnames(w)))
  expect_equal(nrow(w), 8000L)
  expect_equal(sort(unique(w$label)), c(0L, 1L, 2L, 3L))
})

test_that("load_cell_labels loads bundled annotation (case-insensitive, extra cols)", {
  cl <- SimiCviz::load_cell_labels(
    extdata("inputFiles", "treatment_annotation.csv"),
    header = TRUE, sep = ","
  )
  expect_s3_class(cl, "data.frame")
  expect_true(all(c("cell", "label") %in% colnames(cl)))
  expect_equal(nrow(cl), 3000L)
  expect_equal(sort(unique(cl$label)), 0L:3L)
})

test_that("read_auc_csv loads bundled AUC file (case-insensitive Cell column)", {
  a <- SimiCviz::read_auc_csv(extdata("example_auc.csv"))
  expect_s3_class(a, "data.frame")
  expect_true(all(c("cell", "tf", "score") %in% colnames(a)))
  expect_equal(nrow(a), 30000L)
})

test_that("load_from_csv builds SimiCvizExperiment from bundled CSV files", {
  obj <- SimiCviz::load_from_csv(
    weights_file     = extdata("example_weights.csv"),
    cell_labels_file = extdata("inputFiles", "treatment_annotation.csv")
  )
  expect_true(SimiCviz::is.SimiCvizExperiment(obj))
  expect_equal(length(obj@tf_ids), 10L)
  expect_equal(sort(unique(obj@cell_labels$label)), 0L:3L)
})

# ── simic_full RDS ────────────────────────────────────────────────────────────

test_that("simic_full.rds is a valid SimiCvizExperiment with AUC", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_true(SimiCviz::is.SimiCvizExperiment(simic_full))
  expect_equal(length(simic_full@tf_ids), 10L)
  expect_equal(sort(as.integer(names(simic_full@label_names))), 0L:3L)
  expect_false(is.null(simic_full@auc$collected))
})

# ── Dissimilarity ─────────────────────────────────────────────────────────────

test_that("calculate_dissimilarity returns a numeric matrix ranked by MinMax", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  dis <- SimiCviz::calculate_dissimilarity(simic_full)
  expect_true(is.data.frame(dis) || is.matrix(dis))
  expect_equal(nrow(dis), length(simic_full@tf_ids))
  expect_true(all(dis >= 0 & dis <= 1))
  # Scores should be in descending order
  scores <- dis[, 1]
  expect_true(all(diff(scores) <= 0))
})

test_that("calculate_dissimilarity with label subset returns smaller result", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  dis_sub <- SimiCviz::calculate_dissimilarity(simic_full, labels = c(0L, 1L))
  dis_all <- SimiCviz::calculate_dissimilarity(simic_full)
  # Same TFs, different scores
  expect_equal(nrow(dis_sub), nrow(dis_all))
  expect_false(identical(dis_sub[, 1], dis_all[, 1]))
})

# ── ECDF metrics ──────────────────────────────────────────────────────────────

test_that("calculate_ecdf_auc returns a data.frame with one row per TF", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  tfs <- simic_full@tf_ids[1:3]
  m <- SimiCviz::calculate_ecdf_auc(simic_full, tf_names = tfs)
  expect_s3_class(m, "data.frame")
  expect_equal(nrow(m), length(tfs))
})

# ── Network utilities ─────────────────────────────────────────────────────────

test_that("get_tf_network returns a data.frame in wide format (labels as columns)", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  net <- SimiCviz::get_tf_network(simic_full, simic_full@tf_ids[1])
  expect_s3_class(net, "data.frame")
  # Columns are label display names; rows are target genes
  lab_names <- unlist(simic_full@label_names, use.names = FALSE)
  expect_true(all(lab_names %in% colnames(net)))
  expect_gt(nrow(net), 0L)
})

# ── Plotting ──────────────────────────────────────────────────────────────────

test_that("plot_auc_heatmap runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_no_error(SimiCviz::plot_auc_heatmap(simic_full, top_n = 5))
})

test_that("plot_dissimilarity_heatmap runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_no_error(
    SimiCviz::plot_dissimilarity_heatmap(simic_full, top_n = 5, cmap = "viridis")
  )
})

test_that("plot_auc_distributions runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_no_error(
    SimiCviz::plot_auc_distributions(
      simic_full,
      tf_names  = simic_full@tf_ids[1:2],
      fill      = TRUE,
      alpha     = 0.6,
      bw_adjust = 1/8,
      grid      = c(1, 2)
    )
  )
})

test_that("plot_auc_cumulative runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  
  expect_no_error(
    SimiCviz::plot_auc_cumulative(
      simic_full,
      tf_names      = simic_full@tf_ids[1:2],
      grid          = c(1, 2),
      include_table = FALSE
    )
  )
})

test_that("plot_auc_summary_statistics runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_no_error(SimiCviz::plot_auc_summary_statistics(simic_full))
})

test_that("plot_tf_weights runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_no_error(
    SimiCviz::plot_tf_weights(simic_full, tf_names = simic_full@tf_ids[1:2],
                              top_n = 10, grid = c(1, 2))
  )
})

test_that("plot_target_weights runs without error on simic_full", {
  simic_full <- readRDS(extdata("simic_full.rds"))
  expect_no_error(
    SimiCviz::plot_target_weights(simic_full,
                                  target_names = simic_full@target_ids[1:2],
                                  grid = c(1, 2))
  )
})
