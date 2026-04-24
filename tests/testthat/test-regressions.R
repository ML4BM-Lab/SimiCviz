test_that("calculate_activity_scores works without adj_r2 metadata", {
  weights_df <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    weight = c(0.8, 0.6),
    label = c(0, 1),
    stringsAsFactors = FALSE
  )

  cell_labels <- data.frame(
    cell = c("Cell1", "Cell2"),
    label = c(0, 1),
    stringsAsFactors = FALSE
  )

  expr <- matrix(
    c(
      5, 3,
      2, 4,
      1, 1,
      1, 1
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2", "TF1", "TF2"), c("Cell1", "Cell2"))
  )

  simic <- SimiCviz::SimiCvizExperiment(weights = weights_df, cell_labels = cell_labels)
  simic_out <- SimiCviz::calculate_activity_scores(simic, expression = expr, verbose = FALSE)

  expect_true(!is.null(simic_out@auc$collected))
  expect_equal(nrow(simic_out@auc$collected), 2)
})

test_that("get_auc default returns wide format", {
  weights_df <- data.frame(
    tf = c("TF1", "TF1"),
    target = c("Gene1", "Gene2"),
    weight = c(0.7, 0.4),
    label = c(0, 0),
    stringsAsFactors = FALSE
  )

  cell_labels <- data.frame(
    cell = c("Cell1", "Cell2"),
    label = c(0, 0),
    stringsAsFactors = FALSE
  )

  expr <- matrix(
    c(
      6, 2,
      1, 3,
      1, 1
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("Gene1", "Gene2", "TF1"), c("Cell1", "Cell2"))
  )

  proc <- SimiCviz::AUCProcessor(weights = weights_df, expression = expr, cell_labels = cell_labels)
  proc <- SimiCviz::compute_auc(proc, verbose = FALSE)

  auc_wide <- SimiCviz::get_auc(proc)
  expect_true(is.data.frame(auc_wide))
  expect_equal(nrow(auc_wide), 2)
  expect_true("TF1" %in% colnames(auc_wide))

  expect_error(SimiCviz::get_auc(proc, format = "matrix"), "format must be 'wide' or 'long'")
})

test_that("load_from_csv works without cell labels when auc is missing", {
  weights_df <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    weight = c(0.5, 0.9),
    label = c(0, 1),
    stringsAsFactors = FALSE
  )

  weights_file <- tempfile(fileext = ".csv")
  utils::write.csv(weights_df, weights_file, row.names = FALSE)

  obj <- SimiCviz::load_from_csv(weights_file = weights_file)

  expect_s4_class(obj, "SimiCvizExperiment")
  expect_equal(nrow(obj@cell_labels), 0)
})

test_that("plot_r2_distribution uses fallback colors", {
  adjusted_r_squared <- list("0" = stats::runif(10), "1" = stats::runif(12))

  out_dir <- tempdir()
  out_file <- tempfile(pattern = "r2_", fileext = ".pdf", tmpdir = out_dir)

  expect_invisible(
    SimiCviz::plot_r2_distribution(
      adjusted_r_squared = adjusted_r_squared,
      save = TRUE,
      filename = basename(out_file),
      out_dir = dirname(out_file)
    )
  )
  expect_true(file.exists(out_file))
})
