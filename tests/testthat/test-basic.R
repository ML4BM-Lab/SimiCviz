test_that("SimiCvizExperiment constructor works", {
  w <- data.frame(
    tf = c("TF1", "TF2"),
    target = c("Gene1", "Gene2"),
    weight = c(0.5, 0.8),
    label = c(0, 1)
  )
  cell_labels <- data.frame(
    cell = c("Cell1", "Cell2"),
    label = c(0, 1)
  )
  expect_s4_class(SimiCviz::SimiCvizExperiment(weights = w, cell_labels = cell_labels), "SimiCvizExperiment")
})