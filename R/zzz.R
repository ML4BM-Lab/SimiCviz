# .onLoad <- function(libname, pkgname) {
#   reticulate::py_require(
#     python_version = ">=3.8",
#     packages = c("numpy", "pandas")  # whatever you need
#   )
# }
.ensure_python <- function() {
  reticulate::py_require(
    python_version = ">=3.8",
    packages = c("numpy", "pandas", "anndata")
  )
}

#' @importFrom dplyr %>%
#' @importFrom graphics abline grid hist par plot.new title
#' @importFrom methods as new
#' @importFrom utils head read.csv
NULL