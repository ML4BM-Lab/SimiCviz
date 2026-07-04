# Due to warnings in build checks:
# Warning in utils::tar(filepath, pkgname, compression = compression, compression_level = 9L,  :
#    storing paths of more than 100 bytes is not portable
# The original file structure of SimiCPipeline is not mantained in the inst/exdata directory of the SimiCviz package. 
# Therefore the function `load_SimiCPipeline` does not work due to path inconsistency (but does work on real data runs), 
# so we need to generate the output as an RDS object to load in the vignette to make it work properly.

library(SimiCviz)
pickle_path = "./extdata/outputSimic/example_simic_weights.pickle"
auc_path = "./extdata/example_simic_auc_collected.csv"
annot_path = "./extdata/inputFiles/disease_stage_annotation.csv"

# Load objects
simic_weights <- read_weights_pickle(pickle_path)
output <- read_pickle(pickle_path)
adjusted_r_squared <- output$adjusted_r_squared
auc_csv <- read.csv(auc_path, row.names = 1, header = TRUE, sep = ",")
cell_labels <- load_cell_labels(annot_path, header = TRUE, sep = ",")

# Generate the experiment
simic_full <- SimiCvizExperiment(
  weights = simic_weights,
  auc = auc_csv,
  cell_labels = cell_labels,
  label_names = c("NBM","SMM","MM"),
  colors = c("#3B7EA1", "#E66101", "#B2182B"),
  meta = list(adjusted_r_squared = adjusted_r_squared)
  )
saveRDS(simic_full, file = "./extdata/simic_full.rds")
