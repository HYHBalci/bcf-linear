# run_all.R
# This script sequentially runs multiple R scripts and saves their output.

# Define the machine-independent output directory relative to the current working directory
out_dir <- file.path(getwd(), "run_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Run test_subgroup_sechidis.R
cat("Running test_subgroup_sechidis.R...\n")
sink(file.path(out_dir, "test_subgroup_sechidis_output.txt"))
env1 <- new.env()
source("R/test_subgroup_sechidis.R", local = env1)
sink()
save(list = ls(env1), envir = env1, file = file.path(out_dir, "test_subgroup_sechidis_workspace.RData"))

# 2. Run ACTG_complete.R
cat("Running ACTG_complete.R...\n")
plot_dir_actg <- file.path(out_dir, "ACTG_plots")
if (!dir.exists(plot_dir_actg)) dir.create(plot_dir_actg, recursive = TRUE)
sink(file.path(out_dir, "ACTG_complete_output.txt"))
env2 <- new.env()
env2$plot_dir <- plot_dir_actg
source("R/ACTG_complete.R", local = env2)
sink()
save(list = ls(env2), envir = env2, file = file.path(out_dir, "ACTG_complete_workspace.RData"))

# 3. Run ltr_test.R for PG
cat("Running ltr_test.R for PG...\n")
sink(file.path(out_dir, "ltr_test_PG_output.txt"))
env3 <- new.env()
env3$ltr_method <- "PG"
source("R/ltr_test.R", local = env3)
sink()
save(list = ls(env3), envir = env3, file = file.path(out_dir, "ltr_test_PG_workspace.RData"))

# 4. Run ltr_test.R for MSE
cat("Running ltr_test.R for MSE...\n")
sink(file.path(out_dir, "ltr_test_MSE_output.txt"))
env4 <- new.env()
env4$ltr_method <- "MSE"
source("R/ltr_test.R", local = env4)
sink()
save(list = ls(env4), envir = env4, file = file.path(out_dir, "ltr_test_MSE_workspace.RData"))

cat("All runs completed successfully. Outputs saved to:\n", out_dir, "\n")
