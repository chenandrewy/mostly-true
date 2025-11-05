# ABOUTME: Coordinates the full analysis workflow by sourcing project scripts in sequence.
# Inputs:
#   - None directly; each sourced script manages its own dependencies.
# Outputs:
#   - Aggregated side effects from sourced scripts (tables, figures, data artifacts inside `results/` and `data/`).
# How to run:
#   """shell
#   Rscript main.r
#   """

library(here)
here::i_am("main.r")

run_scripts <- function(scripts) {
  for (script in scripts) {
    cat(sprintf("\n=== Running %s ===\n", script))
    source(here(script), chdir = FALSE)
  }
}

data_prep <- c(
  "01-prep-data.r"
)

main_analysis <- c(
  "02-easy-bounds.r",
  "03-visual-bounds.r",
  "04-factor-adjustments-yz.r",
  "04alt-factor-adjustments-clz.r",
  "05-hlz-simpler.r"
)

run_scripts(c(data_prep))
run_scripts(c(main_analysis))

cat("\n=== All scripts completed successfully ===\n")
