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
  "05-factor-adjustments-clz.r",
  "06-literature-reconciliation-hlz.r",
  "07-hlz-simpler.r"
)

simulations <- c(
  "91-run-bootstraps.r",
  "92-bootstrap-validation.r",
  "93-simulations-theory-free.r",
  "94-simulations-publication-bias.r",
  "95-simulations-standard-errors.r"
)

run_scripts(c(data_prep, main_analysis, simulations))

cat("\n=== All scripts completed successfully ===\n")
