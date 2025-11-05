# ABOUTME: Coordinates the full analysis workflow by sourcing project scripts in sequence.
# Inputs:
#   - None directly; each sourced script manages its own dependencies.
# Outputs:
#   - pdf and tex files in results/
# How to run:
#   """shell
#   Rscript main.r
#   """

library(here)
here::i_am("main.r")

cat("\n=== Running 01-prep-data.r ===\n")
source(here("01-prep-data.r"), chdir = FALSE, echo = TRUE)

cat("\n=== Running 02-easy-bounds.r ===\n")
source(here("02-easy-bounds.r"), chdir = FALSE, echo = TRUE)

cat("\n=== Running 03-visual-bounds.r ===\n")
source(here("03-visual-bounds.r"), chdir = FALSE, echo = TRUE)

cat("\n=== Running 04-factor-adjustments-yz.r ===\n")
source(here("04-factor-adjustments-yz.r"), chdir = FALSE, echo = TRUE)

cat("\n=== Running 04alt-factor-adjustments-clz.r ===\n")
source(here("04alt-factor-adjustments-clz.r"), chdir = FALSE, echo = TRUE)

cat("\n=== Running 05-hlz-simpler.r ===\n")
source(here("05-hlz-simpler.r"), chdir = FALSE, echo = TRUE)

cat("\n=== All scripts completed successfully ===\n")
