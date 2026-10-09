# A small version of the two simulation examples, outputs saved to demo_run/.
# The full results (oracle n = 1e7; 1000 replications, 1000 bootstrap draws,
# n = 1000) are in results/ and print_tables.R prints them. Run from this folder.
rm(list = ls())
library(specrobust)
source("helper.R")

n_large <- as.numeric(Sys.getenv("N_LARGE", 1e5))
n_sims  <- as.integer(Sys.getenv("N_SIMS",  50))
n_boot  <- as.integer(Sys.getenv("N_BOOT",  100))
out.dir <- Sys.getenv("OUT_DIR", "demo_run")
cores   <- max(1L, parallel::detectCores() - 1L)

dir.create(out.dir, showWarnings = FALSE)
for (eg in c("1", "2")) {
  tauR <- run_oracle(eg, n_large, out.dir)
  res  <- run_eval(eg, tauR, n_sims, n_boot, cores = cores)
  write.csv(res, file.path(out.dir, sprintf("Eg%s_%dsims_%dboot.csv", eg, n_sims, n_boot)),
            row.names = FALSE)
  print_table(res, eg, tauR)
  cat(readLines(file.path(out.dir, sprintf("one_run_eg%s.txt", eg))), sep = "\n")
}
