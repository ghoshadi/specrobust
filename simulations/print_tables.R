# Tables of the full results: oracle n = 1e7; 1000 replications, 1000 bootstrap
# draws, n = 1000; read from results/. Run from this folder.
source("helper.R")

for (eg in c("1", "2")) {
  tauR <- as.numeric(readLines(file.path("results", sprintf("tauR_eg%s.txt", eg))))
  res  <- as.matrix(read.csv(file.path("results", sprintf("Eg%s_1000sims_1000boot.csv", eg))))
  print_table(res, eg, tauR)
  cat(readLines(file.path("results", sprintf("one_run_eg%s.txt", eg))), sep = "\n")
}
