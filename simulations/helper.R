# Functions for the two simulation examples of Section 6.1.
# run_demo.R runs a small version; print_tables.R prints the full results.

adj_sets <- list("X1", c("X1", "X2"))

draw <- function(eg, n, seed = 123) {
  set.seed(seed)
  if (eg == "1") {
    X1 <- rnorm(n); X2 <- rnorm(n)
    A  <- as.numeric(runif(n) <= 1/(exp(5*X1 + 5*X2) + 1))
    Y  <- A*(1 + X1 - 5*X2) + 4*X2 + rnorm(n)
  } else {
    U1 <- rnorm(n); U2 <- rnorm(n); X1 <- rnorm(n)
    A  <- as.numeric(U1 + X1 > 0); X2 <- U1 + U2
    Y  <- ifelse(A == 1, 1 + 0.5*X1 - U2 + rnorm(n), 5*U2 + rnorm(n))
  }
  list(Y = Y, A = A, X = data.frame(X1, X2))
}

# the paper's overlay: 200 bins, both populations, and a smoothed outline
overlay <- function(data, weights, breaks = 200, hcol = "orange",
                    hcol2 = "skyblue2", xlim = NULL, ylim = NULL, xlab = NULL,
                    legend.text = c("original population", "new target population"),
                    legend.pos = "topleft") {
  h <- hist(data, breaks = breaks, plot = FALSE); bw <- diff(h$breaks)
  nb <- length(bw); i <- findInterval(data, h$breaks); k <- i >= 1L & i <= nb
  wc <- as.numeric(tapply(weights[k], factor(i[k], levels = seq_len(nb)), sum))
  wc[is.na(wc)] <- 0
  d <- list(h$counts/sum(h$counts)/bw, wc/sum(weights)/bw)
  if (is.null(xlim)) xlim <- range(data)
  if (is.null(ylim)) ylim <- c(0, 1.05*max(unlist(d)))
  plot(h$mids, d[[1]], type = "n", xlim = xlim, ylim = ylim, xlab = xlab,
       ylab = "Density", main = NULL)
  cols <- adjustcolor(c(hcol, hcol2), alpha.f = 0.6)
  legend(legend.pos, legend = legend.text, fill = cols, col = cols, bty = "n")
  for (j in 1:2)
    rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], d[[j]], col = cols[j],
         border = adjustcolor(c(hcol, hcol2)[j], alpha.f = 0.2))
  for (j in 1:2)
    lines(h$mids, stats::filter(d[[j]], rep(1/5, 5), sides = 2) + 5e-4,
          col = c(hcol, hcol2)[j], lwd = 3)
  box()
}

lims <- list(
  "1" = list(X1 = list(xlim = c(-4.5, 4.8), ylim = c(0, 0.45), pos = "topleft"),
             X2 = list(xlim = c(-5, 5),     ylim = c(0, 0.45), pos = "topleft")),
  "2" = list(X1 = list(xlim = c(-4.8, 4.8), ylim = c(0, 0.45), pos = "topright"),
             X2 = list(xlim = c(-7, 7),     ylim = c(0, 0.32), pos = "topleft")))

# tau_R and the covariate shift; writes tauR_eg<eg>.txt, one_run_eg<eg>.txt and
# the two figures into out_dir, and returns tau_R.
run_oracle <- function(eg, n_large, out_dir, n_small = 1000, alpha = 0.05) {
  wmean <- function(x, w) sum(w*x)/sum(w)
  wvar  <- function(x, w) sum(w*(x - wmean(x, w))^2)/sum(w)
  txt <- character(0); say <- function(...) txt <<- c(txt, sprintf(...))
  d <- draw(eg, n_small)
  s <- tryCatch(specrobust(d$Y, d$A, d$X, adj_sets, reg_mode = "lm", alpha = alpha,
                           n_boot = 1000, seed = 123), error = function(e) e)
  if (inherits(s, "error")) say("  one replication FAILED: %s", conditionMessage(s)) else {
    say("n = %d, one replication (seed 123)", n_small)
    say("  true ATE  tau_bar               : %.4f", 1)
    for (k in 1:2)
      say("  adj set %d  estimate of tau_bar   : %8.4f   CI [%8.4f, %8.4f]",
          k, s$candidate_estimates[k], s$candidate_ci[k, 1], s$candidate_ci[k, 2])
    say("  convex hull of the above        :            CI [%8.4f, %8.4f]", s$hull_ci[1], s$hull_ci[2])
    say("  proposed estimand  tau_R        : %8.4f   CI [%8.4f, %8.4f]", s$estimate, s$ci[1], s$ci[2])
    say("  lambda*                         : %8.4f", s$lambda)
  }

  d <- draw(eg, n_large)
  o <- tryCatch(specrobust(d$Y, d$A, d$X, adj_sets, reg_mode = "lm", n_boot = 0, seed = 123),
                error = function(e) e)
  if (inherits(o, "error")) {
    say("  oracle FAILED: %s", conditionMessage(o))
    writeLines(txt, file.path(out_dir, sprintf("one_run_eg%s.txt", eg))); return(NA_real_)
  }
  w <- o$weights; X1 <- d$X$X1; X2 <- d$X$X2
  say("")
  say("n = %g, oracle (seed 123)", n_large)
  say("  oracle tau_R                    : %8.4f", o$estimate)
  say("  lambda*                         : %8.4f", o$lambda)
  say("  reweighted estimates per set    : %s", paste(sprintf("%.4f", o$reweighted_estimates), collapse = ", "))
  say("  moments of the shared covariate X1 and of X2:")
  say("    X1  original   mean %8.4f  var %8.4f", mean(X1), var(X1))
  say("    X1  reweighted mean %8.4f  var %8.4f", wmean(X1, w), wvar(X1, w))
  say("    X2  original   mean %8.4f  var %8.4f", mean(X2), var(X2))
  say("    X2  reweighted mean %8.4f  var %8.4f", wmean(X2, w), wvar(X2, w))

  writeLines(sprintf("%.10f", o$estimate), file.path(out_dir, sprintf("tauR_eg%s.txt", eg)))
  writeLines(txt, file.path(out_dir, sprintf("one_run_eg%s.txt", eg)))
  for (v in c("X1", "X2")) {
    f <- file.path(out_dir, sprintf("oracle_shift_eg%s_%s.pdf", eg, v))
    pdf(f, width = 3.5*16/9, height = 3.5)
    par(mar = c(4.0, 4.2, 1.0, 1.0), mgp = c(2.3, 0.8, 0))
    overlay(d$X[[v]], w, 200, xlab = v, xlim = lims[[eg]][[v]]$xlim,
            ylim = lims[[eg]][[v]]$ylim, legend.pos = lims[[eg]][[v]]$pos)
    invisible(dev.off())
  }
  o$estimate
}

# coverage and width over n_sims replications; one row per replication, NA if it failed
run_eval <- function(eg, tauR, n_sims, n_boot, n_samples = 1000, alpha = 0.05,
                     cores = max(1L, parallel::detectCores() - 1L), tau = 1) {
  run <- function(itr) tryCatch({
    d <- draw(eg, n_samples, itr)
    o <- specrobust(d$Y, d$A, d$X, adj_sets, reg_mode = "lm", alpha = alpha,
                    n_boot = n_boot, n_cores = 1, seed = itr)
    lo <- o$candidate_ci_boot[, 1]; up <- o$candidate_ci_boot[, 2]; h <- range(c(lo, up))
    c(as.numeric(lo <= tau & tau <= up), as.numeric(h[1] <= tau & tau <= h[2]),
      as.numeric(o$ci[1] <= tauR & tauR <= o$ci[2]), up - lo, diff(h), diff(o$ci))
  }, error = function(e) rep(NA_real_, 8))
  res <- parallel::mclapply(seq_len(n_sims), run, mc.cores = cores)
  res[!vapply(res, function(z) is.numeric(z) && length(z) == 8L, TRUE)] <- list(rep(NA_real_, 8))
  res <- do.call(rbind, res)
  colnames(res) <- c("ireg1.ecov", "ireg2.ecov", "naive.ecov", "rewt.ecov",
                     "ireg1.len", "ireg2.len", "naive.len", "rewt.len")
  res
}

print_table <- function(res, eg, tauR) {
  rule <- strrep("=", 65)
  cat(sprintf("\n%s\nExample %s (%s), reweighted ATE: tau_R = %.10f\n%s\n",
              rule, eg, c("two confounders", "M-bias")[as.integer(eg)], tauR, rule))
  print(matrix(colMeans(res, na.rm = TRUE), 4, 2,
               dimnames = list(c("C.I. using adj set 1", "C.I. using adj set 2",
                                 "naive C.I. (range of prev C.I.'s)",
                                 "specification-robust C.I."),
                               c("emp coverage", "avg width"))), digits = 3)
  cat("\n")
}
