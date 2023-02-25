# Set working directory to fdpbandsdata/

# # Command-line arguments
# arguments   <- commandArgs(trailingOnly = TRUE)
# c_load      <- as.numeric(arguments[1])
# lambda_load <- as.numeric(arguments[2])
# n_decoys    <- as.numeric(arguments[3])

c_load      <- 0.5
lambda_load <- 0.5
n_decoys    <- 1

# This function calculates the 1 - gamma quantile z(D, c, lambda, gamma)
# of Z(D) := max_{d in D} \hat{U}_d for D = [k], with k varying from 1
# to a pre-specified maximum, m (typically 10,000 or more).

# n_mc: number of Monte-Carlo runs
# D: candidate rejection thresholds (see above)
# c, lambda: competition parameters
# confs: set of desired conf levels
quantile_table = function(n_mc, D, c, lambda, confs) {
  # Probability of a decoy-winning hypothesis, given it was counted
  R_c_lam <- (1 - lambda) / (c + 1 - lambda)

  # Maximum index of consideration
  nds <- length(D)

  cat(details, "\n")
  cat("Beginning computation --", format(Sys.time(), "%a %b %d %X %Y"), "\n")

  # Storage of quantiles
  quantile_table <- matrix(0, nrow = nds, ncol = length(confs))

  if (nds > 1) {
    temp                <- rnbinom(n_mc, 1, R_c_lam)
    htemp               <- (R_c_lam * temp - (1 - R_c_lam) * 1) /
      sqrt((1 - R_c_lam) * 1)
    quantile_table[1, ] <- as.vector(quantile(htemp, 1 - confs))

    for (d in 2:nds) {
      temp                <- temp + rnbinom(n_mc, 1, R_c_lam)
      htemp               <- pmax(
        htemp,
        (R_c_lam * temp - (1 - R_c_lam) * d) /
          sqrt((1 - R_c_lam) * d)
      )
      quantile_table[d, ] <- as.vector(quantile(htemp, 1 - confs))

      if (d %% floor(nds/10) == 0) {
        cat(
          "Current d:", d, "--", "c:", c, "lam:", lambda,
          format(Sys.time(), "%a %b %d %X %Y"),
          "\n"
        )
      }
    }
  } else {
    stop("Argument \"D\" is invalid.")
  }

  # Each column of the output represents a fixed confidence gamma, and each
  # row represents the (1 - gamma)-quantile of Z(l) for a fixed l in
  # {1, 2, ..., k}.
  colnames(quantile_table) <- as.character(1 - confs)
  return(quantile_table)
}

# If it doesn't exist, compute the quantile table for a fixed (c, lambda) pair
filename <- sprintf("data-raw/quantile-tables/qtable_c%.5f_lam%.5f", c_load, lambda_load)
if (!file.exists(filename)) {
  # Choice of confidence levels for quantile
  confs <- c(0.01, 0.025, 0.05, 0.10, 0.20, 0.50)

  # Number of hypotheses
  m <- 5e4

  # Number of Monte-Carlo runs
  n_mc <- 2e6

  # Details to be displayed during computation
  details <- sprintf(
    "(c: %.5f, lam: %.5f, n_decoys: %d, n_mc: %s)",
    c_load, lambda_load, n_decoys, format(n_mc, scientific = T)
  )

  # Adjust c and lambda in case of infinite decimal representations
  d1      <- n_decoys + 1
  fixed_c <- floor(c_load * d1 + 1e-13) / d1
  lambda  <- floor(lambda_load * d1 + 1e-13) / d1

  # D = [k] for k varying from 1 to m
  D <- 1:m

  # Generate quantiles
  qtable <- quantile_table(n_mc, D, fixed_c, lambda, confs)

  # Save quantiles
  save(qtable, file = filename)

  cat(
   "Quantile generation for FDP-SBm complete --",
   format(Sys.time(), "%a %b %d %X %Y")
  )

} else {
  stop("Quantile table for the given (c, lambda) pair already exists.")
}
