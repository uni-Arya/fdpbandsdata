# Set working directory to fdpbandsdata/

# # Command-line arguments
# arguments   <- commandArgs(trailingOnly = TRUE)
# c_load      <- as.numeric(arguments[1])
# lambda_load <- as.numeric(arguments[2])
# n_decoys    <- as.numeric(arguments[3])

c_load      <- 0.5
lambda_load <- 0.5
n_decoys    <- 1

# Monte-Carlo computation of u_gamma(\Delta) for TDC-UB
compute_u_gammas <- function(c_load, lambda_load, n_decoys) {
  filename <- sprintf(
                "data-raw/u_gamma-tables/utable_c%.5f_lam%.5f",
                c_load,
                lambda_load
              )
  if(!file.exists(filename)) {
    # Adjust c and lambda in case of infinite decimal representations
    d1      <- n_decoys + 1
    c       <- floor(c_load * d1 + 1e-13) / d1
    lambda  <- floor(lambda_load * d1 + 1e-13) / d1

    ugamma <- function(n_mc, D, c, lambda, confs) {
      cat(details, "\n") # See below
      cat(
        "Beginning computation --",
        format(Sys.time(), "%a %b %d %X %Y"), "\n"
      )

      # Probability of a decoy-winning hypothesis, given it was counted
      R_c_lam <- (1 - lambda) / (c + 1 - lambda)

      # Maximum index and length of confs
      nds <- length(D)
      ncs <- length(confs)

      # For numerical precision
      eps   <- .Machine$double.eps * 1e4
      index <- n_mc * confs
      hi    <- ceiling(index - eps)

      # Storage of important quantities
      result_list <- as.list(rep(0, 3))
      for (i in 1:3) {
        result_list[[i]] <- matrix(0, nrow = nds, ncol = length(confs))
      }

      if (nds > 1) {
        temp  <- rnbinom(n_mc, 1, R_c_lam)
        gtemp <- 1 - pnbinom(temp - 1, 1, R_c_lam)

        # Find sigma and rho
        stemp <- sort(gtemp, partial = unique(hi))
        sigma <- stemp[hi]
        rho   <- sapply(
                   1:ncs,
                   function(i)
                     max(stemp[1:hi[i]][stemp[1:hi[i]] < sigma[i] - eps])
                 )

        # Coinflip probability
        upper_prop <- sapply(
                        1:ncs,
                        function(i) mean(stemp <= sigma[i])
                      )
        lower_prop <- sapply(
                        1:ncs,
                        function(i) mean(stemp <= rho[i])
                      )
        lower_p    <- (upper_prop - confs) / (upper_prop - lower_prop)

        result_list[[1]][1,] <- rho
        result_list[[2]][1,] <- sigma
        result_list[[3]][1,] <- lower_p

        for (d in 2:nds) {
          temp  <- temp + rnbinom(n_mc, 1, R_c_lam)
          gtemp <- pmin(
                     gtemp,
                     1 - pnbinom(temp - 1, d, R_c_lam)
                   )

          # Find sigma and rho
          stemp <- sort(gtemp, partial = unique(hi))
          sigma <- stemp[hi]
          rho   <- sapply(
                     1:ncs,
                     function(i)
                       max(stemp[1:hi[i]][stemp[1:hi[i]] < sigma[i] - eps])
                   )

          # Coinflip probability
          upper_prop <- sapply(
                          1:ncs,
                          function(i) mean(stemp <= sigma[i])
                        )
          lower_prop <- sapply(
                          1:ncs,
                          function(i) mean(stemp <= rho[i])
                        )
          lower_p    <- (upper_prop - confs) / (upper_prop - lower_prop)

          result_list[[1]][d,] <- rho
          result_list[[2]][d,] <- sigma
          result_list[[3]][d,] <- lower_p

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

      return(result_list)
    }

    # Number of hypotheses
    m <- 5e4

    # Desired confidence levels
    confs <- c(0.01, 0.025, 0.05, 0.10, 0.20, 0.50)

    # For d in D, compute u_gamma(\Delta_d) where \Delta_d = {1, ..., d}
    D <- 1:m

    # Number of Monte-Carlo runs
    n_mc <- 2e6

    details <- sprintf(
      "(c: %.5f, lam: %.5f, n_decoys: %d, n_mc: %s)",
      c_load, lambda_load, n_decoys, format(n_mc, scientific = T)
    )

    # Generate u_gamma tables
    utable <- ugamma(n_mc, D, c, lambda, confs)
    for (i in 1:3) {
      colnames(utable[[i]]) <- 1 - confs
    }

    # Save utable
    save(utable, file = filename)
  } else {
    stop("u_gamma table for the given (c, lambda) pair already exists.")
  }
}

# Compute parameter u_gamma for FDP-UBm
compute_u_gammas(c_load, lambda_load, n_decoys)
