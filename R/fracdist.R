
##################################################
#
# Numerical Distribution Functions of
# Fractional Unit Root and Cointegration Tests
#
# Functions for Critical Values and P-values
#
# Lealand Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business
# University of Central Florida
#
# May 19, 2021
#
##################################################
#
# fracdist.R calculates critical values and p-values
#   for unit root and cointegration tests and
#   for rank tests in the FCVAR model.
#
# Dependencies:
#   None.
#
##################################################



#' Numerical CDFs for
#' Fractional Unit Root and Cointegration Tests
#'
#' A package for calculating numerical distribution functions of fractional
#' unit root and cointegration test statistics. The included functions calculate
#' critical values and P-values used in unit root tests, cointegration tests,
#' and rank tests in the Fractionally Cointegrated Vector
#' Autoregression (FCVAR) model (see Johansen and Nielsen, 2012).
#'
#' Simple tabulation is not a feasible approach for obtaining
#' critical values and P-values because these distributions depend
#' on a real-valued parameter \code{b} that must be estimated.
#' Instead, response surface regressions are used to obtain the numerical
#' distribution functions and combined by model averaging
#' across values taken from a series of tables.
#' As a function of the dimension of the problem, \code{q},
#' and a value of the fractional integration order \code{b},
#' this approach provides either a set of critical values or the asymptotic P-value
#' for any value of the likelihood ratio statistic.
#' The P-values and critical values are calculated by interpolating from the
#' quantiles on a grid of probabilities and values of the fractional integration order,
#' with separate tables for a range of values of cointegrating rank.
#'
#' The functions in this package are based on the functions and subroutines in the
#' Fortran program \code{fracdist.f} to accompany an article by MacKinnon and Nielsen (2014).
#' This program is available from the archive of the \emph{Journal of Applied Econometrics}
#' at \url{http://qed.econ.queensu.ca/jae/datasets/mackinnon004/}.
#' Alternatively, a C++ implementation of this program is also available; see
#' \url{https://github.com/jagerman/fracdist/blob/master/README.md} for details.
#'
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
#' "Likelihood inference for a fractionally cointegrated vector autoregressive model,"
#' \emph{Econometrica} 80, pp.2667-2732.
#'
#' @docType package
#' @name .fracdist
#' @return Returns \code{NULL}. Object included for description only.
NULL


#' Verify that fracdist parameters are valid parameter values
#'
#' This is a helper function for error handling in the fracdist package.
#' It is not intended to be used externally but might help diagnose
#' error messages from improperly chosen arguments.
#'
#' @param fracdist_params a list that may have the following elements:
#' \itemize{
#' \item \code{iq} An integer scalar rank parameter for the test, from 1 through 12.
#' This is often the difference in cointegration rank.
#' \item \code{iscon} An indicator that there is a constant intercept
#' term in the model.
#' \item \code{clevel} The numeric scalar level of significance.
#' \item \code{bb} The fractional integration parameter, which can take on values
#' between 0.0 and 2.0.
#' }
#' @return No return value when checks pass, otherwise execution halts
#' and an error message is printed.
#' @examples
#' # Test with iscon = 7 to see error message:
#' check_fracdist_params(list(iq = 2, iscon = 1))
#' # Test with iq = 13 to see error message:
#' check_fracdist_params(list(iq = 12, iscon = 1))
#' # Test with bb = -0.5 to see error message:
#' check_fracdist_params(list(iq = 12, iscon = 1, bb = 0.75))
#' # Test with bb = 2.5 to see error message:
#' check_fracdist_params(list(iq = 12, iscon = 1, bb = 1.5))
#' # Test with clevel = 2.5 to see error message:
#' check_fracdist_params(list(iq = 1, iscon = 1, clevel = 0.05))
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @export
#'
check_fracdist_params <- function(fracdist_params) {

  # Verify correct form of indicator for constant.
  if ("iscon" %in% names(fracdist_params)) {
    iscon <- fracdist_params$iscon
    if (!(class(iscon) %in% c('integer', 'numeric', 'logical')) |
        !(as.integer(iscon) %in% c(0, 1))) {
      stop("Indicator for constant (iscon) must be either 0 or 1.")
    }
  }

  # Verify that rank of system is within range.
  if ("iq" %in% names(fracdist_params)) {
    iq <- fracdist_params$iq
    if (!(class(iq) %in% c('integer', 'numeric')) |
        !(as.integer(iq %in% seq(1,12)))) {
      stop("Rank of system (iq) must be an integer from 1 through 12.")
    }
  }

  # Verify that the confidence level is in the unit interval.
  if ("clevel" %in% names(fracdist_params)) {
    clevel_list <- fracdist_params$clevel
    for (clevel in clevel_list) {
      if (!(class(clevel) %in% c('numeric')) |
          (clevel < 0 | clevel > 1)) {
        stop("Confidence level must be in the unit interval.")
      }
    }
  }

  # Check for valid fractional differencing order.
  if ("bb" %in% names(fracdist_params)) {
    bb <- fracdist_params$bb
    if (bb > 2.0) {
      stop(sprintf('Specified value of b = %3.2f is out of range 0.0 to 2.0.\n', bb),
           'If b > 2.0, you might consider differencing your series before\n',
           'estimating the model under consideration.')
    }
    if (bb < 0) {
      stop(sprintf('Specified value of b = %3.2f is out of range 0.0 to 2.0.\n', bb),
           'The package does not handle the case b < 0.00, ',
           'since bb is typically positive.\n')
    }

  }


}



#' Obtain Lookup Tables of Probabilities and Quantiles
#'
#' \code{get_fracdist_tab} selects a table of probabilities
#' and quantiles for a list of values of the fractional
#' integration parameter, corresponding to a particular rank
#' or dimension and the specification of a constant term.
#' @param iq An integer scalar rank parameter for the test, from 1 through 12.
#' This is often the difference in cointegration rank.
#' @param iscon An indicator that there is a constant intercept
#' term in the model.
#' @param dir_name A string name of directory in which the tables
#' are stored. This is not normally used, however, a user might want to
#' draw the tables from another location.
#' @param file_ext A string extension indicating the file format of
#' the tables. The default is \code{'rda'}, which is the format of
#' the tables included in the package. The tables may also be stored in
#' \code{'RData'} format. The tables accompanying the
#' original Fortran program, on which this package is based,
#' can be obtained in \code{'txt'} format. See the reference for details.
#' @return A data frame with three columns: \code{bbb} is the fractional
#' integration parameter, \code{probs} are the probability values,
#' and \code{xndf} is the \code{probs} quantile of the distribution.
#' @examples
#' frtab <- get_fracdist_tab(iq = 1, iscon = 0)
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @export
#'
get_fracdist_tab <- function(iq, iscon, dir_name = NULL,
                             file_ext = 'rda') {

  # Check for valid parameter values.
  check_fracdist_params(list(iq = iq, iscon = iscon))
  if (!(file_ext %in% c('rda', 'RData', 'txt'))) {
    stop("File extension (file_ext) must be one of 'rda', 'RData', or 'txt'.")
  }


  # Determine required file name.

  # Depends on whether or not there is a constant term.
  if (iscon == 0) {
    dfirst = 'frmapp'
  } else {
    dfirst = 'frcapp'
  }

  # Depends on the (difference in) cointerating ranks.
  dq <- sprintf('00%d', iq)
  dq <- substr(dq, nchar(dq) - 1, nchar(dq))


  if (file_ext == 'txt') {

    # Load original tables from Fortran version.

    in_file_name <- sprintf('%s/%s%s.txt', dir_name, dfirst, dq)

    # Initialize matrix and open connection.
    frtab <- data.frame(bbb = numeric(221*31),
                        probs = numeric(221*31),
                        xndf = numeric(221*31))

    frtab_lines <- readLines(con = in_file_name)
    lines_read <- 0

    # Loop on fractional integration parameter
    # to collect table from separate pieces.
    for (ib in 1:31) {

      # Read the value of b.
      bbb_str <- frtab_lines[lines_read + 1]
      lines_read <- lines_read + 1
      bbb <- as.numeric(substr(bbb_str, 6, 9))

      # Read the probabilities and quantiles.
      frtab_sub <- data.frame(bbb = rep(bbb, 221),
                              probs = numeric(221),
                              xndf = numeric(221))

      frtab_lines_sub <- frtab_lines[seq(lines_read + 1, lines_read + 221)]
      lines_read <- lines_read + 221
      frtab_sub[, 'probs'] <- as.numeric(substr(frtab_lines_sub, 1, 6))
      frtab_sub[, 'xndf'] <- as.numeric(substr(frtab_lines_sub, 9, 25))

      # Append to the full table.
      frtab[seq((ib - 1)*221 + 1, ib*221), ] <- frtab_sub

    }

  } else if (file_ext == 'RData') {

    # Load compressed files for R.
    in_file_name <- sprintf('%s/%s%s.RData', dir_name, dfirst, dq)
    frtab <- get(load(file = in_file_name))

  } else if (file_ext == 'rda') {

    if (is.null(dir_name)) {

      # Read the table from within the data folder of the R package.
      in_file_name <- sprintf('%s%s', dfirst, dq)
      frtab <- get(in_file_name) # Not necessary, since already assigned.

    } else {

      # Load compressed files for R.
      in_file_name <- sprintf('%s/%s%s.rda', dir_name, dfirst, dq)
      frtab <- get(load(file = in_file_name))

    }

  } else {
    stop('File extension not supported.')
  }


  return(frtab)
}


#' Interpolate Critical Values Local to Fractional Integration Parameter
#'
#' \code{blocal} calculates an approximation to the CDF for a particular
#' value of the fractional integration parameter b.
#' @param nb The integer number of values of the fractional integration
#' parameter \eqn{b}. The default is 31, which is the number of gridpoints
#' between 0.5 and 2.0 in the lookup tables.
#' @param bb The fractional integration parameter, which can take on values
#' between 0.5 and 2.0.
#' @param estcrit A numeric scalar that is an estimate of the critical value,
#' a quantile, corresponding to a particular probability.
#' @param bval A numeric vector of the fractional integration parameter \eqn{b}
#' ranging between 0.5 and 2.0 to match those in the lookup tables.
#' @return A numeric scalar that is the sum of regression coefficients
#' from a response surface regression.
#' @examples
#' frtab <- get_fracdist_tab(iq = 3, iscon = 0)
#' bval <- unique(frtab[, 'bbb'])
#' probs <- unique(frtab[, 'probs'])
#' estcrit <- frtab[frtab[, 'probs'] == probs[201], 'xndf']
#' bedf_i <- blocal(nb = 31, bb = 0.75, estcrit, bval)
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @note The fractional integration parameter \code{b} must lie within
#' the interval \eqn{[0.50, 2]}. If \code{b} is less than 0.5,
#' the relevant distribution is chi-squared. A value above 2 is less common
#' but might be avoided by differencing the data before estimating the model
#' under consideration.
#' @export
#'
blocal <- function(nb = 31, bb, estcrit, bval) {

  if (bb < 0.51 | bb > 2.0) {
    stop(sprintf('Specified value of b = %3.2f is out of range 0.51 to 2.0.\n', bb),
         'If b < 0.50, the chi-squared distribution might be appropriate.\n',
         'If b > 2.0, you might consider differencing your series before\n',
         'estimating the model under consideration.')
  }

  # Weights on neighboring quantiles are based on trapezoidal kernel.
  weight <- 1.0 - 5.0*abs(bval - bb)
  weight[weight < 0] <- 0
  # Includes at most 9 observations.

  # Determine the quantiles to be used for interpolation.
  # jbot will be index of lowest b that gets positive weight.
  jbot <- which(weight > 0)[1]
  nobs <- sum(weight > 0)
  weight_sel <- seq(jbot, jbot + nobs - 1)

  # Create a dataset for interpolation by regression.
  yx_mat <- data.frame(y = numeric(9),
                       x1 = numeric(9),
                       x2 = numeric(9),
                       x3 = numeric(9))

  yx_mat[1:nobs, 'y'] <- weight[weight_sel]*estcrit[weight_sel]
  yx_mat[1:nobs, 'x1'] <- weight[weight_sel]
  yx_mat[1:nobs, 'x2'] <- weight[weight_sel]*bval[weight_sel]
  yx_mat[1:nobs, 'x3'] <- weight[weight_sel]*bval[weight_sel]^2

  lm_olsqc <- stats::lm(formula = y ~ x1 + x2 + x3 - 1, data = yx_mat)

  bfit <- sum(lm_olsqc$coefficients*bb^seq(0,2))

  return(bfit)
}



#' Calculate P-values for Fractional Unit Root and Cointegration Tests
#'
#' \code{fpval} calculates P-values for a particular value of the observed
#' statistic and a set of intermediate calculations
#' that are output from other functions in the \code{fracdist} package.
#' @param npts An integer number of points for local approximation
#' of the EDF near the observed value of \code{stat}.
#' It is usually 9, unless near a boundary.
#' @param iq An integer scalar rank parameter for the test, from 1 through 12.
#' This is often the difference in cointegration rank.
#' @param stat A numeric scalar value of the test statistic.
#' @param probs A numeric vector of probabilities over which an approximating
#' empirical distribution function is obtained, taken from precalculated tables.
#' @param bedf A numeric vector of quantiles of numerical distribution for specified
#' value of fractional integration order \eqn{b} or values of \eqn{b} and \eqn{d},
#' depending on the particular model. Each element is the output of the function \code{blocal}.
#' @param ginv A numeric vector of quantiles of the approximating chi-squared distribution.
#' @return A numeric scalar P-value.
#' @examples
#' frtab <- get_fracdist_tab(iq = 3, iscon = 0)
#' bval <- unique(frtab[, 'bbb'])
#' probs <- unique(frtab[, 'probs'])
#' bedf <- rep(NA, length(probs))
#' for (i in 1:length(probs)) {
#'     estcrit <- frtab[frtab[, 'probs'] == probs[i], 'xndf']
#'     bedf[i] <- blocal(nb = 31, bb = 0.75, estcrit, bval)
#' }
#' fpval(npts = 9, iq = 3, stat = 3.84, probs, bedf, ginv = qchisq(probs, df = 3^2))
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @seealso \code{fracdist_pvalues} for the calculation of P-values including any
#' intermediate calculations.
#' @export
#'
fpval <- function(npts = 9, iq, stat, probs, bedf, ginv) {

  # Check for valid parameter values.
  check_fracdist_params(list(iq = iq))

  # nomax <- 25
  # nvmax <- 3
  ndf <- iq**2

  # Deal with extreme cases.
  btiny <- 0.5*bedf[1]
  bbig <- 2.0*bedf[221]
  if (stat < btiny) {
    pval <- 1.0
    return(pval)
  }
  if (stat > bbig) {
    pval <- 0.0
    return(pval)
  }


  # Find critical value closest to test statistic.
  diff <- abs(stat - bedf)
  imin <- which.min(diff)[1]
  diffm <- diff[imin]

  # nph <- npts/2 # Division by integer does floor.
  nph <- floor(npts/2)
  nptop <- 221 - nph


  # Create a dataset for interpolation by regression.
  yx_mat <- data.frame(y = numeric(npts),
                       x1 = numeric(npts),
                       x2 = numeric(npts),
                       x3 = numeric(npts))
  # The form depends on the proximity to the endpoints.

  if (imin > nph & imin < nptop) {
    # imin is not too close to the end. Use np points around stat.


    # Populate the dataset for interpolation by regression.
    ic <- imin - nph - 1 + seq(1, npts)

    yx_mat[, 'y'] <- ginv[ic]
    yx_mat[, 'x1'] <- 1.0
    yx_mat[, 'x2'] <- bedf[ic]
    yx_mat[, 'x3'] <- bedf[ic]^2

  } else {
    # imin is close to one of the ends. Use points from imin +/- nph to end.

    if (imin < nph) {

      np1 <- imin + nph
      np1 <- max(np1, 5)

      # Populate the dataset for interpolation by regression.
      ic <- seq(1, np1)
      yx_mat[1:np1, 'y'] <- ginv[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- bedf[ic]
      yx_mat[1:np1, 'x3'] <- bedf[ic]^2

    } else {

      np1 <- 222 - imin + nph
      np1 <- min(max(np1, 5), 221)

      # Populate the dataset for interpolation by regression.
      ic <- 222 - seq(1, np1)

      yx_mat[1:np1, 'y'] <- ginv[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- bedf[ic]
      yx_mat[1:np1, 'x3'] <- bedf[ic]^2

    }

  }

  # Run regression and obtain p-value.
  lm_olsqc <- stats::lm(formula = y ~ x1 + x2 + x3 - 1, data = yx_mat)

  crfit <- sum(lm_olsqc$coefficients*stat^seq(0,2))
  crfit <- max(crfit, 10^(-6))

  pval <- stats::pchisq(crfit, df = ndf)
  pval <- 1.0 - pval


  # Accuracy only claimed up to the fourth-fifth decimal place.
  pval <- round(pval, 4)

  return(pval)
}


#' Calculate P-values for Fractional Unit Root and Cointegration Tests
#'
#' \code{fracdist_pvalues} calculates P-values for a particular value of the observed
#' statistic.
#' @param iq An integer scalar rank parameter for the test, from 1 through 12.
#' This is often the difference in cointegration rank.
#' @param iscon An indicator that there is a constant intercept
#' term in the model.
#' @param dir_name A string name of directory in which the approximating tables
#' are stored. This is not normally used, since sufficient tables are included in the package.
#' However, a user might want to draw the tables from another location.
#' @param bb The fractional integration parameter, which can take on values
#' between 0.5 and 2.0.
#' @param stat A numeric scalar value of the test statistic.
#' @return A numeric scalar P-value.
#' @examples
#' fracdist_pvalues(iq = 1, iscon = 0, bb = 0.73, stat = 3.84)
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @seealso Calls \code{fpval} for the calculation of P-values after
#' performing some intermediate calculations.
#' @export
#'
fracdist_pvalues <- function(iq, iscon, dir_name = NULL, bb, stat) {

  # Check for valid parameter values.
  check_fracdist_params(list(iq = iq, iscon = iscon, bb = bb))

  # Obtain relevant table of statistics.
  frtab <- get_fracdist_tab(iq, iscon, dir_name)
  bval <- unique(frtab[, 'bbb'])
  # bval <- bval[order(bval)]
  probs <- unique(frtab[, 'probs'])
  # probs <- ginv_tab[, 'probs']
  # probs <- probs[order(probs)]

  # Calculate the approximate CDF for this particular value of b.
  nb <- 31
  np <- 221
  bedf <- rep(NA, np)
  for (i in 1:np) {

    prob_i <- probs[i]
    estcrit <- frtab[frtab[, 'probs'] == prob_i, 'xndf']

    bedf[i] <- blocal(nb, bb, estcrit, bval)

  }

  # Obtain inverse CDF of the Chi-squared distribution.
  # ginv <- stats::qchisq(probs, df = iq^2)
  ginv <- ginv_tab[, sprintf('iq_%d', iq)]
  # In the interpolation regressions, this is the regressand
  # for fpval and the regressors for fpcrit.

  # Calculate p-values.
  pval <- fpval(npts = 9, iq, stat, probs, bedf, ginv)

  return(pval)
}


#' Calculate Critical Values for Fractional Unit Root and Cointegration Tests
#'
#' \code{fcrit} calculates critical values for a particular level of significance
#' and a set of intermediate calculations
#' that are output from other functions in the \code{fracdist} package.
#' @param npts An integer number of points for local approximation
#' of the EDF near the level of significance \code{clevel}.
#' It is usually 9, unless near a boundary.
#' @param iq An integer scalar rank parameter for the test, from 1 through 12.
#' This is often the difference in cointegration rank.
#' @param clevel The numeric scalar level of significance.
#' @param probs A numeric vector of probabilities over which an approximating
#' empirical distribution function is obtained, taken from precalculated tables.
#' @param bedf A numeric vector of quantiles of numerical distribution for specified
#' value of fractional integration order \eqn{b} or values of \eqn{b} and \eqn{d},
#' depending on the particular model. Each element is the output of the function \code{blocal}.
#' @param ginv A numeric vector of quantiles of the approximating chi-squared distribution.
#' @return A numeric scalar critical value, a quantile of the distribution.
#' @examples
#' frtab <- get_fracdist_tab(iq = 3, iscon = 0)
#' bval <- unique(frtab[, 'bbb'])
#' probs <- unique(frtab[, 'probs'])
#' bedf <- rep(NA, length(probs))
#' for (i in 1:length(probs)) {
#'     estcrit <- frtab[frtab[, 'probs'] == probs[i], 'xndf']
#'     bedf[i] <- blocal(nb = 31, bb = 0.75, estcrit, bval)
#' }
#' fpcrit(npts = 9, iq = 3, clevel = 0.05, probs, bedf, ginv = qchisq(probs, df = 3^2))
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @seealso fracdist_values for the calculation of critical values and
#' P-values including any intermediate calculations.
#' @export
#'
fpcrit <- function(npts = 9, iq, clevel, probs, bedf, ginv) {

  # Check for valid parameter values.
  check_fracdist_params(list(iq = iq, clevel = clevel))

  # nomax <- 25
  # nvmax <- 3
  ndf <- iq**2

  # Handle extreme cases.
  ptiny <- 0.0001
  pbig <- 0.9999
  if (clevel < ptiny) {
    ccrit <- bedf[221]
    return(ccrit)
  }
  if (clevel > pbig) {
    ccrit <- bedf[221]
    return(ccrit)
  }


  # Find probability closest to test level.
  cquant <- 1.0 - clevel
  diff <- abs(cquant - probs)
  imin <- which.min(diff)[1]
  diffm <- diff[imin]

  gcq <- stats::qchisq(cquant, df = ndf)

  # nph <- npts/2
  nph <- floor(npts/2)
  nptop <- 221 - nph


  # Create a dataset for interpolation by regression.
  yx_mat <- data.frame(y = numeric(npts),
                       x1 = numeric(npts),
                       x2 = numeric(npts),
                       x3 = numeric(npts))
  # The form depends on the proximity to the endpoints.

  if (imin > nph & imin < nptop) {
    # imin is not too close to the end. Use np points around stat.

    # Populate the dataset for interpolation by regression.
    ic <- imin - nph - 1 + seq(1, npts)
    yx_mat[, 'y'] <- bedf[ic]
    yx_mat[, 'x1'] <- 1.0
    yx_mat[, 'x2'] <- ginv[ic]
    yx_mat[, 'x3'] <- ginv[ic]^2

  } else {
    # imin is close to one of the ends. Use points from imin +/- nph to end.

    if (imin < nph) {

      np1 <- imin + nph
      np1 <- max(np1, 5)

      # Populate the dataset for interpolation by regression.
      ic <- seq(1, np1)
      yx_mat[1:np1, 'y'] <- bedf[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- ginv[ic]
      yx_mat[1:np1, 'x3'] <- ginv[ic]^2

    } else {

      np1 <- 222 - imin + nph
      np1 <- max(np1, 5)

      # Populate the dataset for interpolation by regression.
      ic <- 222 - seq(1, np1)
      yx_mat[1:np1, 'y'] <- bedf[ic]
      yx_mat[1:np1, 'x1'] <- 1.0
      yx_mat[1:np1, 'x2'] <- ginv[ic]
      yx_mat[1:np1, 'x3'] <- ginv[ic]^2

    }

  }

  # Run regression and obtain critical value.
  lm_olsqc <- stats::lm(formula = y ~ x1 + x2 + x3 - 1, data = yx_mat)

  ccrit <- sum(lm_olsqc$coefficients*gcq^seq(0,2))

  # Accuracy only claimed up to the fourth-fifth decimal place.
  ccrit <- round(ccrit, 4)

  return(ccrit)
}


#' Calculate Critical Values or P-values for Fractional Unit Root and Cointegration Tests
#'
#' \code{fracdist_values} calculates either critical Values or P-values for
#' for fractional unit root and cointegration tests
#' @param iq An integer scalar rank parameter for the test, from 1 through 12.
#' This is often the difference in cointegration rank.
#' @param iscon An indicator that there is a constant intercept
#' term in the model.
#' @param dir_name A string name of directory in which the approximating tables
#' are stored. This is not normally used, since sufficient tables are included in the package.
#' However, a user might want to draw the tables from another location.
#' @param bb The fractional integration parameter, which can take on values
#' between 0.0 and 2.0.
#' @param stat A numeric scalar value of the test statistic. This is only
#' used if P-values are required.
#' @param ipc A logical indicator to calculate a P-value.
#' If \code{ipc == FALSE} critical values are calculated instead.
#' @param clevel The numeric scalar level of significance. The default is to
#' calculate critical values for the conventional levels of significance:
#' \code{clevel = c(0.01, 0.05, 0.10)}.
#' @return Either a numeric scalar P-value, if \code{ipc == TRUE},
#' otherwise, a numeric vector of critical values, the same length as \code{clevel}.
#' @examples
#' # Calculate P-values:
#' fracdist_values(iq = 1, iscon = 0, bb = 0.43, stat = 3.84)
#' fracdist_values(iq = 1, iscon = 0, bb = 0.73, stat = 3.84)
#' # Calculate critical values:
#' fracdist_values(iq = 1, iscon = 0, bb = 0.73, ipc = FALSE, clevel = 0.05)
#' fracdist_values(iq = 1, iscon = 0, bb = 0.43, ipc = FALSE, clevel = 0.05)
#' fracdist_values(iq = 1, iscon = 0, bb = 0.73, ipc = FALSE)
#' @references James G. MacKinnon and Morten \enc{Ø}{O}rregaard Nielsen,
#' "Numerical Distribution Functions of Fractional Unit Root and Cointegration Tests,"
#' \emph{Journal of Applied Econometrics}, Vol. 29, No. 1, 2014, pp.161-171.
#' @references Johansen, S. and M. \enc{Ø}{O}. Nielsen (2012).
#' "Likelihood inference for a fractionally cointegrated vector autoregressive model,"
#' \emph{Econometrica} 80, pp.2667-2732.
#' @seealso Calls \code{fpval} to calculate P-values
#' or \code{fpcrit} to calculate critical values,
#' after performing some intermediate calculations.
#' @note For fractional integration orders between 0 and 0.5, the chi-square
#' distribution is used. See Johansen and Nielsen 2012 for details.
#' @export
#'
fracdist_values <- function(iq, iscon, dir_name = NULL, bb, stat,
                            ipc = TRUE, clevel = c(0.01, 0.05, 0.10)) {

  # Check for valid parameter values.
  check_fracdist_params(list(iq = iq, iscon = iscon, bb = bb, clevel = clevel))

  # Use chi-square distribution for fractional integration orders between 0 and 0.5.
  # Strictly speaking, the chi-squared distribution only applies for b <= 0.5.
  # However, the response surface might be unreliable near the boundary and
  # the chi-squared distribution is a more reliable approximation.
  if (bb > 0 & bb < 0.51) {

    if (ipc == TRUE) {

      # Calculate p-values.
      outval <- 1 - stats::pchisq(q = stat, df = iq^2)

    } else {

      # Calculate critical values.
      outval <- stats::qchisq(p = 1 - clevel, df = iq^2)

    }

  } else {

    # Obtain relevant table of statistics.
    frtab <- get_fracdist_tab(iq, iscon, dir_name)
    bval <- unique(frtab[, 'bbb'])
    # bval <- bval[order(bval)]
    probs <- unique(frtab[, 'probs'])
    # probs <- ginv_tab[, 'probs']
    # probs <- probs[order(probs)]

    # Calculate the approximate CDF for this particular value of b.
    nb <- 31
    np <- 221
    bedf <- rep(NA, np)
    for (i in 1:np) {

      prob_i <- probs[i]
      estcrit <- frtab[frtab[, 'probs'] == prob_i, 'xndf']

      bedf[i] <- blocal(nb, bb, estcrit, bval)

    }

    # Obtain inverse CDF of the Chi-squared distribution.
    # ginv <- stats::qchisq(probs, df = iq^2)
    ginv <- ginv_tab[, sprintf('iq_%d', iq)]
    # In the interpolation regressions, this is the regressand
    # for fpval and the regressors for fpcrit.

    # Perform required calculation.
    if (ipc == TRUE) {

      # Calculate p-values.
      outval <- fpval(npts = 9, iq, stat, probs, bedf, ginv)

    } else {

      # Calculate critical values.

      # Might specify a list of critical values.
      outval <- rep(NA, length(clevel))
      for (i in 1:length(clevel)) {
        clevel_i <- clevel[i]
        outval[i] <- fpcrit(npts = 9, iq, clevel = clevel_i, probs, bedf, ginv)

      }

    }

  }

  return(outval)
}



##################################################
# End
##################################################
