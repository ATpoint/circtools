#' Detection of rhythmicity using limma with a cosinor model.
#'
#' The function detects periodic rhythmicity based on the cosinor model as described in the limma design matrix guide:
#' \url{https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html#cyclical-models}
#' For this it uses the limma framework to first fit per-observation cosinor models with \code{limma::lmFit()} and then
#' obtains moderated pvalues and FDRs with \code{limma::eBayes()}. It also returns estimates for mesor, amplitude and acrophase
#' based on the model coefficients. The function assumes that input data are "ready-to-use",  meaning that normalization, prefiltering (if applicable) appropriate 
#' for this type of data has been done. 
#'
#' Since the input data `x` (see argument description) can be both a numeric matrix or an `EList`, the user has a great flexibility over the analysis:
#' For an ordinary limma analysis (default settings) `x` could be a numeric matrix with data normalized and on log2 scale.
#' For a limma-trend analysis (\code{eBayes_args=list(trend=TRUE)}) `x` could be a numeric matrix of logCPMs from an RNA-seq experiment. 
#' For a limma-voom analysis `x` could be the output of \code{limma::voom()}.
#' Since limma automatically uses **weights** if provided in the EList as x$weights, one can put them there in case weighted regression is desired.
#' Weights could be estimated using \code{limma::arrayWeights} based on the cosinor model. A `model.matrix()` for the cosinor model can be obtained with \code{make_design()}.
#' Alternatively, `x` could be the output of \code{limma::voomWithQualityWeights()}. 
#' 
#'
#' @param x A numeric matrix or data.frame with observations in rows and samples in columns. Alternatively, an \code{EList} object.
#' @param time A numeric vector with the time information for each column of x.
#' @param period The numeric period of the rhythmicity, default 24.
#' @param lmFit_args A list with arguments for \code{limma::lmFit()}.
#' @param eBayes_args A list with arguments for \code{limma::eBayes()}.
#'
#' @examples
#' period <- 24
#' n_reps <- 1
#' time <- rep(seq(0, 24, by = 2), each = n_reps)
#' acrophases <- c(1, 10, 21)
#' y <- lapply(acrophases, function(ac) make_data(time = time, acrophase = ac, error_sd = .1)$value)
#'
#' # x is a numeric matrix
#' x <- matrix(unlist(y), byrow = TRUE, ncol = length(time))
#' rownames(x) <- paste0("acrophase_", acrophases)
#' run_cosinor(x = x, time = time)
#'
#' # x is an EList
#' xe <- new("EList", list(E = x))
#' run_cosinor(x = xe, time = time)
#'
#' # x is an EList with weights
#' set.seed(1)
#' xew <- new("EList", list(E = x, weights = rnorm(ncol(x), 1)))
#' run_cosinor(x = xew, time = time)
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom methods is new
#'
#' @export
#'
run_cosinor <- function(x, time, period = 24, lmFit_args = list(), eBayes_args = list()) {
  is_class <- is(x, "EList") | is(x, "data.frame") | is(x, "matrix")
  if (!is_class) stop("x must either be a numeric matrix or data.frame, or a limma EList")

  # Ensure consistent input format and row/colnames are set
  if (is(x, "EList")) {
    x$E <- as.matrix(x$E)

    rownames(x$E) <- if (is.null(rownames(x$E))) paste0("row_", 1:nrow(x$E)) else rownames(x$E)
    colnames(x$E) <- if (is.null(colnames(x$E))) paste0("sample_", 1:ncol(x$E)) else colnames(x$E)

    y <- x$E
  } else {
    x <- as.matrix(x)

    rownames(x) <- if (is.null(rownames(x))) paste0("row_", 1:nrow(x)) else rownames(x)
    colnames(x) <- if (is.null(colnames(x))) paste0("sample_", 1:ncol(x)) else colnames(x)

    y <- x
  }

  # Check input is numeric and finite
  is_numeric <- is.numeric(as.vector(y))
  has_NAs <- sum(!is.finite(as.vector(y))) > 0
  if (!is_numeric | has_NAs) stop("x contains non-numeric entries and/or NAs/NaNs/infinite values")
  rm(y)

  # Check that time is numeric and of same length as samples in x
  if (!is.numeric(time)) stop("time must be numeric")
  len_time <- length(time)
  len_x <- ncol(x)
  if (len_time != len_x) stop("time must have the same length as ncol(x) -- one time entry per column")

  if (!is.numeric(period) | period < 0) stop("period must be numeric and positive")
  if (!is.list(lmFit_args)) stop("lmFit_args must be a list")
  if (!is.list(eBayes_args)) stop("eBayes_args must be a list")

  # Fit cosinor model and extract stats
  design <- make_design(time = time, period = period)
  fit <- do.call(lmFit, c(list(x, design = design), lmFit_args))
  fit <- do.call(eBayes, c(list(fit), eBayes_args))
  coef_pos <- grep("^sinphase$|^cosphase$", colnames(fit$coefficients))
  tt <- topTable(fit = fit, coef = coef_pos, number = Inf, sort.by = "none")
  tt <- data.frame(id = rownames(x), tt)

  tt$amplitude <- sqrt(rowSums(tt[, c("sinphase", "cosphase")]^2))
  tt$acrophase <- (atan2(tt$sinphase, tt$cosphase) / 2 / pi * period + period) %% period
  tt <- tt[, c("id", "P.Value", "adj.P.Val", "amplitude", "acrophase", "AveExpr", "sinphase", "cosphase", "F")]

  return(tt)
}

#' Simulate cosinor data
#'
#' Accept numeric timepoints and cosinor parameters and returns a data.frame with time and values.
#'
#' @param time A numeric vector with time information
#' @param amplitude The amplitude of the cosinor curve
#' @param acrophase The acrophase of the cosinor curve
#' @param period The period of the time series
#' @param mesor The mesor of the cosinor curve
#' @param error_sd An error value to simulate noise
#' @param use_seed A random seed to make error_sd reproducible
#'
#' @examples
#' tt <- seq(0, 23.9, by = .1)
#' data <- make_data(time = tt, acrophase = 5)
#' plot(data$time, data$value, pch = 20)
#' abline(v = 5)
#'
#' @importFrom stats rnorm
#'
#' @export
#'
make_data <- function(time, amplitude = 1, acrophase = 1, period = 24, mesor = 0, error_sd = .1, use_seed = 1) {
  if (!is.numeric(time) | sum(time < 0) > 0 | length(time) == 1) stop("time must be a numeric vector with non-negative positive values")

  set.seed(use_seed)
  error <- rnorm(length(time), mean = 0, sd = error_sd)
  v <- mesor + amplitude * cos(2 * pi / period * (time - acrophase)) + error

  data <- data.frame(time = time, value = v)

  return(data)
}

#' Make a design matrix using a cosinor model.
#'
#' @param time A numeric vector with timepoints
#' @param period The period of the time series
#'
#' @importFrom stats model.matrix
#' @export
#'
make_design <- function(time, period = 24) {
  if (!is.numeric(time)) stop("time must be numeric")
  if (!is.numeric(period) | period < 0) stop("period must be numeric and positive")

  d <- data.frame(
    sinphase = sin(2 * pi * time / period),
    cosphase = cos(2 * pi * time / period)
  )

  design <- model.matrix(~ sinphase + cosphase, d)

  return(design)
}

#' Calculate densities of periodic values such as acrophases between user-defined start and end points, respecting the circular nature of the data.
#' 
#' Density calculation is prone to skewed density estimates at the start and end of the data range if the method is not aware of the circular nature of the data.
#' For that this function appends the data to itself at the front and rear end, and then uses \code{stats::density} with a default bandwidth of 0.5 to calculate densities
#' that can be used for histogram/ridge-like plots. The user should provide a start point `from` representing the smallest possible value of the periodoc values, and an end point
#' `to` which represents the period of the circular data / the time series.
#'
#' @param x Numeric vector with input value.
#' @param from,to The left and right-most points of the grid at which the density is to be estimated.
#'   The defaults are 0 and 24, assuming a 24h clock cycle. The "to" argument is typically the period of the time series.
#' @param density_args arguments for \code{stats::density}
#'
#' @examples
#'
#' set.seed(1)
#' data <- abs(c(rnorm(1000, 12, 5), rnorm(500, 23))) %% 24
#'
#' # linear density estimation incorrectly has low density at around x=0, as it does not know that at
#' # around x=24 density is high
#' library(ggplot2)
#' ggplot(data = data.frame(y = data), aes(x = y)) +
#'   geom_density()
#'
#' # circular density estimation is aware that x=0 is almost identical to x=24
#' data_circ <- circular_density(data)
#' ggplot(data = data_circ, aes(x = time, y = density)) +
#'   geom_line()
#'
#' @importFrom stats density
#'
#' @export
#'
circular_density <- function(x, from = 0, to = 24, density_args = list()) {
  if (!is.numeric(x)) stop("x must be numerical")
  ft <- c(from, to)
  if (!is.numeric(ft) | sum(ft < 0) > 0 | sum(is.na(ft)) > 0) stop("from/to must be non-negative & numeric")

  # Set a default bandwidth if none is set via density_args
  if (!is.list(density_args)) stop("density_args must be a list")
  if (!"bw" %in% names(density_args)) {
    density_args$bw <- .5
  }

  # Append data to itself to respect circular nature
  y <- c(x - to, x, x + to)

  h <- do.call(density, c(list(y, from = from, to = to), density_args))
  d <- data.frame(time = h$x, density = h$y)
  d <- d[d$time >= from & d$time <= to, ]

  return(d)
}

#' Convert radians hours to period hours.
#'
#' @param x a radians time value, from 0 to 2*pi.
#' @param period The numeric period.
#'
#' @export
#'
rad2period <- function(x, period = 24) {
  if (x < 0 | x > (2 * pi)) stop("x must be between 0 and 2*pi")
  if (!is.numeric(period) | period < 0) stop("period must be numeric and not be negative")

  period_hours <- (x / (2 * pi)) * period

  return(period_hours)
}

#' Circular distance between two values, given a period.
#'
#' Calculates circular distance and direction with respect to a period.
#' Positive values mean that x is ahead of y in circular space.
#' For example, x=1 and y=23 gives +2 because the circular distance is 2 and
#' on a clock with 24h x is two hours later than y. Likewise,
#' x=23 and y=1 gives negative 2 because distance again is 2 but on a 24h clock x
#' is earlier than y.
#'
#' @param x,y The two values to calculate circular distance for.
#' @param period THe period of the time series
#' 
#' @examples
#' # result is +2 because on a circular clock x is 2h ahead of y
#' circular_distance(x=1, y=23, period=24)
#' 
#' # result is -5 because on a circular clock x is 5h behind y
#' circular_distance(x=15, y=20, period=24)
#'
#' @export
#'
circular_distance <- function(x, y, period = 24) {
  xyp <- c(x, y, period)

  if (!is.numeric(xyp)) stop("x, y, period must be numeric and >= 0")
  if (sum(!is.finite(xyp)) > 0) stop("x, y, period must be finite")
  if (length(xyp) != 3) stop("x, y, period must all have length of 1")
  if (sum(xyp < 0) > 0) stop("x, y, period must not be negative")

  diff <- abs(x - y)
  d <- min(diff, period - diff)

  steps1 <- ((x - y) %% period + period) %% period
  steps2 <- ((y - x) %% period + period) %% period

  z <- if (steps1 <= steps2) d else d * (-1)

  return(z)
}

#' Calculate absolute circular distance between one query value and a vector with one or many values given a fixed period.
#'
#' @param query A numeric value
#' @param reference A vector with numeric values
#' @param period A period for the circle
#'
#' @export
#'
circular_distance_reference <- function(query, reference, period = 24) {
  qr <- c(query, reference)
  qrp <- c(query, reference, period)

  if (!is.numeric(qrp) | sum(!is.finite(qrp)) > 0 | sum(qrp < 0) > 0) stop("query, reference and period must all be numeric and finite and >= 0")
  if (sum(qr < 0) > 0) stop()

  diff <- abs(reference - query)
  d <- pmin(diff, period - diff)

  return(d)
}

#' Circular-aware intersection of a query value with a time window.
#'
#' Intersect a single timepoint with a timepoint window in a circular-aware fashion given a period,
#' returning TRUE or FALSE if there is an intersection or not. The lower limit of the window is inclusive,
#' the upper one is not. Meaning, a query of 2 or a window of 2-4 would return TRUE but for a window 1-2 would
#' return FALSE.
#'
#' @param query The query value.
#' @param window A vector of two values representing the lower and upper bounds of the time window.
#'
#' @examples
#' query <- c(2)
#' window1 <- c(0, 3) # TRUE
#' window2 <- c(23, 3) # TRUE
#' window3 <- c(3, 5) # FALSE
#' window4 <- c(23, 2) # FALSE
#'
#' circular_intersect(query, window1)
#' circular_intersect(query, window2)
#' circular_intersect(query, window3)
#' circular_intersect(query, window4)
#'
#' @export
#'
circular_intersect <- function(query, window) {
  cc <- c(query, window)

  if (!is.numeric(cc) > 0 | sum(!is.finite(cc)) > 0 | sum(cc < 0) > 0) stop("Both query and window must be numeric, finite and non-negative")
  if (length(query) > 1) stop("query must be a single value")
  if (length(window) != 2) stop("window must be of length 2")
  if (length(unique(window)) == 1) stop("The two window values are the same")

  if (window[1] > window[2]) {
    x <- query >= window[1] | query < window[2]
  } else {
    x <- query >= window[1] & query < window[2]
  }

  return(x)
}
