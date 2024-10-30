#' Classify rhythmicity between two groups using model selection as suggested in the compareRhythms package.
#'
#' The function accepts a observation by sample matrix with circadian data and then uses a model fitting approach
#' inspired and adopted from the \url{https://github.com/bharathananth/compareRhythms} package, to classify each gene
#' into a category and calculate a probability score. The categories are "change", "gain", "loss", and "same",
#' and the probability describes whether the best model is significantly #' better than the 2nd-best model.
#' The compareRhythms authors suggest a cutoff of 0.6 and higher to "trust" the classification.
#'
#' Since the available models do not contain an "arrhythmic model" where both conditions are entirely arrhythmic,
#' it is recommended to do a preselection to ensure that at least one condition per gene shows significant cosinor-defined rhythmicity,
#' for example using the \code{run_cosinor()} function in this package.
#'
#' Input data should be normalized and transformed (for example log2) to be suitable for linear models.
#'
#' @param x a numeric matrix or data frame where rows are observations and columns are samples
#' @param condition a character or factor vector of length \code{ncol(x)} that contains the condition membership information
#' @param time vector of length \code{ncol(x)} that contains the time information
#' @param period Period for the timecourse, typically 24 for circadian data
#' @param criterion BIC (Bayesian) or AIC (Akaike) information criterion to select the best model fit
#' @param sample_weights Optional vector of length \code{ncol(x)} with sample weights
#'
#' @importFrom stats model.matrix
#' @importFrom limma selectModel
#'
#' @examples
#' y <- rbind(
#'   # amplitude change
#'   make_data(time = seq(1, 21, 4), amplitude = 2),
#'   make_data(time = seq(1, 21, 4), amplitude = 1),
#'   # acrophase change
#'   make_data(time = seq(1, 21, 4), acrophase = 1),
#'   make_data(time = seq(1, 21, 4), acrophase = 13),
#'   # same rhythmicity
#'   make_data(time = seq(1, 21, 4), acrophase = 1),
#'   make_data(time = seq(1, 21, 4), acrophase = 1),
#'   # loss of rhythmicity
#'   make_data(time = seq(1, 21, 4), amplitude = 1),
#'   make_data(time = seq(1, 21, 4), amplitude = 0),
#'   # gain of rhythmicity
#'   make_data(time = seq(1, 21, 4), amplitude = 0),
#'   make_data(time = seq(1, 21, 4), amplitude = 1)
#' )
#' x <- matrix(y$value, byrow = TRUE, ncol = 12)
#' rownames(x) <- paste0("gene_", 1:nrow(x))
#' time <- y$time[1:12]
#' condition <- rep(c("A", "B"), each = 6)
#' run_model_selection(x = x, condition = condition, time = time, add_arrhy = TRUE)
#' @export
run_model_selection <- function(x, condition, time, period = 24,
                                criterion = c("BIC", "AIC"),
                                sample_weights = NULL, add_arrhy = FALSE) {
  # Checks
  if (!ncol(x) == length(condition) | !ncol(x) == length(time)) {
    stop("condition and time must have length as ncol(x)")
  }

  if (!is.null(sample_weights)) {
    if (ncol(x) != length(sample_weights)) {
      stop("sample_weights must have length as ncol(x)")
    }
  }

  criterion <- match.arg(criterion, choices = c("BIC", "AIC"))

  if (!is.numeric(period) | period < 0) {
    stop("cutoff must be numeric and >= 0")
  }

  w <- if (!is.null(sample_weights)) sample_weights else NULL
  v <- new("EList", list(E = x, weights = w))

  condition <- factor(condition)
  metadata <- data.frame(
    time = time, condition = condition,
    cosphase = cos(2 * pi * time / period),
    sinphase = sin(2 * pi * time / period)
  )

  design_list <- list()

  # The full or "change" model so both conditions are rhythmic but parameters acrophase/amplitude changes
  design_list$change <- model.matrix(~ condition + condition:sinphase + condition:cosphase, data = metadata)

  # The "same" model, no difference, meaning no interaction terms
  design_list$same <- model.matrix(~ condition + sinphase + cosphase, data = metadata)

  # The loss and gain models so rhythmic parameters only for one group or the other.
  lvls <- levels(condition)
  design_list$loss <- design_list$change[, !grepl(paste0(lvls[2], ":"), colnames(design_list$change))]
  design_list$gain <- design_list$change[, !grepl(paste0(lvls[1], ":"), colnames(design_list$change))]

  # Calculate the conditional probability as in compareRhythms
  sm <- selectModel(y = v, designlist = design_list, criterion = tolower(criterion))
  prob <- apply(sm$IC, 1, function(x) base::max(exp(-0.5 * x) / sum(exp(-0.5 * x))))

  summarized <- data.frame(
    row.names = rownames(x),
    model_category = sm$pref,
    model_probability = prob
  )

  attr(summarized, "criterion") <- criterion
  summarized$model_probability[is.nan(summarized$model_probability)] <- NA

  return(summarized)
}
