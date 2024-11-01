% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_selection.R
\name{run_model_selection}
\alias{run_model_selection}
\title{Classify rhythmicity between two groups using model selection as suggested in the compareRhythms package.}
\usage{
run_model_selection(
  x,
  condition,
  time,
  period = 24,
  criterion = c("BIC", "AIC"),
  sample_weights = NULL
)
}
\arguments{
\item{x}{a numeric matrix or data frame where rows are observations and columns are samples}

\item{condition}{a character or factor vector of length \code{ncol(x)} that contains the condition membership information}

\item{time}{vector of length \code{ncol(x)} that contains the time information}

\item{period}{Period for the timecourse, typically 24 for circadian data}

\item{criterion}{BIC (Bayesian) or AIC (Akaike) information criterion to select the best model fit}

\item{sample_weights}{Optional vector of length \code{ncol(x)} with sample weights}
}
\description{
The function accepts a observation by sample matrix with circadian data and then uses a model fitting approach
inspired and adopted from the \url{https://github.com/bharathananth/compareRhythms} package, to classify each gene
into a category and calculate a probability score. The categories are "change", "gain", "loss", and "same",
and the probability describes whether the best model is significantly #' better than the 2nd-best model.
The compareRhythms authors suggest a cutoff of 0.6 and higher to "trust" the classification.
}
\details{
Since the available models do not contain an "arrhythmic model" where both conditions are entirely arrhythmic,
it is recommended to do a preselection to ensure that at least one condition per gene shows significant cosinor-defined rhythmicity,
for example using the \code{run_cosinor()} function in this package.

Input data should be normalized and transformed (for example log2) to be suitable for linear models.
}
\examples{
y <- rbind(
  # amplitude change
  make_data(time = seq(1, 21, 4), amplitude = 2),
  make_data(time = seq(1, 21, 4), amplitude = 1),
  # acrophase change
  make_data(time = seq(1, 21, 4), acrophase = 1),
  make_data(time = seq(1, 21, 4), acrophase = 13),
  # same rhythmicity
  make_data(time = seq(1, 21, 4), acrophase = 1),
  make_data(time = seq(1, 21, 4), acrophase = 1),
  # loss of rhythmicity
  make_data(time = seq(1, 21, 4), amplitude = 1),
  make_data(time = seq(1, 21, 4), amplitude = 0),
  # gain of rhythmicity
  make_data(time = seq(1, 21, 4), amplitude = 0),
  make_data(time = seq(1, 21, 4), amplitude = 1)
)
x <- matrix(y$value, byrow = TRUE, ncol = 12)
rownames(x) <- paste0("gene_", 1:nrow(x))
time <- y$time[1:12]
condition <- rep(c("A", "B"), each = 6)
run_model_selection(x = x, condition = condition, time = time)
}
