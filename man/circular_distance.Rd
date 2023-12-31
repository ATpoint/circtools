% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{circular_distance}
\alias{circular_distance}
\title{Circular distance between two values, given a period.}
\usage{
circular_distance(x, y, period = 24)
}
\arguments{
\item{x, y}{The two values to calculate circular distance for.}

\item{period}{THe period of the time series}
}
\description{
Calculates circular distance and direction with respect to a period.
Positive values mean that x is ahead of y in circular space.
For example, x=1 and y=23 gives +2 because the circular distance is 2 and
on a clock with 24h x is two hours later than y. Likewise,
x=23 and y=1 gives negative 2 because distance again is 2 but on a 24h clock x
is earlier than y.
}
\examples{
# result is +2 because on a circular clock x is 2h ahead of y
circular_distance(x=1, y=23, period=24)

# result is -5 because on a circular clock x is 5h behind y
circular_distance(x=15, y=20, period=24)

}
