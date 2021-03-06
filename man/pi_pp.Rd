% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pi_pp.R
\name{pi_pp}
\alias{pi_pp}
\title{Probability of Passing PPQ Test using Prediction Interval}
\usage{
pi_pp(Llim, Ulim, mu, sigma, n, n.batch, alpha)
}
\arguments{
\item{Llim}{lower specification limit}

\item{Ulim}{upper specification limit}

\item{mu}{hypothetical mean of the attribute}

\item{sigma}{hypothetical standard deviation of the attribute}

\item{n}{sample size (number of locations) per batch}

\item{n.batch}{number of batches for passing PPQ during validation}

\item{alpha}{significant level for constructing the prediction interval.}
}
\value{
A numeric value of the passing/acceptance probability
}
\description{
The function for calculating the probability of passing critical quality attributes (CQA) PPQ test .
}
\examples{
\dontrun{
pi_pp(sigma=0.5, mu=2.5, n=10, n.batch=1, Llim=1.5, Ulim=3.5, alpha=0.05)

sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = pi_pp, mu=97, n=10, Llim=95, Ulim=105,
n.batch=1, alpha=0.05)
sapply(X=c(0.1,0.5, 1,2,3,4,5,10), FUN = pi_pp, mu=100, n=10, Llim=95, Ulim=105,
n.batch=1, alpha=0.05)
}
}
\references{
Burdick, R. K., LeBlond, D. J., Pfahler, L. B., Quiroz, J., Sidor, L., Vukovinsky, K., & Zhang, L. (2017).
Statistical Applications for Chemistry, Manufacturing and Controls (CMC) in the Pharmaceutical Industry.
\emph{Springer}.
}
\seealso{
\code{rl_pp}.
}
\author{
Yalin Zhu
}
