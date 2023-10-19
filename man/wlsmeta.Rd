\name{Weighted least squares meta analysis}
\alias{wlsmeta}
\title{
Weighted least squares meta analysis
}

\description{
Weighted least squares meta analysis.
}

\usage{
wlsmeta(yi, vi)
}

\arguments{
\item{yi}{
The observations.
}
\item{vi}{
The variances of the observations.
}
}

\details{
It implements weighted least squares (WLS) meta analysis. See references for this.
}

\value{
A vector with many elements. The fixed effects mean estimate, the \eqn{\bar{v}}
estimate, the \eqn{I^2}, the \eqn{H^2}, the Q test statistic and it's p-value,
the \eqn{\tau^2} estimate and the random effects mean estimate.
}


\references{
Stanley T. D. and Doucouliagos H. (2015).
Neither fixed nor random: weighted least squares meta-analysis.
Statistics in Medicine, 34(13): 2116--2127.

Stanley, T. D. and Doucouliagos, H. (2017).
Neither fixed nor random: Weighted least squares meta-regression.
Research synthesis methods, 8(1): 19--42.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colwlsmeta} }
}

\examples{
y <- rnorm(30)
vi <- rexp(30, 3)
wlsmeta(y, vi)
}

