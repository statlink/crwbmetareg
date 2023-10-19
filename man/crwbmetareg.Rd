\name{Meta regression using cluster robust wild bootstrap}
\alias{crwbmetareg}
\title{
Meta regression using cluster robust wild bootstrap
}

\description{
Meta regression using cluster robust wild bootstrap.
}

\usage{
crwbmetareg(target, se, dataset, cluster, weights, boot.reps = 1000,
prog.bar = FALSE, seed = NULL)
}

\arguments{
\item{target}{
A vector with the effect sizes.
}
\item{se}{
A vector with the standard errors, or the variances, of the effect sizes.
}
\item{dataset}{
A matrix or data.frame with the independent variables.
}
\item{cluster}{
A vector indicating the clusters.
}
\item{weights}{
A vector with the inverse of the the variances of the effect sizes.
}
\item{boot.reps}{
The number of bootstrap re-samples to generate.
}
\item{prog.bar}{
If you want the progress bar to appear set this equal to TRUE.
}
\item{seed}{
IF you want the results to be rerpoducible set this equal to TRUE.
}
}

\details{
It implements metaregression using cluster robust wild bootstrap to compute the p-values.
See references for this.

The function uses a modification of the function "cluster.wild.glm()" of the package "clusterSEs".
}

\value{
A vector with two p-values. One for the constant and one for the cofficient of the "se".
}


\references{
Oczkowski, E. and Doucouliagos, H. (2015).
Wine prices and quality ratings: a meta-regression analysis.
American Journal of Agricultural Economics, 97(1): 103--121.

Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2008).
Bootstrap-based improvements for inference with clustered errors.
The Review of Economics and Statistics, 90(3): 414--427.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{fatpet} }
}

\examples{
y <- rnorm(50)
se <- rexp(50, 3)
cluster <- sample(1:20, 50, replace = TRUE)
dataset <- matrix( rnorm(50 * 2), ncol = 2 )
fatpet(y, se, dataset, cluster, weights = se^2, boot.reps = 100)
}

