\name{binomCI}
\alias{binomCI}

\title{
Confidence Intervals for a Binomial Proportion.
}

\description{
Confidence Intervals for a Binomial Proportion.
}

\usage{
binomCI(x, n, a = 0.05)
}

\arguments{
\item{x}{
The number of successes.
}
\item{n}{
The number of trials.
}
\item{a}{
The significance level to compute the \eqn{(1-\alpha)\%} confidence intervals.
}
}

\details{
The confidence intervals are:

\emph{Jeffreys}:
\deqn{
\left[ F(\alpha/2; x+0.5, n-x+0.5), F(1-\alpha/2; x+0.5, n-x+0.5) \right],
}
where \eqn{F(\alpha, a, b)} denotes the \eqn{\alpha} quantile of the Beta distribution with parameters \eqn{a} and \eqn{b}, \eqn{Be(a, b)}.

\emph{Wald}:
\deqn{
\left[ \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{p}(1-\hat{p})}{n}}, \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{p}(1-\hat{p})}{n}} \right],
}
where \eqn{\hat{p}=\frac{x}{n}} and \eqn{Z_{1-\alpha/2}} denotes the \eqn{1-\alpha/2} quantile of the standard normal distribution. If \eqn{\hat{p}=0} the interval becomes \eqn{(0 , 1 - e^{\frac{1}{n}\log({\alpha}{2})})} and if \eqn{\hat{p}=1} the interval becomes \eqn{(e^{\frac{1}{n}\log({\alpha}{2})}, 1)}.

\emph{Wald corrected}:
\deqn{
\left[ \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{p}(1-\hat{p})}{n}} - \frac{0.5}{n}, \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{p}(1-\hat{p})}{n}} + \frac{0.5}{n} \right],
}
and if \eqn{\hat{p}=0} or \eqn{\hat{p}=1} the previous (Wald) adjustment applies.

\emph{Wald BS}:
\deqn{
\left[ \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{p}(1-\hat{p})}{n-Z_{1-\alpha/2}-2Z_{1-\alpha/2}/n-1/n}} - \frac{0.5}{n}, \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{p}(1-\hat{p})}{n-Z_{1-\alpha/2}-2Z_{1-\alpha/2}/n-1/n}} + \frac{0.5}{n} \right],
}
and if \eqn{\hat{p}=0} or \eqn{\hat{p}=1} the previous (Wald) adjustment applies.

\emph{Agresti and Coull}:
\deqn{
\left[ \hat{\theta} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{\theta}(1-\hat{\theta})}{n+4}}, \hat{p} - Z_{1-\alpha/2} \times \sqrt{\frac{\hat{\theta}(1-\hat{\theta})}{n+4}} \right],
}
where \eqn{\hat{\theta}=\frac{x+2}{n+4}}.

\emph{Wilson}:
\deqn{
\left[ \frac{x_b}{n_b} - \frac{Z_{1-\alpha/2}\sqrt{n}}{n_b} \times \sqrt{\hat{p}(1-\hat{p})+Z_{1-\alpha/2}/4}, \frac{x_b}{n_b} + \frac{Z_{1-\alpha/2}\sqrt{n}}{n_b} \times \sqrt{\hat{p}(1-\hat{p})+Z_{1-\alpha/2}/4} \right],
}
where \eqn{x_b=x+Z_{1-\alpha/2}^2/2} and \eqn{n_b=n+Z_{1-\alpha/2}^2}.

\emph{Score}:
\deqn{
\left[ \frac{x+Z_{1-\alpha/2}^2-c}{n+Z_{1-\alpha/2}^2} , \frac{x+Z_{1-\alpha/2}^2+c}{n+Z_{1-\alpha/2}^2} \right],
}
where \eqn{c=Z_{1-\alpha/2}\sqrt{x-x^2/n+Z_{1-\alpha/2}^2/4}}.

\emph{Score corrected}:
\deqn{
\left[ \frac{\ell_1}{n+Z_{1-\alpha/2}} , \frac{\ell_2}{n+Z_{1-\alpha/2}} \right],
}
where \eqn{\ell_1=b_1+0.5Z_{1-\alpha/2}^2-Z_{1-\alpha/2}\sqrt{b_1-b_1^2/n+0.25Z_{1-\alpha/2}^2}}, \eqn{\ell_2=b_2+0.5Z_{1-\alpha/2}^2+Z_{1-\alpha/2}\sqrt{b_2-b_2^2/n+0.25Z_{1-\alpha/2}^2}} and \eqn{b_1=x-0.5}, \eqn{b_2=x+0.5}.

\emph{Wald-logit}:
\deqn{
\left[ 1-(1+e^{b-c})^{-1}, 1-(1+e^{b+c})^{-1} \right],
}
where \eqn{b=\log(\frac{x}{n-x})} and \eqn{c=\frac{Z_{1-\alpha/2}}{\sqrt{n\hat{p}(1-\hat{p})}}}. If \eqn{\hat{p}=0} or \eqn{\hat{p}=1} the previous (Wald) adjustment applies.

\emph{Wald-logit corrected}:
\deqn{
\left[ 1-(1+e^{b-c})^{-1}, 1-(1+e^{b+c})^{-1} \right],
}
where \eqn{b=\log(\frac{\hat{p}_b}{\hat{q}_b})}, \eqn{\hat{p}_b=x+0.5}, \eqn{\hat{q}_b=n-x+0.5} and \eqn{c=\frac{Z_{1-\alpha/2}}{\sqrt{(n+1)\frac{\hat{p}_b}{n+1}(1-\frac{\hat{p}_b}{n+1})}}}.

\emph{Arcsine}:
\deqn{
\left\lbrace \sin^2\left[sin^{-1}(\sqrt{\hat{p}})-0.5\frac{Z_{1-\alpha/2}}{\sqrt{n}}\right], \sin^2\left[sin^{-1}(\sqrt{\hat{p}})+0.5\frac{Z_{1-\alpha/2}}{\sqrt{n}}\right] \right\rbrace.
}
If \eqn{\hat{p}=0} or \eqn{\hat{p}=1} the previous (Wald) adjustment applies.

\emph{Exact binomial}:
\deqn{
\left[ (1+\frac{a_1}{d_1})^{-1}, (1+\frac{a_2}{d_2})^{-1} \right],
}
where \eqn{a_1=n-x+1}, \eqn{a_2=a_1-1}, \eqn{d_1=x-F(\alpha/2,2x,2a_1)}, \eqn{d_2=(x+1)F(1-\alpha/2,2(x+1),2a_2)} and \eqn{F(\alpha,a,b)} denotes the \eqn{\alpha} quantile of the F distribution with degrees of freedom \eqn{a} and \eqn{b}, \eqn{F(a, b)}.
}

\value{
A list including:
\item{prop}{The proportion.
}
\item{ci}{
A matrix with 12 rows containing the 12 different \eqn{(1-\alpha)\%} confidence intervals.
}
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{binomCIs}
}
}

\examples{
binomCI(45, 100)
}


