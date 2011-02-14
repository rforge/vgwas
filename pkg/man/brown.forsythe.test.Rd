\name{brown.forsythe.test}
\alias{brown.forsythe.test}
\title{Brown-Forsythe's Test of Equality of Variances}
\description{
The function performs the robust Brown-Forsythe test using the group medians. Instead of the ANOVA statistic, the Kruskal-Wallis ANOVA may also be applied using this function.
}
\usage{
brown.forsythe.test(y, group, kruskal.test=FALSE)
}
\arguments{
  \item{y}{a numeric vector of data values.}
  \item{group}{factor of the data.}
  \item{kruskal.test}{a \code{logical} value specifying whether to use Kruskal-Wallis statistic. The default option is \code{FALSE}, i.e., the usual ANOVA statistic is used in place of Kruskal-Wallis statistic.}
}
\details{
Levene (1960) proposed a test for homogeneity of variances in \emph{k} groups which is based on the ANOVA statistic applied to absolute deviations of observations from the corresponding group mean. The robust Brown-Forsythe version of the Levene-type test substites the group mean by the group median in the classical Levene statistic.
}
\value{
A list with the following numeric components.
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{type of test performed.}
  \item{data.name}{a character string giving the name of the data.}
}

\note{Modified from the \code{lawstat} package.}

\references{
Brown, M. B. and Forsythe, A.B. (1974). \bold{Robust tests for equality of 
variances.} \emph{Journal of the American Statistical Association}, \bold{69}, 364-367.\cr

Levene, H. (1960). \bold{Robust Tests for Equality of Variances}, \emph{in Contributions 
to Probability and Statistics}, ed. I. Olkin, Palo Alto, CA: Stanford Univ. Press.\cr
}

\author{Xia Shen}

\section{Acknowledgement}{The authors of package \code{lawstat} is ackownledged for their source code under free GPL license.}

\examples{
\dontrun{

data(pheno)
data(geno)
brown.forsythe.test(pheno, geno[,911])

}
}
\keyword{htest}
