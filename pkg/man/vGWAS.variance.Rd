\name{vGWAS.variance}
\alias{vGWAS.variance}
\title{
Calculating Variance Explained by A Single Marker
}
\description{
The function calculates and reports the variance explained for a single marker by fitting a double generalized linear model. It gives both the variance explained by the mean and variance parts of model.
}
\usage{
vGWAS.variance(phenotype, marker.genotype, print = TRUE)
}
\arguments{
  \item{phenotype}{
a \code{numeric} vector of the phenotyic values. See \bold{Examples}.
}
  \item{marker.genotype}{
a \code{numeric} or \code{character} or \code{factor} vector of the genotypes of a single marker. See \bold{Examples}.
}
  \item{print}{
a \code{logical} value. If \code{FALSE}, the heritability values will be returned for storage.
}
}
\details{
The \bold{Value} will only be available if \code{only.print = FALSE}.
}
\value{
\item{variance.mean}{the variance explained by the mean part of model.}
\item{variance.disp}{the variance explained by the variance part of model.}
}
\references{
Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O. (2011): \bold{Inheritance beyond plain heritability: variance-controlling genes in \emph{Arabidopsis thaliana}}. \emph{Submitted}.\cr
}
\author{
Xia Shen
}
\seealso{
\code{\link{vGWAS-package}}, \code{\link{vGWAS}}, \code{\link{plot.vGWAS}}
}
\examples{
\dontrun{

# ----- load data ----- #

data(pheno)
data(geno)
data(chr)
data(map)

# ----- variance GWA scan ----- #

vgwa <- vGWAS(phenotype = pheno, geno.matrix = geno,
              marker.map = map, chr.index = chr)
              
# ----- visualize the scan ----- #

plot(vgwa)

# ----- calculate the variance explained by strongest the marker ----- #

vGWAS.variance(phenotype = pheno, 
                   marker.genotype = geno[,vgwa$p.value == min(vgwa$p.value)])
                   
}
}