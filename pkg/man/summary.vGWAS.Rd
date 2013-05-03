\name{summary.vGWAS}
\alias{summary.vGWAS}
\title{
Variance GWA Summary
}
\description{
The function summarized the variance GWA result for the given scan object.
}
\usage{
\method{summary}{vGWAS}(object, nrMarkers = 10, ...)
}
\arguments{
  \item{object}{
a result object from \code{vGWAS} scan. It can be any \code{list} or \code{data.frame} that contains \code{chromosome}, \code{marker.map}, and \code{p.value}, with \code{class = 'vGWAS'}. See \code{\link{vGWAS}}.
}
  \item{nrMarkers}{
a numeric value giving the number of top markers to be summarized.
}
  \item{...}{
Not in use.
}
}
\value{
a summary for viewing vGWAS result.
}
\references{
Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O. (2011): \bold{Inheritance beyond plain heritability: variance-controlling genes in \emph{Arabidopsis thaliana}}. \emph{Submitted}.\cr
}
\author{
Simon Forsberg, Xia Shen
}
\seealso{
\code{\link{vGWAS-package}}, \code{\link{vGWAS}}
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

summary(vgwa)

# ----- calculate the variance explained by strongest the marker ----- #

vGWAS.heritability(phenotype = pheno, 
                   marker.genotype = geno[,vgwa$p.value == min(vgwa$p.value)])
                   
# ----- genomic control ----- #

vgwa2 <- vGWAS.gc(vgwa)
                                
plot(vgwa2) 

summary(vgwa2)
}
}