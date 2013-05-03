\name{plot.vGWAS}
\alias{plot.vGWAS}
\title{
Variance GWA Manhattan Plot
}
\description{
The function plots the variance GWA result for the given scan object.
}
\usage{
\method{plot}{vGWAS}(x, sig.threshold = NULL, low.log.p = 0, pch = 16, 
           cex = 0.6, col.manhattan = c("slateblue4", "olivedrab"), 
           col.sig.threshold = "darkgoldenrod", ...)
}
\arguments{
  \item{x}{
a result object from \code{vGWAS} scan. It can be any \code{list} or \code{data.frame} that contains \code{chromosome}, \code{marker.map}, and \code{p.value}, with \code{class = 'vGWAS'}. See \code{\link{vGWAS}}.
}
  \item{sig.threshold}{
a numeric value giving the significance threshold for \code{-log(pvalues, 10)}. If \code{NULL}, Bonferroni correction will be used.
}
  \item{low.log.p}{
a numeric value giving the lower limit of the \code{-log(pvalues, 10)} to plot.
}
  \item{pch}{
point character. See \code{\link{par}}.
}
  \item{cex}{
size of points. See \code{\link{par}}.
}
  \item{col.manhattan}{
two colors as a vector for the Manhattan plot.
}
  \item{col.sig.threshold}{
one color for the significance threshold.
}
  \item{...}{
not used.
}
}
\value{
a plot for viewing vGWAS result.
}
\references{
Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O. (2011): \bold{Inheritance beyond plain heritability: variance-controlling genes in \emph{Arabidopsis thaliana}}. \emph{Submitted}.\cr
}
\author{
Xia Shen
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