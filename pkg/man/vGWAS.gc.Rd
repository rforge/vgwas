\name{vGWAS.gc}
\alias{vGWAS.gc}
\title{
Genomic Control for vGWAS
}
\description{
The function does genomic control for the variance GWA result object.
}
\usage{
vGWAS.gc(object, plot = TRUE, proportion = 1, ...)
}
\arguments{
  \item{object}{
a result object from \code{vGWAS} scan. It can be any \code{list} or \code{data.frame} that contains \code{p.value}. See \code{\link{vGWAS}}.
}
  \item{plot}{
a logical value turning on/off the QQ plot for genomic control.
}
  \item{proportion}{
a numeric value between 0 and 1 giving the proportion of obtained p-values to be used for genomic control.
}
  \item{...}{
not used.
}
}
\value{
\item{lambda}{
estimated inflation ratio.
}
\item{lambda.se}{
standard error of the estimated inflation ratio.
}
\item{gc.p.value}{
p-values after genomic control.
}
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

# ----- genomic control ----- #

vgwa2 <- vGWAS.gc(vgwa)
                                
plot(vgwa2) 
}
}