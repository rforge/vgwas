\name{vGWAS.heritability}
\alias{vGWAS.heritability}
\title{
Calculating Heritability for A Single Marker
}
\description{
The function calculates and reports the heritability for a single marker by fitting a double generalized linear model. It gives both the heritability from the mean and variance parts of model.
}
\usage{
vGWAS.heritability(phenotype, marker.genotype, only.print = TRUE)
}
\arguments{
  \item{phenotype}{
a \code{numeric} vector of the phenotyic values. See \bold{Examples}.
}
  \item{marker.genotype}{
a \code{numeric} or \code{character} or \code{factor} vector of the genotypes of a single marker. See \bold{Examples}.
}
  \item{only.print}{
a \code{logical} value. If \code{FALSE}, the heritability values will be returned for storage.
}
}
\details{
The \bold{Value} will only be available if \code{only.print = FALSE}.
}
\value{
\item{heritability.mean}{the heritability from the mean part of model.}
\item{heritability.disp}{the heritability from the variance part of model.}
}
\references{
Shen et al. (2010)
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

# ----- calculate the heritability of strongest the marker ----- #

vGWAS.heritability(phenotype = pheno, 
                   marker.genotype = geno[,vgwa$p.value == min(vgwa$p.value)])
                   
}
}