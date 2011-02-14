\name{vGWAS.heva}
\alias{vGWAS.heva}
\title{
Removing Population Confounding from Phenotype
}
\description{
The function corrects the population stratification using hierarchical generalized linear model (HGLM).
}
\usage{
vGWAS.heva(phenotype, geno.matrix = NULL, kinship = NULL, 
           family = gaussian())
}
\arguments{
  \item{phenotype}{
a \code{numeric} or \code{logical} vector of the phenotyic values. See \bold{Examples}.
}
  \item{geno.matrix}{
a \code{matrix} or \code{data.frame} with individuals as rows and markers as columns. The marker genotypes for each marker are coded as one column. See \bold{Examples}.
}
  \item{kinship}{
a \code{matrix} with size number of individuals times number of individuals, giving the kinship between each pair of the individuals, e.g. IBD or IBS matrix. Only useful when \code{heva = TRUE}. If \code{NULL}, a simple genomic kinship matrix is created using the genotype matrix for HEVA. See \bold{Examples}.
}
  \item{family}{
a \code{family} function specifying the distribution of \code{phenotype}. See the depended R package \bold{hglm} for more information.
}
}
\value{
a \code{numeric} vector of corrected phenotype on a continuous scale, which is the studentized deviance residuals of HGLM. See the depended R package \bold{hglm} for more information.
}
\references{
Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O. (2011): \bold{Variance Controlling Genes: Significant Contributors to The Missing Genetic Variation for Complex Traits}. \emph{Submitted}.\cr

Ronnegard, L., Shen, X. and Alam, M. (2011): \bold{hglm: A Package for Fitting Hierarchical Generalized Linear Models}. \emph{The R Journal}, \bold{2}(2), 20-28.\cr
}
\author{
Xia Shen
}
\seealso{
\code{\link{vGWAS}}, \code{\link{vGWAS-package}}
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

# ----- vGWAS via HEVA ----- #
 
vgwa2 <- vGWAS(phenotype = pheno, geno.matrix = geno, heva = TRUE,
              marker.map = map, chr.index = chr)
              
plot(vgwa2)
}
}
