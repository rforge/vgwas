\name{vGWAS}
\alias{vGWAS}
\title{
Variance Genome-wide Association
}
\description{
Variance Genome-wide association for using nonparametric variance test
}
\usage{
vGWAS(phenotype, geno.matrix, kruskal.test = FALSE, 
      marker.map = NULL, chr.index = NULL)
}
\arguments{
  \item{phenotype}{
a \code{numeric} or \code{logical} vector of the phenotyic values. See \bold{Examples}.
}
  \item{geno.matrix}{
a \code{matrix} or \code{data.frame} with individuals as rows and markers as columns. The marker genotypes for each marker are coded as one column. See \bold{Examples}.
}
  \item{kruskal.test}{
a \code{logical} value specifying whether to use Kruskal-Wallis statistic. The default option is \code{FALSE}, i.e., the usual ANOVA statistic is used in place of Kruskal-Wallis statistic.
}
  \item{marker.map}{
a \code{numeric} vector giving the marker map positions for each chromosome. See \bold{Examples}.
}
  \item{chr.index}{
a \code{numeric} vector giving the chromosome index for each marker. See \bold{Examples}.
}
}
\value{
a \code{data.frame} containing columns of \code{marker} names, \code{chromosome} indices, \code{marker.map} positions, test \code{statistic} values, and \code{p.value} for each position.
}
\references{
Shen, X., Pettersson, M., Ronnegard, L. and Carlborg, O. (2011): \bold{Inheritance beyond plain heritability: variance-controlling genes in \emph{Arabidopsis thaliana}}. \emph{Submitted}.\cr

Ronnegard, L., Shen, X. and Alam, M. (2010): \bold{hglm: A Package for Fitting Hierarchical Generalized Linear Models}. \emph{The R Journal}, \bold{2}(2), 20-28.\cr
}
\author{
Xia Shen
}
\seealso{
\code{\link{vGWAS-package}}
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