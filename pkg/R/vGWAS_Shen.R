`brown.forsythe.test` <-
function(y, group, kruskal.test = FALSE)
{
 # ----- stop the code if the length of y does not match the length of group ----- #
 if (length(y) != length(group))
 {
  stop("The length of the data does not match the length of the group.")
 }
 # ----- assign stuffs ----- #
 DNAME <- deparse(substitute(y))
 y <- y[!is.na(y)]
 group <- group[!is.na(y)] 
 # ----- sort the order just in case the input is not sorted by group ----- # 
 reorder <- order(group)
 group <- group[reorder]
 y <- y[reorder]
 gr <- group
 group <- as.factor(group) # precautionary
 # ----- define the measure of central tendency (median) ----- # 
 means <- tapply(y, group, median)
 METHOD <- "Brown-Forsythe test based on the absolute deviations from the median"
 # ----- calculate the sample size of each group and absolute deviation from center ----- #  
 n <- tapply(y, group, length)
 resp.mean <- abs(y - means[group])
 ngroup <- n[group] 
 # ----- set d ----- # 
 d <- group
 # ----- if the Kruskal-Wallis test is not used ----- #
 if (kruskal.test == FALSE)
 {
  statistic <- anova(lm(resp.mean ~ d))[1, 4]
  p.value <- anova(lm(resp.mean ~ d))[1, 5]
 }
 # ----- if the Kruskal-Wallis test is used ----- #
 else
 {
  METHOD <- paste("Rank-based (Kruskal-Wallis)", METHOD)
  ktest <- kruskal.test(resp.mean,d)
  statistic <- ktest$statistic
  p.value <- ktest$p.value
 } 
 # ----- display output ----- #
 STATISTIC <- statistic
 names(STATISTIC) = "Test Statistic" 
 structure(list(statistic = STATISTIC, p.value = p.value, method = METHOD, data.name = DNAME), class = "htest")
}





`vGWAS` <- 
function(phenotype, geno.matrix, kruskal.test = FALSE, marker.map = NULL, chr.index = NULL)
{
 Call <- match.call()
 n <- length(phenotype)
 m <- ncol(geno.matrix)
 # ----- check phenotype ----- #
 if (!is.numeric(phenotype) & !is.logical(phenotype))
 {
  stop('phenotype has to be numeric or logical.')
 }
 # ----- check genotypes ----- #
 if (!is.matrix(geno.matrix) & !is.data.frame(geno.matrix))
 {
  stop('geno.matrix has to be a matrix or a data frame.')
 }
 # ----- check if data sizes match ----- #
 if (n != nrow(geno.matrix))
 {
  stop('size of phenotype and geno.matrix do not match.')
 }
 if (!is.null(chr.index))
 {
  if (m != length(chr.index))
  {
   stop('size of chr.index and geno.matrix do not match.')
  }
 }
 else
 {
  chr.index <- rep(1, m)
 }
 if (!is.null(marker.map))
 {
  if (m != length(marker.map))
  {
   stop('size of marker.map and geno.matrix do not match.')
  }
 }
 else
 {
  tab.chr <- table(chr.index)
  marker.map <- c()
  for (i in 1:length(tab.chr)) {
   marker.map <- c(marker.map, 1:tab.chr[i])
  }
 }
 # ----- preallocation ----- #
 p.values <- statistics <- numeric(m)
 pb <- txtProgressBar(style = 3)
 # ----- scan using Brown-Forsythe test -----#
 for (j in 1:m)
 {
  test <- try(brown.forsythe.test(phenotype, as.factor(geno.matrix[,j]), kruskal.test = kruskal.test), silent  = TRUE)
  if (!inherits(test, 'try-error'))
  {
   p.values[j] <- test$p.value
   statistics[j] <- test$statistic
  }
  else
  {
   p.values[j] <- 1
   statistics[j] <- 0
  }
  setTxtProgressBar(pb, j/m)
 }
 cat('\n')
 marker.names <- names(as.data.frame(geno.matrix))
 res <- data.frame(marker = marker.names, chromosome = chr.index, marker.map = marker.map, statistic = statistics, p.value = p.values)
 class(res) <- 'vGWAS'
 return(res)
}





`plot.vGWAS` <-
function(x, sig.threshold = NULL, low.log.p = 0, pch = 16, cex = .6, col.manhattan = c('slateblue4', 'olivedrab'), col.sig.threshold = 'darkgoldenrod', ...)
{
 tab.chr <- table(x$chromosome)
 chr <- as.numeric(names(tab.chr))
 ends <- cumsum(tab.chr)
 cumpos <- numeric(length(x$marker.map))
 cumpos[1:ends[1]] <- x$marker.map[1:ends[1]]
 logp <- -log(x$p.value, 10)
 if (length(chr) > 1)
 {
  for (i in 2:length(chr)) 
  {
   cumpos[(ends[i - 1] + 1):ends[i]] <- cumpos[ends[i - 1]] + x$marker.map[(ends[i - 1] + 1):ends[i]] 
  }
 }
 if (is.null(sig.threshold)) 
 {
  sig.threshold <- -log(.05/length(x$marker.map), 10)
  cat('nominal significance threshold with Bonferroni correction for', length(x$marker.map), 'tests are calculated.\n')
 }
 cutp <- logp > low.log.p
 plot(cumpos, logp, type = 'n', ann = FALSE, axes = FALSE)
 at1 <- cumpos[1]
 points(cumpos[(1:ends[1])[cutp[1:ends[1]]]], logp[(1:ends[1])[cutp[1:ends[1]]]], pch = pch, cex = cex, col = col.manhattan[1])
 at1 <- c(at1, cumpos[ends[1]])
 at2 <- (cumpos[1] + cumpos[ends[1]])/2
 if (length(chr) > 1)
 {
  for (i in 2:length(chr)) 
  {
   points(cumpos[((ends[i - 1] + 1):ends[i])[cutp[(ends[i - 1] + 1):ends[i]]]], logp[((ends[i - 1] + 1):ends[i])[cutp[(ends[i - 1] + 1):ends[i]]]], pch = pch, cex = cex, col = col.manhattan[2 - i%%2])
   at1 <- c(at1, cumpos[ends[i]])
   at2 <- c(at2, (cumpos[ends[i - 1]] + cumpos[ends[i]])/2)
  }
 }
 abline(h = sig.threshold, lty = 2, col = col.sig.threshold)
 axis(1, at = at1, labels = FALSE)
 mtext(chr, 1, 0, at = at2)
 mtext('Chromosome', 1, 2)
 axis(2)
 mtext(expression(-log[10]~'('~italic(P)~-value~')'), 2, 3)
}





`vGWAS.heritability` <- 
function(phenotype, marker.genotype, only.print = TRUE)
{
 # ----- check phenotype ----- #
 if (!is.numeric(phenotype))
 {
  stop('phenotype has to be numeric.')
 }
 # ----- check marker genotype ----- #
 if (!is.numeric(marker.genotype) & !is.character(marker.genotype) & !is.factor(marker.genotype))
 {
  stop('marker genotype has a wrong format.')
 }
 # ----- check if data sizes match ----- #
 if (length(phenotype) != length(marker.genotype))
 {
  stop('size of phenotype and marker genotype do not match.')
 }
 dm <- dglm(phenotype ~ as.factor(marker.genotype), ~ as.factor(marker.genotype))
 heritability.mean <- (dm$null.deviance - dm$deviance)/dm$null.deviance
 heritability.disp <- (dm$dispersion.fit$null.deviance - dm$dispersion.fit$deviance)/dm$dispersion.fit$null.deviance
 cat('heritability explained by the mean part of model:\n')
 cat(round(heritability.mean*100, digits = 2), '%\n')
 cat('heritability explained by the variance part of model:\n')
 cat(round(heritability.disp*100, digits = 2), '%\n')
 cat('heritability in total:\n')
 cat(round((heritability.mean + heritability.disp)*100, digits = 2), '%\n')
 if (!only.print) return(list(heritability.mean = heritability.mean, heritability.disp = heritability.disp))
}





.onAttach <- 
function(...)
{
 cat("vGWAS: Variance Genome-wide Association\n")
 cat('Version 2010.10.12 installed\n')
 cat('Correspondence to: Xia Shen (xia.shen@lcb.uu.se)\n')
}