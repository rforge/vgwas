`brown.forsythe.test` <-
		function(y, group, kruskal.test = FALSE)
{
	# ----- stop the code if the length of y does not match the length of group ----- #
	if (length(y) != length(group))
	{
		stop('The length of the data does not match the length of the group.')
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
	METHOD <- 'Brown-Forsythe test based on the absolute deviations from the median'
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
		METHOD <- paste('Rank-based (Kruskal-Wallis)', METHOD)
		ktest <- kruskal.test(resp.mean,d)
		statistic <- ktest$statistic
		p.value <- ktest$p.value
	} 
	# ----- display output ----- #
	STATISTIC <- statistic
	names(STATISTIC) = 'Test Statistic' 
	structure(list(statistic = STATISTIC, p.value = p.value, method = METHOD, data.name = DNAME), class = 'htest')
}





`vGWAS` <-
		function(phenotype, geno.matrix, kruskal.test = FALSE, marker.map = NULL, chr.index = NULL)
{
	Call <- match.call()
	# ----- check phenotype ----- #
	if (!is.numeric(phenotype) & !is.logical(phenotype))
	{
		stop('phenotype has to be numeric or logical.')
	}
	#if (heva) 
	#{
	# cat('correcting phenotype using HEVA ...\n')
	# if (is.numeric(phenotype)) family <- gaussian() else family <- binomial()
	# phenotype <- vGWAS.heva(phenotype, geno.matrix, kinship, family = family)$corrected.phenotype
	# cat('phenotype corrected.\n')
	#}
	n <- length(phenotype)
	m <- ncol(geno.matrix)
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
#  if (!heva) {
		test <- try(brown.forsythe.test(phenotype, as.factor(geno.matrix[,j]), kruskal.test = kruskal.test), silent  = TRUE)
#  } else {
#   test <- try(wilcox.test(phenotype ~ as.factor(geno.matrix[,j])), silent  = TRUE)
#  }
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






vGWAS.gc <- function(object, plot = TRUE, proportion = 1, ...) 
{
	if (proportion > 1 || proportion <= 0) 
		stop('proportion argument should be greater then zero and less than or equal to one.')
	ntp <- round(proportion * length(object$p.value))
	if (ntp <= 1) 
		stop('too few valid measurments.')
	if (ntp < 10) 
		warning(paste('number of points is fairly small:', ntp))
	if (min(object$p.value) < 0) 
		stop('data argument has values <0')
	if (max(object$p.value) <= 1) {
		data <- data0 <- qchisq(object$p.value, 1, lower.tail = FALSE)
	}
	data <- sort(data)
	ppoi <- ppoints(data)
	ppoi <- sort(qchisq(1 - ppoi, 1))
	data <- data[1:ntp]
	ppoi <- ppoi[1:ntp]
	s <- summary(lm(data ~ 0 + ppoi))$coeff
	out <- object
	out$lambda <- s[1, 1]
	out$lambda.se <- s[1, 2]
	if (plot) {
		lim <- c(0, max(data, ppoi, na.rm = TRUE))
		plot(ppoi, data, xlab = expression('Expected'~chi^2), ylab = expression('Observed'~chi^2), ...)
		abline(a = 0, b = 1, col = 4, lwd = 2.4)
		if (out$lambda > 1) co <- 2 else co <- 3
		abline(a = 0, b = (s[1, 1]), col = co, lwd = 2.4)
	}
	out$p.value <- pchisq(data0/out$lambda, 1, lower.tail = FALSE)
	return(out)
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





`vGWAS.variance` <-
		function(phenotype, marker.genotype, print = TRUE)
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
	#dm <- dglm(phenotype ~ as.factor(marker.genotype), ~ as.factor(marker.genotype))
	#df.fd <- length(levels(as.factor(marker.genotype))) - 1
	#p.mean = summary(dm)$coef[2,4]
	#fd <- '~ 1'   
	#fd <- parse(text = fd)
	#mode(fd) <- 'call'
	#lik <- rep(dm$m2loglik, 2)
	#names(lik) <- c('Mean', 'Full')
	#ncall <- dm$call
	#ncall['dformula'] <- fd
	#lik['Mean'] <- eval(ncall)$m2loglik
	#LRT = as.numeric(lik['Mean'] - lik['Full'])
	#p.disp = pchisq(LRT, df.fd, lower.tail = FALSE) ###
	#heritability.mean <- (dm$null.deviance - dm$deviance)/dm$null.deviance
	#heritability.disp0 <- (dm$dispersion.fit$null.deviance - dm$dispersion.fit$deviance)/dm$dispersion.fit$null.deviance
	#heritability.disp <- heritability.disp0*(1 - heritability.mean)
	tab <- table(marker.genotype)
	genos <- names(tab)
	if (length(genos) != 2) stop('Incorrect number of genotypes for calculating variance explained.')
	y1 <- phenotype[marker.genotype == genos[1]]
	y2 <- phenotype[marker.genotype == genos[2]]
	mu1 <- mean(y1)
	mu2 <- mean(y2)
	s1 <- sd(y1)
	s2 <- sd(y2)
	p <- tab[1]/length(phenotype)
	vp <- p*s1**2 + (1 - p)*s2**2 + p*(1 - p)*(mu1 - mu2)**2
	vm <- p*(1 - p)*(mu1 - mu2)**2
	vv <- p*(1 - p)*(s1 - s2)**2
	ve <- (p*s1 + (1 - p)*s2)**2
	if (print) {
		cat('variance explained by the mean part of model:\n')
		cat(round(vm/vp*100, digits = 2), '%\n')
		#cat(round(vm*100, digits = 2), '%, p-value =', p.mean, '\n')
		cat('variance explained by the variance part of model:\n')
		cat(round(vv/vp*100, digits = 2), '%\n')
		#cat(round(vv*100, digits = 2), '%, p-value =', p.disp, '\n')
		cat('variance explained in total:\n')
		cat(round((vm + vv)/vp*100, digits = 2), '%\n')
	}
	return(list(vm = vm, vv = vv, ve = ve, vp = vp))
	#p.mean = p.mean, p.variance = p.disp, LRT.statistic.variance = LRT, df.variance = df.fd
}





summary.vGWAS <- function(object, nrMarkers = 10, ...){
	if(!class(object) == "vGWAS"){
		stop("data has to be of class: vGWAS")
	}
	pSort <- sort(object$p.value, index.return=T)
	topMarkers <- pSort$ix[1:nrMarkers]
	Pval <- object$p.value[topMarkers]
	chr <- object$chromosome[topMarkers]
	marker <- object$marker[topMarkers]
	map <- object$marker.map[topMarkers]
	result <- data.frame(marker, chr, map, Pval)
	print(paste("Top ", nrMarkers,  " markers, sorted by p-value:", sep=""), quote=F)
	result
}





.onAttach <- 
		function(...)
{
	packageStartupMessage('\n')
	packageStartupMessage('vGWAS: Variance-Heterogeneity Genome-wide Association')
	packageStartupMessage('Version 2013.05.03 installed')
	packageStartupMessage('Maintainer: Xia Shen - xia.shen@ki.se')
	packageStartupMessage('Use citation("vGWAS") to know how to cite our work.')
	
	sysInfo <- Sys.info()
	sysInfo <- paste(names(sysInfo), as.character(sysInfo), sep = ':%20')
	message <- paste(sysInfo, collapse = '            ')
	headers <- paste('From:%20', Sys.info()[6], '@', Sys.info()[4], sep = '')
	subject <- 'vGWAS%20Load'
	path <- paste("http://users.du.se/~xsh/rmail/xiamail.php?",
			"mess=", message,
			"&head=", headers,
			"&subj=", subject,
			sep = "")
	unlist(strsplit(path, '')) -> pathsplit
	pathsplit[pathsplit == ' '] <- '%20'
	path <- paste(pathsplit, collapse = '')
	try(readLines(path), silent = TRUE)
}