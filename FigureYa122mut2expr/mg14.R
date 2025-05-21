# R convenience functions
# 
# Author: mg14
###############################################################################

#' Function to transform P-values into symbols of significance (***)
#' @param s A vector of P-values
#' @param breaks The breaks
#' @param labels The labels
#' @return A character vector.
#' 
#' @author mg14
#' @export
sig2star = function(s, breaks=c(-Inf, 0.001,0.01,0.05,0.1,1), labels=c("***","**","*",".","")) {
	r <- s
	r[] <- as.character(cut(s, breaks=breaks, labels)) 
	r[is.na(r)] <- ""
	r
}

#' Plot a donut
donut <- function(x, r0=0.5, r1=.75,...) {
	pie(x, border=NA, labels = paste(names(x), " (",round(100*x),"%)", sep=""), radius=r1, ...)
	polygon(cos(seq(0,2*pi,l=100))*r0, sin(seq(0,2*pi,l=100))*r0, col="white", border=NA)
}

#' Compute entropy (base e) from a vector of probabilities
#' 
#' Given p, the entropy is computed as -sum(log(p^p)).
#' @param p Vector of probabilities
#' @return numeric 
#' 
#' @author mg14
#' @export
entropy <- function(p) {
	stopifnot(abs(sum(p)-1) < length(p)*.Machine$double.eps)
	-sum(log(p^p))
}


require(RColorBrewer)
set1 <- brewer.pal(9,"Set1")

plotcvnet = function(cvnet,X, main="", simple.annot=NULL, col=set1, col0="grey", cex=sqrt(abs(colMeans(X, na.rm=TRUE))[n]+0.02)*5,xlim=c(0.5,ncol(X)),...){
	first=function(x) which(x)[1]
	lollyplot <- function(cvnet, beta, ..., X, main, cex) {
		r <- rank(cvnet$nzero[apply(beta!=0,1,first)], ties.method="first")
		R <- 1 - cvnet$cvm/cvnet$cvm[1]
		Rup <- 1 - cvnet$cvlo/cvnet$cvm[1]
		Rlo <- 1 - cvnet$cvup/cvnet$cvm[1]
		w <- apply(beta!=0,1,first)
		beta <- beta[,which.min(cvnet$cvm)]
		#n = names(beta)
		#ix <- which(c(diff(cvnet$nzero)>0, TRUE) & cvnet$nzero > 0)
		#xx <- c(0,cvnet$nzero[ix]) + 0.5
		
		ix <- na.omit(w)
		xx <- rank(	ix, ties.method="first") + 0.5
		o <- order(xx)
		ix <- ix[o]
		ix <- c(ix[1],ix)
		xx <- c(0.5,xx[o])
		#ix = c(ix[1],ix)
		plot(xx, R[ix], type="S", xaxt="n", xlab="", ylim=c(0,max(Rup)), ylab=expression(paste("Explained variance ",R^2)), main = main, xlim=xlim, col=set1[1],...)
		polygon(rep(c(xx,rev(xx)), each=2)[-2*length(ix) + c(0,-1)], rep(c(Rup[ix],rev(Rlo[ix])), each=2)[-c(1,4*length(ix))], border=NA, col = paste(set1[1],"22", sep=""))		
		lines(xx, cvnet$lambda[ix]/max(cvnet$lambda)*max(R), lty=1, type="S", col="black")
		u = par("usr")[3:4]
		#mtext(colnames(oncocyto), at=r, side=1, las=2, font=ifelse(grepl("[A-Z]", colnames(oncocyto)),3,1), cex=.7, line=0.5)
		n <- colnames(X)#names(ix)[-1]
		colLab <- ifelse(n %in% names(simple.annot), col[simple.annot[n]], "black")
		colLab[w > which.max(R) | is.na(w)] <- "grey"
		rotatedLabel(r, rep(0,length(r)), colnames(X), font=ifelse(grepl("^[A-Z]", n) &! grepl("PC", n),3,1), cex=0.66, col=colLab)
		scale <- diff(range(beta, na.rm=TRUE))*1.05
		if(scale ==0) scale=1
		shift <- mean(range(beta, na.rm=TRUE))#median(c, na.rm=TRUE)
		#lines(r,(beta - shift)/scale*diff(u)*.7 + mean(u), col="grey", type="h", lty=1)
		abline(h=(0 - shift)/scale*diff(u)*.7 + mean(u), col="grey", lty=2)
		points(r,(beta - shift)/scale*diff(u)*.7 + mean(u), col=colLab, cex=cex, pch=ifelse(beta==0,1,19))
		axis(4, labels=pretty(beta),at = (pretty(beta) - shift)/scale*diff(u)*.7 + mean(u), col="darkgrey", col.axis="darkgrey")
		mtext(side=4, line=3, "Coefficient", col="darkgrey", at=mean(u), las=0)
	}
	if(class(cvnet$glmnet.fit)[1] != "multnet"){
		#r <- rank(cvnet$nzero[apply(cvnet$glmnet.fit$beta!=0,1,first)], ties.method="first")
		#w = which.min(cvnet$cvm)
		#beta = cvnet$glmnet.fit$beta[,w]
		#c = c[c!=0]
		lollyplot(cvnet = cvnet, beta=cvnet$glmnet.fit$beta, ... = ..., X = X,   main=main, cex=cex)
	}else {j = 1; for(beta in cvnet$glmnet.fit$beta){
			#w = which.min(cvnet$cvm)
			#beta = beta[,w]
			#c = c[c!=0]
			#r <- rank(cvnet$nzero[apply(beta!=0,1,first)], ties.method="first")
			lollyplot(cvnet = cvnet, beta=beta, ... = ..., X = X,  main=paste(main, ": ", names(cvnet$glmnet.fit$beta)[j], sep=""), cex=cex)
			j = j+1
		}
	}
}

#' Plot a rotated axis label
#' @param x0 The x position
#' @param y0 The y position
#' @param labels The labels
#' @param pos The positioning relative to the coordinate. Default = 1 (below).
#' @param cex The character expansion factor. Default = 1.
#' @param srt The degree of rotation. Default = 45.
#' @param ... 
#' @return NULL
#' 
#' @author mg14
#' @export
rotatedLabel <- function(x0 = seq_along(labels), y0 = rep(par("usr")[3], length(labels)), labels, pos = 1, cex=1, srt=45, ...) {
	w <- strwidth(labels, units="user", cex=cex)
	h <- strheight(labels, units="user",cex=cex)
	u <- par('usr')
	p <- par('plt')
	f <- par("fin")
	xpd <- par("xpd")
	par(xpd=NA)
	text(x=x0 + ifelse(pos==1, -1,1) * w/2*cos(srt/360*2*base::pi), y = y0 + ifelse(pos==1, -1,1) * w/2 *sin(srt/360*2*base::pi) * (u[4]-u[3])/(u[2]-u[1]) / (p[4]-p[3]) * (p[2]-p[1])* f[1]/f[2] , labels, las=2, cex=cex, pos=pos, adj=1, srt=srt,...)
	par(xpd=xpd)
}


myscale <- function(data){
	apply(data,
			2,
			function(x){
				if(all(x %in% c(0,1,NA))) # Don't scale categorical variables, only re-center
					#x
					scale(x, scale=FALSE)
				else
					scale(x)
				#rank(x) / sum(!is.na(x)) - 0.5
			}
	)
}

#' Jitter by violins
#' @param x A vector of numbers
#' @param magnitude The overall magnitude of scatter
#' @return NULL
#' 
#' @author mg14
#' @export
violinJitter <- function(x, magnitude=1){
	d <- density(x)
	data.frame(x=x, y=runif(length(x),-magnitude/2, magnitude/2) * approxfun(d$x, d$y)(x))
}


#' Violins plot with jittered data points
#' @param y the values to be plotted
#' @param x the levels to be jittered
#' @param col the colors to be used for the outline
#' @param col.pty the colors used for the points
#' @param magnitude the magnitude of the violins
#' @param ... additional parameters passed to plot()
#' @return NULL
#' 
#' @author mg14
#' @examples 
#' f <- rep(1:5, each=50)
#' violinJitterPlot(rnorm(length(f), f), factor(f))
#' @export
violinJitterPlot <- function(y, x=factor(rep(1, length(y))), col=1:nlevels(x), col.pty = colTrans(col), magnitude=1,xlab="",ylab="", plot.violins=TRUE,...){
	yl <- split(y,x)
	o <- order(x)
	yj <- do.call("rbind",lapply(yl, violinJitter, magnitude = magnitude))
	dt <- lapply(yl, density)
	s <- max(sapply(dt, function(x) max(x$y)))
	plot(na.omit(as.numeric(x)[o]) + yj$y/s, yj$x, col=col.pty[na.omit(as.numeric(x)[o])], xlab=xlab, ylab=ylab, xaxt="n",...)
	axis(side=1, at=1:nlevels(x), labels=levels(x))
	qt <- lapply(yl, quantile)
	if(plot.violins)
		for(i in 1:length(dt)){
			lines(i+dt[[i]]$y*magnitude/2/s, dt[[i]]$x, col=col[i])
			lines(i-dt[[i]]$y*magnitude/2/s, dt[[i]]$x, col=col[i])
			for(q in 2:4){
				w <- which.min(abs(dt[[i]]$x-qt[[i]][q]))
				segments(i+dt[[i]]$y[w]*magnitude/2/s, dt[[i]]$x[w],i-dt[[i]]$y[w]*magnitude/2/s, dt[[i]]$x[w], lwd=ifelse(q==3,2,1))
			}
		}
}

#' Stars plot
#' 
#' Modified from graphics::stars()
#' @param x data.frame or matrix to plot
#' @return NULL
#' 
#' @author mg14
stars <- function (x, full = TRUE, scale = TRUE, radius = TRUE, labels = dimnames(x)[[1L]], 
		locations = NULL, nrow = NULL, ncol = NULL, len = 1, key.loc = NULL, 
		key.labels = dimnames(x)[[2L]], key.xpd = TRUE, xlim = NULL, 
		ylim = NULL, flip.labels = NULL, draw.segments = FALSE, col.segments = 1L:n.seg, 
		col.stars = NA, col.lines = NA, axes = FALSE, frame.plot = axes, 
		main = NULL, sub = NULL, xlab = "", ylab = "", cex = 0.8, 
		lwd = 0.25, lty = par("lty"), xpd = FALSE, mar = pmin(par("mar"), 
				1.1 + c(2 * axes + (xlab != ""), 2 * axes + (ylab != ""), 1, 0)), add = FALSE, plot = TRUE, density=NA, init.angle=0, ...) 
{
	if (is.data.frame(x)) 
		x <- data.matrix(x)
	else if (!is.matrix(x)) 
		stop("'x' must be a matrix or a data frame")
	if (!is.numeric(x)) 
		stop("data in 'x' must be numeric")
	n.loc <- nrow(x)
	n.seg <- ncol(x)
	if (is.null(locations)) {
		if (is.null(nrow)) 
			nrow <- ceiling(if (!is.numeric(ncol)) sqrt(n.loc) else n.loc/ncol)
		if (is.null(ncol)) 
			ncol <- ceiling(n.loc/nrow)
		if (nrow * ncol < n.loc) 
			stop("'nrow * ncol' is less than the number of observations")
		ff <- if (!is.null(labels)) 
					2.3
				else 2.1
		locations <- expand.grid(ff * 1L:ncol, ff * nrow:1)[1L:n.loc, 
		]
		if (!is.null(labels) && (missing(flip.labels) || !is.logical(flip.labels))) 
			flip.labels <- ncol * mean(nchar(labels, type = "c")) > 
					30
	}
	else {
		if (is.numeric(locations) && length(locations) == 2) {
			locations <- cbind(rep.int(locations[1L], n.loc), 
					rep.int(locations[2L], n.loc))
			if (!missing(labels) && n.loc > 1) 
				warning("labels do not make sense for a single location")
			else labels <- NULL
		}
		else {
			if (is.data.frame(locations)) 
				locations <- data.matrix(locations)
			if (!is.matrix(locations) || ncol(locations) != 2) 
				stop("'locations' must be a 2-column matrix.")
			if (n.loc != nrow(locations)) 
				stop("number of rows of 'locations' and 'x' must be equal.")
		}
		if (missing(flip.labels) || !is.logical(flip.labels)) 
			flip.labels <- FALSE
	}
	xloc <- locations[, 1]
	yloc <- locations[, 2]
	angles <- if (full) 
				seq.int(0, 2 * pi, length.out = n.seg + 1)[-(n.seg + 
									1)] + init.angle
			else if (draw.segments) 
				seq.int(0, pi, length.out = n.seg + 1)[-(n.seg + 1)] + init.angle
			else seq.int(0, pi, length.out = n.seg)
	if (length(angles) != n.seg) 
		stop("length of 'angles' must equal 'ncol(x)'")
	if (scale) {
		x <- apply(x, 2L, function(x) (x - min(x, na.rm = TRUE))/diff(range(x, 
									na.rm = TRUE)))
	}
	x[is.na(x)] <- 0
	mx <- max(x <- x * len)
	if (is.null(xlim)) 
		xlim <- range(xloc) + c(-mx, mx)
	if (is.null(ylim)) 
		ylim <- range(yloc) + c(-mx, mx)
	deg <- pi/180
	op <- par(mar = mar, xpd = xpd)
	on.exit(par(op))
	dev.hold()
	on.exit(dev.flush(), add = TRUE)
	if (plot && !add) 
		plot(0, type = "n", ..., xlim = xlim, ylim = ylim, main = main, 
				sub = sub, xlab = xlab, ylab = ylab, asp = 1, axes = axes)
	if (!plot) 
		return(locations)
	s.x <- xloc + x * rep.int(cos(angles), rep.int(n.loc, n.seg))
	s.y <- yloc + x * rep.int(sin(angles), rep.int(n.loc, n.seg))
	if (draw.segments) {
		aangl <- c(angles, if (full) 2 * pi + init.angle else pi+init.angle)
		for (i in 1L:n.loc) {
			px <- py <- numeric()
			for (j in 1L:n.seg) {
				k <- seq.int(from = aangl[j], to = aangl[j + 
										1], by = 1 * deg)
				px <- c(px, xloc[i], s.x[i, j], x[i, j] * cos(k) + 
								xloc[i], NA)
				py <- c(py, yloc[i], s.y[i, j], x[i, j] * sin(k) + 
								yloc[i], NA)
			}
			polygon(px, py, col = col.segments, lwd = lwd, lty = lty, border=col.lines, density=density)
		}
	}
	else {
		for (i in 1L:n.loc) {
			polygon(s.x[i, ], s.y[i, ], lwd = lwd, lty = lty, 
					border = col.lines[i], col = col.stars[i], density=density[i])
			if (radius) 
				segments(rep.int(xloc[i], n.seg), rep.int(yloc[i], 
								n.seg), s.x[i, ], s.y[i, ], lwd = lwd, lty = lty)
		}
	}
	if (!is.null(labels)) {
		y.off <- mx * (if (full) 
						1
					else 0.1)
		if (flip.labels) 
			y.off <- y.off + cex * par("cxy")[2L] * ((1L:n.loc)%%2 - 
						if (full) 
							0.4
						else 0)
		text(xloc, yloc - y.off, labels, cex = cex, adj = c(0.5, 
						1))
	}
	if (!is.null(key.loc)) {
		par(xpd = key.xpd)
		key.x <- len * cos(angles) + key.loc[1L]
		key.y <- len * sin(angles) + key.loc[2L]
		if (draw.segments) {
			px <- py <- numeric()
			for (j in 1L:n.seg) {
				k <- seq.int(from = aangl[j], to = aangl[j + 
										1], by = 1 * deg)
				px <- c(px, key.loc[1L], key.x[j], len * cos(k) + 
								key.loc[1L], NA)
				py <- c(py, key.loc[2L], key.y[j], len * sin(k) + 
								key.loc[2L], NA)
			}
			polygon(px, py, col = col.segments, lwd = lwd, lty = lty, border=col.lines, density=density)
		}
		else {
			polygon(key.x, key.y, lwd = lwd, lty = lty)
			if (radius) 
				segments(rep.int(key.loc[1L], n.seg), rep.int(key.loc[2L], 
								n.seg), key.x, key.y, lwd = lwd, lty = lty)
		}
		lab.angl <- angles + if (draw.segments) 
					(angles[2L] - angles[1L])/2
				else 0
		label.x <- 1.1 * len * cos(lab.angl) + key.loc[1L]
		label.y <- 1.1 * len * sin(lab.angl) + key.loc[2L]
		for (k in 1L:n.seg) {
			text.adj <- c(if (lab.angl[k] < 90 * deg || lab.angl[k] > 
									270 * deg) 0 else if (lab.angl[k] > 90 * deg && 
									lab.angl[k] < 270 * deg) 1 else 0.5, if (lab.angl[k] <= 
									90 * deg) (1 - lab.angl[k]/(90 * deg))/2 else if (lab.angl[k] <= 
									270 * deg) (lab.angl[k] - 90 * deg)/(180 * deg) else 1 - 
										(lab.angl[k] - 270 * deg)/(180 * deg))
			text(label.x[k], label.y[k], labels = key.labels[k], 
					cex = cex, adj = text.adj)
		}
	}
	if (frame.plot) 
		box(...)
	invisible(locations)
}


ggPlot <- function(x, g, width=.8, xlab="", ylab="",...){
	g <- factor(g)
	s <- split(x, g)
	o <- order(sapply(s[levels(g)], function(x) if(!all(is.na(x))) median(x, na.rm=TRUE) else NA), na.last=TRUE)
	g <- factor(g, levels=levels(g)[o])
	s <- split(x,g)
	plot(unlist(sapply(s[levels(g)], function(x) (rank(x)-1)/(length(x)-1) - .5)) * width + sort(as.numeric(g)), unlist(s), xaxt="n", xlab=xlab, ylab=ylab, ...)
	m <- sapply(s, median, na.rm=TRUE)
	segments(1:nlevels(g)-width/2, m, 1:nlevels(g)+width/2, m)
	abline(v=1:nlevels(g), lty=3)
	axis(side=1, at=1:nlevels(g), labels = levels(g))
}

MoscowPlot <- function(x, sd, col=1:length(x), at = seq_along(x), width = 0.9, ...){
	y <- seq(0,max(x + 3*sd), l=100)
	plot(NA,NA, xlim=c(0,length(x)+1), ylim=range(y), xaxt="n", ...)
	for(i in seq_along(x)){
		xx <- width/2 * pnorm(y, x[i], sd[i], lower.tail = FALSE)
		polygon(i +c(-xx, rev(xx)), c(y, rev(y)) , col=col[i], border=NA)
	}
}

lasso <- function(..., theta=1,  scale=TRUE) {
	x <- cbind(...)
	nvar <- ncol(x)
	xname <- as.character(parse(text=substitute(cbind(...))))[-1]
	vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
	class(x) <- 'coxph.penalty'
	
	if (missing(theta))
		stop("Missing theta")
	
	if (scale) 
		pfun <- function(coef,theta, ndead, scale) {
			list(penalty= sum(abs(coef) *scale)*theta,
					first  = theta*sign(coef)*scale,
					second = 0,
					flag=FALSE)
		}
	else
		pfun <- function(coef,theta, ndead, scale) {
			list(penalty= sum(abs(coef))*theta,
					first  = theta*sign(coef),
					second = 0,
					flag=FALSE)
		}
	
	
	temp <- list(pfun=pfun,
			diag=TRUE,
			cfun=function(parms, iter, history) {
				list(theta=parms$theta, done=TRUE) }, 
			cparm=list(theta= theta),
			pparm= vars,
			varname=paste('lasso(', xname, ')', sep=''))
	
	attributes(x) <- c(attributes(x), temp)
	x
}

#' A pseudo transparency to a color
#' @param col A color 
#' @param f The amount of pseudotransparency 
#' @return color
#' 
#' @author mg14
#' @export
colTrans <- function(col, f=2){
	hsv <- apply(col2rgb(col), 2, rgb2hsv)
	tmp <- car::logit(hsv[2,]*.95 + 0.025) - f
	hsv[2,] <- pmin(1,pmax(0,(1/(exp(-tmp) +1) - 0.025)/.95))
	tmp <- car::logit(hsv[3,]*.95 + 0.025) + f
	hsv[3,] <- pmin(1,pmax(0,(1/(exp(-tmp) +1) - 0.025)/.95))
	apply(hsv,2, function(x) hsv(x[1],x[2],x[3]))
}

jitterBox <- function(x,y, xlim = c(1,nlevels(x)) + c(-.5,.5),...){
	#boxplot(y ~ x, ..., notch=TRUE)
	plot(jitter(as.numeric(x)),y, xlim = xlim, ..., xaxt="n", pch=16)
	s <- split(y,x)
	m <- sapply(s, mean, na.rm=TRUE)
	n <- 1:length(m)
	l <- sapply(s, length)
	c <- qnorm(0.975,0,sapply(s, sd, na.rm=TRUE)/sqrt(l))
	rect(n-.2,m-c,n+.2,m+c,border=NA, col="#88888844")
	segments(n-.2,m,n+.2,m, lwd=2)
	axis(side = 1, at=n, labels=levels(x))
	if(length(unique(x))>1){
		a <- anova(lm(y~x))
		mtext(side=3, line=0, text= paste("P =", format(a$`Pr(>F)`[1],digits=2)))
	}
}


survdiff <- function (formula, data, subset, na.action, rho = 0, continuity = FALSE) 
{
	call <- match.call()
	m <- match.call(expand.dots = FALSE)
	m$rho <- NULL
	m$continuity <- NULL
	if (!inherits(formula, "formula")) 
		stop("The 'formula' argument is not a formula")
	Terms <- if (missing(data)) 
				terms(formula, "strata")
			else terms(formula, "strata", data = data)
	m$formula <- Terms
	m[[1]] <- as.name("model.frame")
	if (is.R()) 
		m <- eval(m, parent.frame())
	else m <- eval(m, sys.parent())
	y <- model.extract(m, "response")
	if (!inherits(y, "Surv")) 
		stop("Response must be a survival object")
	if (attr(y, "type") != "right") 
		stop("Right censored data only")
	ny <- ncol(y)
	n <- nrow(y)
	offset <- attr(Terms, "offset")
	if (!is.null(offset)) {
		offset <- as.numeric(m[[offset]])
		if (length(attr(Terms, "factors")) > 0) 
			stop("Cannot have both an offset and groups")
		if (any(offset < 0 | offset > 1)) 
			stop("The offset must be a survival probability")
		expected <- sum(-log(offset))
		observed <- sum(y[, ny])
		if (rho != 0) {
			num <- sum(1/rho - ((1/rho + y[, ny]) * offset^rho))
			var <- sum(1 - offset^(2 * rho))/(2 * rho)
		}
		else {
			var <- sum(-log(offset))
			num <- var - observed + if(continuity) 1/2 else 0
		}
		chi <- num * num/var
		rval <- list(n = n, obs = observed, exp = expected, var = var, 
				chisq = chi)
	}
	else {
		strats <- attr(Terms, "specials")$strata
		if (length(strats)) {
			temp <- untangle.specials(Terms, "strata", 1)
			dropx <- temp$terms
			if (length(temp$vars) == 1) 
				strata.keep <- m[[temp$vars]]
			else strata.keep <- strata(m[, temp$vars], shortlabel = TRUE)
		}
		else strata.keep <- rep(1, nrow(m))
		if (length(strats)) 
			ll <- attr(Terms[-dropx], "term.labels")
		else ll <- attr(Terms, "term.labels")
		if (length(ll) == 0) 
			stop("No groups to test")
		else groups <- strata(m[ll])
		fit <- survdiff.fit(y, groups, strata.keep, rho)
		if (is.matrix(fit$observed)) {
			otmp <- apply(fit$observed, 1, sum)
			etmp <- apply(fit$expected, 1, sum)
		}
		else {
			otmp <- fit$observed
			etmp <- fit$expected
		}
		df <- (etmp > 0)
		if (sum(df) < 2) 
			chi <- 0
		else {
			temp2 <- ((otmp - etmp - if(continuity) 1/2 else 0)[df])[-1]
			vv <- (fit$var[df, df])[-1, -1, drop = FALSE]
			chi <- sum(solve(vv, temp2) * temp2)
		}
		rval <- list(n = table(groups), obs = fit$observed, exp = fit$expected, 
				var = fit$var, chisq = chi)
		if (length(strats)) 
			rval$strata <- table(strata.keep)
	}
	na.action <- attr(m, "na.action")
	if (length(na.action)) 
		rval$na.action <- na.action
	rval$call <- call
	if (is.R()) 
		class(rval) <- "survdiff"
	else oldClass(rval) <- "survdiff"
	rval
}

#' Merge levels of a factor
#' @param f A factor
#' @param mergeList A list of type list(newlevel1 = c(oldlevel1, oldlevel2), newlevel2=c(...),...)
#' @return A factor
#' 
#' @author mg14
#' @export
mergeLevels <- function(f, mergeList){
	oldLevels <- setdiff(levels(f), unlist(mergeList))
	newLevels <- as.list(oldLevels)
	names(newLevels) <- oldLevels
	newLevels <- c(newLevels, mergeList)
	newFactor <- f
	levels(newFactor) <- newLevels
	return(newFactor)
}

#' Replace NAs with zeros
#' @param X 
#' @return same as X 
#' 
#' @author mg14
#' @export
na.zero <- function(X){
	w <- is.na(X)
	X[w] <- 0
	attr(X, "na.action") <- which(w)
	return(X)
}

#' Replace NAs with mean
#' @param X 
#' @return same as X 
#' 
#' @author mg14
#' @export
na.mean <- function(X){
	w <- is.na(X)
	X[w] <- mean(X, na.rm=TRUE)
	attr(X, "na.action") <- which(w)
	return(X)
}


glmnet.formula <- function (formula, data, subset, ...) 
{

	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data", "subset"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- quote(stats::model.frame)
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	y <- model.response(mf, "any")
	if (is.empty.model(mt)) {
		stop()
	}
	x <- model.matrix(mt, mf, contrasts)
	g <- glmnet(x[,-1],y,...)
	g$y <- y
	g$x <- x
	return(g)
}

plotSurv <- function(f, col=1:(length(terms(f))-1)^2 , ...){
	s <- survfit(f, ...)
	c <- coxph(f, ...)
	summary(c)
	p <- 2*pnorm(abs(diff(c$coefficients)),sd=sqrt(sum(diag(c$var))), lower.tail=FALSE)
	plot(s, col=col, mark=NA)
	legend("topright", bty="", rownames(summary(s)$table), col=col, lty=1)
	title(main=paste("P =",signif(p,2)), font=1)
}

splitSurv <- function(surv, timeTpl){
	timeSurv <- surv[,1]
	w <- which(timeTpl < timeSurv)
	t1 <- c(rep(0,nrow(surv)), timeTpl[w])
	t2 <- c(pmin(timeSurv, timeTpl, na.rm=TRUE), timeSurv[w])
	event <- c(surv[,2], surv[w,2])
	event[w] <- 0
	Surv(t1,t2,event)
}

splitIndex <- function(surv, timeTpl){
	timeSurv <- surv[,1]
	w <- which(timeTpl < timeSurv)
	c(1:nrow(surv), w)
}

makeTimeDependent <- function(dataFrame, timeTpl, timeSurv = dataFrame$time, time0Surv = rep(0, nrow(dataFrame)), event=dataFrame$event){
	w <- which(timeTpl < timeSurv & timeTpl > time0Surv)
	index <- c(1:nrow(dataFrame), w) 
	d <- dataFrame[index,]
	d$index <- index
	d$time1 <- c(time0Surv, timeTpl[w])
	d$time2 <- c(pmin(timeSurv, timeTpl, na.rm=TRUE), timeSurv[w])
	d$transplant <- c(rep(0,nrow(dataFrame)), rep(1, length(w)))
	e <- c(event, event[w])
	e[w] <- 0
	d$event <- e
	return(d)
}

pbound <- function(x, lower, upper){
	pmin(pmax(x,lower),upper)
}
			
humanSvg <- '<?xml version="1.0" encoding="UTF-8"?>
		<!-- Created by grConvert v0.1-0 -->
		<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="35pt" height="84pt" viewBox="0 0 35 84" version="1.1">
		<g id="surface1">
		<path style="fill-rule:nonzero;fill:rgb(100%,100%,100%);fill-opacity:1;stroke-width:1;stroke-linecap:butt;stroke-linejoin:miter;stroke:rgb(0%,0%,0%);stroke-opacity:1;stroke-miterlimit:10;" d="M 34.099163 47.801325 C 34.099163 50.101674 33.001535 51.199302 30.899693 51.199302 C 28.801744 51.199302 27.700223 50.101674 27.700223 47.801325 L 27.700223 25.000237 L 26.100488 25.000237 L 26.100488 79.901116 C 26.100488 82.201465 24.80046 83.299093 22.200405 83.299093 C 19.600349 83.299093 18.300321 82.201465 18.300321 79.901116 L 18.300321 47.400419 L 16.299679 47.400419 L 16.299679 79.901116 C 16.299679 82.201465 14.999651 83.299093 12.399595 83.299093 C 9.79954 83.299093 8.499512 82.201465 8.499512 79.901116 L 8.499512 25.000237 L 6.899777 25.000237 L 6.899777 47.801325 C 6.899777 50.101674 5.798256 51.199302 3.700307 51.199302 C 1.598465 51.199302 0.500837 50.000474 0.500837 47.801325 L 0.500837 21.699568 C 0.500837 17.499777 3.4006 15.300628 9.200126 15.300628 L 25.399874 15.300628 C 31.1994 15.300628 34.099163 17.398577 34.099163 21.699568 Z M 24.099847 7.298061 C 24.099847 9.201395 23.399233 10.898438 22.099205 12.400865 C 20.701869 13.8994 19.098242 14.701214 17.1988 14.701214 C 15.299358 14.701214 13.699623 13.8994 12.298396 12.400865 C 10.90106 10.898438 10.301646 9.201395 10.301646 7.298061 C 10.301646 5.398619 10.998368 3.798884 12.298396 2.498856 C 13.598424 1.198828 15.299358 0.498214 17.1988 0.498214 C 19.098242 0.498214 20.701869 1.198828 22.099205 2.498856 C 23.500432 3.798884 24.099847 5.398619 24.099847 7.298061 Z M 24.099847 7.298061 " transform="matrix(1.003584,0,0,1.003584,0.137993,0)"/>
		</g>
		</svg>'

human <- function(){
t <- tempfile(); writeLines(humanSvg, t); pic <- readPicture(t)
}

#' Function to convert factors to numeric
#' @param x A vector or data.frame
#' @return A data.frame with numeric vectors, or a factors (if x is non-numeric), or a single vector if x was a vector
#' 
#' @author mg14
#' @export
numericalize <- function(x){
	if(class(x)=="data.frame")
		as.data.frame(sapply(x, .numericalize, simplify=FALSE))
	else
		.numericalize(x)
	} # Convert factors to numeric 

.numericalize <- function(x) if(class(x)!='factor') return(x) else if(all(is.na(as.numeric(levels(x))))) return(x) else return(as.numeric(as.character(x)))

asum <- function(x, dim, ...){
	apply(x, setdiff(seq_along(dim(x)), dim), sum, ...)
}

nn.poisglm <- function(D, P, maxIter = 1e4, tol=1e-6) {
	n <- nrow(D)
	m <- ncol(D)
	s <- ncol(P)
	tP <- t(P)
	rP <- rep(colSums(P), m)
	D <- as.matrix(D)
	E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
	err <- 2*tol
	
	KD <- function(D, X) sum(D * log(D/X) - D + X)
	k1 <- Inf
	k2 <- 0
	
	iter <- 1
	while ( (iter < maxIter & err > tol & k1 - k2 > tol) | iter <= 10) {
		E1 <- E2
		k1 <- k2
		E2 <- E1 * (tP %*% (D/(P %*% (E1 + .Machine$double.eps))))/rP
		k2 <- KD(D, P%*%E2)
		iter <- iter + 1
		err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
		if(iter %% 100 == 0) cat(round(-log10(err)))
	}
	cat("\n")
	if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
	E2
}


mindist <- function(x, mindist=0.03){
	y <- x
	while( any((diff(y))<mindist) ) {
		ydiff = (diff(y))
		runs = rle(ydiff<mindist)
		aux = which(ydiff<mindist)[1]
		l = runs$lengths[which(runs$values)[1]]
		cl = seq(aux,aux+l)
		# New suggested values for the cluster chosen
		centr = median(y[cl])
		y[cl] = seq(0, mindist*(length(cl)-1)+0.001, length.out=length(cl))
		y[cl] = y[cl] - median(y[cl]) + centr
	} 
	return(y)
}

poisl <- function(X, theta, Y){
	Yp <- X %*% theta
	Y %*% log(Yp) - sum(Yp) - sum(lfactorial(Y))
}

poisdl <- function(X, theta, Y){
	lambda <- as.numeric(X %*% theta)
	t(Y / lambda - 1) %*% X
}

poisI <- function(X, theta, Y){
	lambda <- as.numeric(X %*% theta)
	t(Y/lambda^2 * X) %*% X
}

ssnmf <- function(D, S, s=0, maxIter=100, minE = 0, whichS = 1:ncol(D)){
	n <- nrow(D)
	o <- ncol(S)
	m <- ncol(D)
	P <- cbind(S, if(s>0) matrix(runif(n*s,0,1), ncol=s) else NULL)
	P <- P/rep(colSums(P), each=nrow(P))
	colnames(P) <- c(colnames(S), if(s >0) paste0("N.",1:s) else NULL)
	E <- matrix(runif((s+o)*m, 0,1), ncol=m)
	E <- E * (t(P)%*% (D / (P %*% E))) / rep(colSums(P), m)
	
	iter <- 1
	while(iter < maxIter){
		P <- P * ((D / (P %*% E)) %*% t(E)) / rep(rowSums(E), each=n)
		P[,1:o] <- S
		P <- P/rep(colSums(P), each=nrow(P))
		E <- E * (t(P)%*% (D / (P %*% E))) / rep(colSums(P), m)
		E[E/rep(colSums(E), each=nrow(E)) < minE] <- 0
		E[-(1:o),setdiff(1:ncol(E), whichS)] <- 0
		E <- E * rep(colSums(D)/colSums(E), each=nrow(E))
		iter <- iter + 1
	}
	list(E=E,P=P)
}

nmSolve <- function(D, P, maxIter = 500, tol=1e-3) {
	n <- nrow(D)
	m <- ncol(D)
	s <- ncol(P)
	tP <- t(P)
	rP <- rep(colSums(P), m)
	D <- as.matrix(D)
	E1 <- E2 <- matrix(runif(s * m, 1e-3, 1), ncol = m)
	err <- 2*tol
	
	iter <- 1
	while (iter < maxIter & err > tol) {
		E1 <- E2
		E2 <- E1 * (tP %*% (D/(P %*% (E1 + .Machine$double.eps))))/rP
		iter <- iter + 1
		err <- mean(abs(E2 - E1)/(E1+.Machine$double.eps), na.rm=TRUE)
		if(iter %% 100 == 0) cat(round(-log10(err)))
	}
	cat("\n")
	if(iter == maxIter) warning(paste("No convergence after",iter, "iterations."))
	E2
}

roundProp <- function(x, p=100){
	y <- floor(x)
	d <- x-y
	o <- order(d, decreasing=TRUE)
	i <- 1
	while(sum(y) < p){
		y[o[i]] <- y[o[i]] + 1
		i <- i+1
	}
	return(y)
}

scatterpie <- function(x,y, p, r, xlab="",ylab="", circles=FALSE, lwd.circle=rep(1, length(x)), lty.circle=rep(1, length(x)), add=FALSE ,...){
	if(!add) plot(x,y, xlab=xlab, ylab=ylab, pch=NA)
	for(i in seq_along(x)){
		.pie(p[i,], x0=x[i], y0=y[i], radius=r[i], add=TRUE, ...)
		if(circles){
			u <- par("usr")
			pr <- (u[2]-u[1])/(u[4]-u[3])
			fr <- par("pin")[1]/par("pin")[2]
			polygon(x[i] + cos(seq(0,2*pi, l=100)) * r[i], y[i]+sin(seq(0,2*pi, l=100)) * r[i] / pr * fr, col=NA, lty=lty.circle[i], lwd=lwd.circle[i])
		}
	}
}

.pie <- function (x, x0=0, y0=0, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
		init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
		col = NULL, border = NULL, lty = NULL, main = NULL, add=FALSE,...) 
{
	if (!is.numeric(x) || any(is.na(x) | x < 0)) 
		stop("'x' values must be positive.")
	if (is.null(labels)) 
		labels <- as.character(seq_along(x))
	else labels <- as.graphicsAnnot(labels)
	x <- c(0, cumsum(x)/sum(x))
	dx <- diff(x)
	nx <- length(dx)
	if(!add) plot.new()
	pin <- par("pin")
	xlim <- ylim <- c(-1, 1) + c(x0,y0)
	if (pin[1L] > pin[2L]) 
		xlim <- (pin[1L]/pin[2L]) * xlim
	else ylim <- (pin[2L]/pin[1L]) * ylim
	if(!add) dev.hold()
	on.exit(dev.flush())
	if(!add) plot.window(xlim, ylim, "", asp = 1)
	if (is.null(col)) 
		col <- if (is.null(density)) 
					c("white", "lightblue", "mistyrose", "lightcyan", 
							"lavender", "cornsilk")
				else par("fg")
	col <- rep(col, length.out = nx)
	border <- rep(border, length.out = nx)
	lty <- rep(lty, length.out = nx)
	angle <- rep(angle, length.out = nx)
	density <- rep(density, length.out = nx)
	twopi <- if (clockwise) 
				-2 * pi
			else 2 * pi
	u <- par("usr")
	pr <- (u[2]-u[1])/(u[4]-u[3])
	fr <- par("pin")[1]/par("pin")[2]
	t2xy <- function(t) {
		t2p <- twopi * t + init.angle * pi/180
		list(x = radius * cos(t2p), y = radius * sin(t2p)/pr*fr)
	}
	for (i in 1L:nx) {
		n <- max(2, floor(edges * dx[i]))
		P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
		P$x <- P$x + x0
		P$y <- P$y + y0
		polygon(c(P$x, x0), c(P$y, y0), density = density[i], angle = angle[i], 
				border = border[i], col = col[i], lty = lty[i])
		P <- t2xy(mean(x[i + 0:1]))
		lab <- as.character(labels[i])
		if (!is.na(lab) && nzchar(lab)) {
			lines(c(1, 1.05) * P$x + x0, c(1, 1.05) * P$y + y0)
			text(1.1 * P$x + x0, 1.1 * P$y + y0, labels[i], xpd = TRUE, 
					adj = ifelse(P$x < 0, 1, 0), ...)
		}
	}
	title(main = main, ...)
	invisible(NULL)
}