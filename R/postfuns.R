library("geiger")
library("redist")
library("ggplot2")
library(grid)
library(gridExtra)

##########################################################################
# map.summary
# This function produces a list of summary statistics for
# a redist.mcmc output object and a vector of party votes
# as well as total votes and population
##########################################################################

#algobj <- c1
#votedem <- map10@data$DEMOCRAT
#voterep <- map10@data$REPUBLICAN
#popvec <- map10@data$POP100
#algobj <- list(c1, c2, c3, c4, c5)

map.summary <- function(algobj, votedem, voterep, popvec, diagplots=TRUE) {

	if (class(algobj)=="list") {
		# diagnostics on multiple chains
		seglist <- vector("list", length(algobj))
		for(i in 1:length(algobj)) {
			eval(parse(text=paste0("c", i, "_seg <- redist.segcalc(algobj[[", i, "]], votedem, votedem+voterep)")))
			eval(parse(text=paste0("seglist[[", i, "]] <- mcmc(c", i, "_seg[!is.na(c", i, "_seg)])")))
		}
		seglist <- mcmc.list(seglist)
		gelman_diag <- gelman.diag(seglist)
		
		# print diagnostic plots
		if (diagplots==TRUE) {
			print(gelman.plot(seglist))
			print(traceplot(seglist))
		}
		
		# combine chains for subsequent analyses
		combined <- combinechains(algobj)
		algobj <- combined
	
	} else {
		gelman_diag <- NULL
	}
	
	ndists <- length(unique(algobj$partitions[,1]))
	votetot <- votedem + voterep

	############################# MAPWIDE #####################################
	################## republican dissimilarity index #########################
	dissim.rep <- redist.segcalc(algobj, voterep, popvec)

	############## republican seats and votes across simulations ##############
	repubseats <- rep(NA, ncol(algobj$partitions))
	voteshare <- matrix(NA, nrow=ndists, ncol=ncol(algobj$partitions))
	rownames(voteshare) <- sort(unique(algobj$partitions[,1]))
	for(i in 1:ncol(algobj$partitions)) {
	  voteshare[,i] <- tapply(voterep, algobj$partitions[,i], sum, na.rm=TRUE)/
										 tapply(votetot, algobj$partitions[,i], sum, na.rm=TRUE)
		repubseats[i] <- sum(voteshare[,i]>.5)
	}

	########################### partisan bias #################################
	statebase <- sum(voterep)/sum(votetot)
	equal <- .5 - statebase
	range <- .1
	inc <- seq(0, range, by = .01)
	bias <- unique(c(rev(equal - inc), equal + inc))
	algout <- algobj$partitions[,!is.na(dissim.rep)]
	repseats <- matrix(NA, nrow = ncol(algout), ncol = length(bias))
	for(j in 1:length(bias)){
		repseats[,j] <- pBias(votedem, voterep, algout, bias[j])
		# print(j)
	}

	## Make bias correspond to deviation from 50-50; convert to seat share
	bias <- seq(-1 * range, range, length = length(bias))
	repseats <- repseats/ndists
	repseats <- 1 - repseats

	## Plot step function
	xmin <- -1 * range
	xmax <- range

	## Calculate bias - get change points
	storebias <- rep(NA, nrow(repseats))

	for(j in 1:nrow(repseats)) {
		swing <- repseats[j,]
	
		## Get the bias
		mod <- lm(swing ~ bias)
		null <- predict(mod, data.frame(bias = bias))
	
		## Calculate the area
		gt0_area <- geiger:::.area.between.curves(bias[which(bias > 0)],
																						null[which(bias > 0)],
																						swing[which(bias > 0)],
																						xrange = c(-1,1))
		lt0_area <- geiger:::.area.between.curves(bias[which(bias <= 0)],
																						null[which(bias <= 0)],
																						swing[which(bias <= 0)],
																						xrange = c(-1,1))
		bias_area <- gt0_area + lt0_area
		storebias[j] <- bias_area
		#if(j %% 1000 == 0){
		#	print(j)
		#}
	} 

	##################### competitiveness (tam cho) ###########################
	xax <- 1 - algobj$distance_original[!is.na(dissim.rep)]

	t_p <- apply(
		algout, 2, function(x) {
			unq <- length(unique(x))
			return(sum(abs(tapply(votedem, x, sum) /
									 	(tapply(votedem, x, sum) +
									 	 tapply(voterep, x, sum)) - .5)) / unq)
		}
	)
	t_e <- apply(
		algout, 2, function(x){
			unq <- length(unique(x))
			n_d <- tapply(votedem, x, sum)
			n_r <- tapply(voterep, x, sum)
			b_r <- sum(n_d > n_r)
			return(abs(b_r / unq - .5))
		}
	)
	alpha <- 1
	beta <- 4/3
	f <- t_p * (1 + alpha * t_e) * beta ## which part to save from competitiveness??

	win <- ifelse(voteshare>.5, "R", "D")

	# diagnostics
	mh_acceptance <- mean(algobj$mhdecisions)

		# put it together 
	mapwide <- list(republican_dissimilarity=dissim.rep,
									republican_seats=repubseats,
									partisan_bias=storebias, 
									competitiveness=f)
	
	districtwide <- list(republican_voteshare=voteshare,
											 republican_win=win)
	
	if (is.null(gelman_diag)) {
		diag <- list(mh_acceptance=mh_acceptance)
	} else {
		diag <- list(mh_acceptance=mh_acceptance,
								 gelman_diag=gelman_diag)
	}
	
	final <- list(diagnostics=diag,
								mapwide=mapwide,
								districtwide=districtwide)
	return(final)
}


###########################################################################
# animate.sims
# This function takes a map object and redist.mcmc output and 
# creates an animated map over the simulations 
###########################################################################

animate.sims <- function(map, algobj) {
	frames <- 100
	ivec <- seq(1, ncol(algobj$partitions), by=ncol(algobj$partitions)/frames)
	for(i in ivec) {
	  pref <- deparse(substitute(algobj))
		if (i < 10) {name = paste(pref, '000',i,'plot.png',sep='')}
		if (i < 100 && i >= 10) {name = paste(pref, '00',i,'plot.png', sep='')}
		if (i >= 100) {name = paste(pref, '0', i,'plot.png', sep='')}
		png(name)
		parts <- algobj$partitions[,i]	
		plot(map, col=parts+2)
		dev.off()
	}
}

###########################################################################
# plot.compete
# Competitiveness (Tam Cho)
###########################################################################
	
plot.compete <- function(sumobj, algobj) {
  ndists <- length(unique(algobj$partitions[,1]))
	range <- .1
	bias_min <- - (range - -1 * range) * 1 / 2
	bias_max <- (range - -1 * range) * 1 / 2

	xax <- 1 - algobj$distance_original[!is.na(sumobj$mapwide$republican_dissimilarity)]
	fplot <- tapply(sumobj$mapwide$competitiveness, round(xax, 4), mean) 
	
	test <- (1 - -1) * (tapply(sumobj$mapwide$partisan_bias, round(xax, 4), mean) - bias_min) / 
	  (bias_max - bias_min) + -1
	x <- as.numeric(names(test)) 
	n <- table(round(xax, 4))
	colrep <- rep("black", length(test[x<.05]))
	
	x_cp <- as.numeric(names(fplot))
	
	plot(x_cp[x_cp<.05], 1-fplot[x_cp<.05], pch = 16,
	     main = "Competitiveness of Simulated Plans",
	     xlab = "% of Precincts Switched From Original District",
	     ylab = "Electoral Competitiveness towards Democrats",
	     xaxt = "n",
	     cex.lab = 1.2,
	     cex.axis = 1.2,
	     cex.main = 1.2,
	     col = colrep
	)
	axis(1, seq(0, 0.05, by = 0.01), c("0%", "1%", "2%", "3%", "4%", "5%"), cex.axis = 1.7)
	abline(h = 1-sumobj$mapwide$competitiveness[1])
}

###########################################################################
# plot.seats
# Republican seats 
###########################################################################

plot.seats <- function(sumobj, algobj) {
  print(qplot(sumobj$mapwide$republican_seats, geom="histogram", 
	            main="Distribution of Republican seats over simulations",
	            xlab="Republican seats",
	            ylab="Frequency"))
}

###########################################################################
# plot.voteshare
# Distribution of Democratic vote share for each district across simulations
###########################################################################	

plot.voteshare <- function(sumobj, algobj) {
  ndists <- length(unique(algobj$partitions[,1]))
  for(i in 1:ndists) {
	  title <- paste("District", i)
	  print(qplot(sumobj$districtwide$republican_voteshare[i,], geom="histogram",
	              main=title,
	              xlab="Republican vote share",
	              ylab="Frequency",
	              bins=length(sumobj$districtwide$republican_voteshare[i,])/1000) +
	  geom_vline(xintercept=mean(sumobj$districtwide$republican_voteshare[i,]), col="red"))
  }
}
	