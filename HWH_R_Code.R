# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

# R code for:

			## HARD-WORKING HELPERS CONTRIBUTE TO LONG BREEDER LIFESPANS IN COOPERATIVE BIRDS ##

# Philip A. Downing (philip.downing@biol.lu.se), Ashleigh S. Griffin, Charlie K. Cornwallis

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

## Packages ##

library(ape)
library(MCMCglmm)
library(metafor)
library(doBy)

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

## Functions ##

# lnRR
Calc.lnRR <- function(P1, P2){ES <- log(P1) - log(P2); return(ES)}

# var.lnRR
Calc.var.lnRR <- function(P1, N1, P2, N2){S2 <- ((1-P1)/(N1*P1)) + ((1-P2)/(N2*P2)); return(S2)}

# Hedges D
Calc.d <- function(CMean, CVAR, CN, EMean, EVAR, EN, adjusted=T){
                   sPooled <- sqrt(((EN - 1)*EVAR + (CN - 1)*CVAR) / (EN + CN - 2))  
                   d <- (EMean - CMean) / sPooled
                   H.d <- d * (1 - (3 / (4 * (EN + CN - 2) - 1)))
                   if(adjusted==T){return(H.d)}
                   if(adjusted==F){return(d)}}

# var D
Calc.SE.d <- function(CN, EN, d){
                      SE <- sqrt(( (CN + EN) / (CN * EN) ) + ( (d^2) / (2 * (CN + EN - 2) ) ) )
                      return(SE)}

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### DATA MANIPULATION ###

nLong <- read.csv(".../longData.csv")	# 23 species

# turn survival values into proportions
nLong$sxP2 <- nLong$sxP/100
nLong$sxG2 <- nLong$sxG/100

# create effect sizes
nLong$lnRR <- Calc.lnRR(nLong$sxG2, nLong$sxP2)		# since group first, + values means groups have higher sx
nLong$var.lnRR <- Calc.var.lnRR(nLong$sxG2, nLong$sxGN, nLong$sxP2, nLong$sxPN)

nLong$Dbreed <- Calc.d(nLong$brPFeed, nLong$brPVar, nLong$brPN, nLong$brGFeed, nLong$brGVar, nLong$brGN)
nLong$var.Dbreed <- Calc.SE.d(nLong$brPN, nLong$brGN, nLong$Dbreed)

nLong$Dhelp <- Calc.d(nLong$brGFeed, nLong$brGVar, nLong$brGN, nLong$hpFeed, nLong$hpVar, nLong$hpN)
nLong$var.Dhelp <- Calc.SE.d(nLong$brGN, nLong$hpN, nLong$Dhelp)

# create sample sizes for d
nLong$n.Dbreed <- nLong$brPN + nLong$brGN
nLong$n.Dhelp <- nLong$hpN

# to estimate study-level variance in metafor
nLong$ID <- paste(nLong$animal, 1:length(nLong$animal), sep="")			

# get rid of NAs in fixed effects #
nLongB <- nLong[-which(is.na(nLong$Dbreed)),]		# 32 effect sizes from 16 species
nLongC <- nLong[-which(is.na(nLong$Dhelp)),]		# 20 effect sizes from 10 species

## phylogenetic trees ##
trees <- read.nexus(".../BirdZilla1300.nex")
is.ultrametric(trees[[1]])

## trim the trees ##
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, nLong$animal)))]
myTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
myTrees <- lapply(myTrees, makeNodeLabel, method = "number")
nLong$animal[which((nLong$animal %in% myTrees[[1]]$tip.label) == FALSE)]
myTrees[[1]]$tip.label[which((myTrees[[1]]$tip.label %in% nLong$animal) == FALSE)]
# reduced to 23 species

## variance-covariance matrices ##
birdCor <- vcv.phylo(myTrees[[sample(1:1300, 1)]], cor=TRUE)
dropTipB <- myTrees[[1]]$tip.label[which(is.na(match(myTrees[[1]]$tip.label, nLongB$animal)))]
myTreesB <- lapply(myTrees, drop.tip, dropTipB, trim.internal=T)
myTreesB <- lapply(myTreesB, makeNodeLabel, method = "number")
birdCorB <- vcv.phylo(myTreesB[[sample(1:1300, 1)]], cor=TRUE)
dropTipC <- myTrees[[1]]$tip.label[which(is.na(match(myTrees[[1]]$tip.label, nLongC$animal)))]
myTreesC <- lapply(myTrees, drop.tip, dropTipC, trim.internal=T)
myTreesC <- lapply(myTreesC, makeNodeLabel, method = "number")
birdCorC <- vcv.phylo(myTreesC[[sample(1:1300, 1)]], cor=TRUE)
# note that because sample at random for the 'birdCor' VCV
# the metafor models will give slightly different parameter estimates to those reported below

## priors ##
prior <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=1,nu=0.002), G2=list(V=1,nu=0.002)))

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

										### PUBLICATION BIAS ###

## lnRR ##
lnRRModel <- rma(lnRR ~ 1, vi=var.lnRR, data=nLong)
summary(lnRRModel)   		# B = 0.10, lwr = 0.05, upr = 0.15 (Q = 157.6, p < 0.001)
forest.rma(lnRRModel, order="obs", showweights=TRUE)
trimfill(lnRRModel)			# 0 studies missing from the left-hand side
nLong$invSE <- 1/sqrt(nLong$var.lnRR)
summary(lm(scale(lnRR) ~ invSE, data=nLong))	# Egger's test: intercept = 0.36 (p = 0.16), slope = -0.03 (p = 0.09)
# within-study variance
wj <- 1/nLong$var.lnRR									# inverse measurement error / sampling error variance
k <- length(nLong$var.lnRR)								# number of studies
within <- sum(wj * (k-1)) / (sum(wj)^2 - sum(wj^2))		# 0.005
# between-study variance
0.0171 / (0.0171 + within)	# I2.between = 78 %


## d (breeders) ##
dModel.1 <- rma(Dbreed ~ 1, vi=var.Dbreed, data=nLongB)
summary(dModel.1)   		# B = -0.55, lwr = -0.72, upr = -0.37 (Q = 27.3, p = 0.66)
forest.rma(dModel.1, order="obs", showweights=TRUE)
trimfill(dModel.1)			# 3 studies missing from the right-hand side
nLongB$DbreedinvSE <- 1/sqrt(nLongB$var.Dbreed)
summary(lm(scale(Dbreed) ~ DbreedinvSE, data=nLongB))	# Egger's test: intercept = -0.83 (p = 0.11), slope = 0.45 (p = 0.09)
# within-study variance
wj <- 1/nLongB$var.Dbreed			# inverse measurement error / sampling error variance
k <- length(nLongB$var.Dbreed)		# number of studies
withinBreeders <- sum(wj * (k-1)) / (sum(wj)^2 - sum(wj^2))			# 0.261
# between-study variance
0.0000 / (0.0000 + withinBreeders)	# I2.between = 00 %


## d (helpers) ##
dModel.2 <- rma(Dhelp ~ 1, vi=var.Dhelp, data=nLongC)
summary(dModel.2)   		# B = 0.23, lwr = -0.16, upr = -0.62 (Q = 41.5, p = 0.002)
forest.rma(dModel.2, order="obs", showweights=TRUE)
trimfill(dModel.2)			# 2 studies missing from the left-hand side
nLongC$DhelpinvSE <- 1/sqrt(nLongC$var.Dhelp)
summary(lm(scale(Dhelp) ~ DhelpinvSE, data=nLongC))	# Egger's test: intercept = 0.86 (p = 0.40), slope = -0.52 (p = 0.38)
# within-study variance
wj <- 1/nLongC$var.Dhelp				# inverse measurement error / sampling error variance
k <- length(nLongC$var.Dhelp)		# number of studies
withinHelpers <- sum(wj * (k-1)) / (sum(wj)^2 - sum(wj^2))			# 0.350
# between-study variance
0.3950 / (0.3950 + withinHelpers)	# I2.between = 53 %


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

					### QUESTION A. DOES BREEDER SURVIVAL INCREASE WITH GROUP SIZE? ###


## POOLED BREEDER RESPONSE ##

# metafor #
modApooled.1 <- rma.mv(lnRR ~ 1, V=var.lnRR, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCor), data=nLong)
summary(modApooled.1)   				# B = 0.10, lwr = 0.04, upr = 0.16	(Q = 157.6, p < 0.001)
I2.between <- 0.0056 / (0.0056 + 0.0133 + 0.0000 + within) * 100		# 24 %
I2.sex <- 0.0133 / (0.0056 + 0.0133 + 0.0000 + within) * 100			# 56 %
I2.phylo <- 0.0000 / (0.0056 + 0.0133 + 0.0000 + within) * 100		# 0 %


# MCMCglmm #
MEV <- nLong$var.lnRR
INtree <- inverseA(myTrees[[1]], nodes="TIPS")
model.start <- MCMCglmm(lnRR ~ 1, random = ~ animal+commonName, data=nLong, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTrees[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(lnRR ~ 1, random = ~ animal+commonName, data=nLong, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modApooled.2.a <- model.start
modApooled.2.b <- model.start
modApooled.2.c <- model.start
save(modApooled.2.a, file=".../modApooled.2.a")
save(modApooled.2.b, file=".../modApooled.2.b")
save(modApooled.2.c, file=".../modApooled.2.c")
# chain convergence
hist(modApooled.2.a$Liab)
plot(modApooled.2.a$VCV)     		# animal and units close to 0
plot(modApooled.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modApooled.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modApooled.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modApooled.2.Sols <- mcmc.list(list(modApooled.2.a$Sol, modApooled.2.b$Sol, modApooled.2.c$Sol))
plot(modApooled.2.Sols)
gelman.diag(modApooled.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modApooled.2.a$VCV)    # units passed halfwidth
heidel.diag(modApooled.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modApooled.2.a)
posterior.mode(modApooled.2.a$Sol)		# intercept = 0.13
HPDinterval(modApooled.2.a$Sol)			# intercept = 0.01 to 0.24
# I^2 values
between <- posterior.mode(modApooled.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modApooled.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modApooled.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + within) * 100		# 36 %
I2.sex <- sex / (between + sex  + animal + within) * 100			# 9 %
I2.phylo <- animal / (between + sex  + animal + within) * 100		# 16 %



## SEX-SPECIFIC RESPONSES ##

# metafor #
modAsex.1 <- rma.mv(lnRR ~ sex-1, V=var.lnRR, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCor), data=nLong)
summary(modAsex.1)   				# female = 0.03 < 0.10 < 0.18; male = 0.03 < 0.10 < 0.17 (Q = 157.5, p < 0.001)
I2.between <- 0.0062 / (0.0062 + 0.0126 + 0.0000 + within) * 100		# 26 %
I2.sex <- 0.0129 / (0.0062 + 0.0126 + 0.0000 + within) * 100			# 55 %
I2.phylo <- 0.0000 / (0.0062 + 0.0126 + 0.0000 + within) * 100		# 0 %

# MCMCglmm #
MEV <- nLong$var.lnRR
INtree <- inverseA(myTrees[[1]], nodes="TIPS")
model.start <- MCMCglmm(lnRR ~ sex-1, random = ~ animal+commonName, data=nLong, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTrees[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(lnRR ~ sex-1, random = ~ animal+commonName, data=nLong, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modAsex.2.a <- model.start
modAsex.2.b <- model.start
modAsex.2.c <- model.start
save(modAsex.2.a, file=".../modAsex.2.a")
save(modAsex.2.b, file=".../modAsex.2.b")
save(modAsex.2.c, file=".../modAsex.2.c")
# chain convergence
hist(modAsex.2.a$Liab)
plot(modAsex.2.a$VCV)     	# animal and units close to 0
plot(modAsex.2.a$Sol)	   	# intercept estimate well mixed
autocorr(modAsex.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modAsex.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modAsex.2.Sols <- mcmc.list(list(modAsex.2.a$Sol, modAsex.2.b$Sol, modAsex.2.c$Sol))
plot(modAsex.2.Sols)
gelman.diag(modAsex.2.Sols)   	 	# upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modAsex.2.a$VCV)   		# units passed halfwidth
heidel.diag(modAsex.2.a$Sol)    		# intercept passed halfwidth
# model parameters
summary(modAsex.2.a)
posterior.mode(modAsex.2.a$Sol)		# female = 0.11; male = 0.11
HPDinterval(modAsex.2.a$Sol)			# female = 0.01 to 0.25; male = 0.01 to 0.25
# I^2 values
between <- posterior.mode(modAsex.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modAsex.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modAsex.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + within) * 100		# 44 %
I2.sex <- sex / (between + sex  + animal + within) * 100			# 11 %
I2.phylo <- animal / (between + sex  + animal + within) * 100		# 9 %


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

							### MEAN EFFECTS FOR D (BREEDERS) ###


## POOLED BREEDER RESPONSE ##

# metafor #
modXpooled.1 <- rma.mv(Dbreed ~ 1, V=var.Dbreed, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCorB), data=nLongB)
summary(modXpooled.1)   			# B = -0.62, lwr = -0.94, upr = -0.30 (Q = 27.3, p = 0.66)
I2.between <- 0.0000 / (0.0000+ 0.0000 + 0.0614 + withinBreeders) * 100		# 0 %
I2.sex <- 0.0000 / (0.0000 + 0.0000 + 0.0614 + withinBreeders) * 100			# 0 %
I2.phylo <- 0.0656 / (0.0000 + 0.0000 + 0.0614 + withinBreeders) * 100		# 20 %


# MCMCglmm #
MEV <- nLongB$var.Dbreed
INtree <- inverseA(myTreesB[[1]], nodes="TIPS")
model.start <- MCMCglmm(Dbreed ~ 1, random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTreesB[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(Dbreed ~ 1, random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modXpooled.2.a <- model.start
modXpooled.2.b <- model.start
modXpooled.2.c <- model.start
save(modXpooled.2.a, file=".../modXpooled.2.a")
save(modXpooled.2.b, file=".../modXpooled.2.b")
save(modXpooled.2.c, file=".../modXpooled.2.c")
# chain convergence
hist(modXpooled.2.a$Liab)
plot(modXpooled.2.a$VCV)     		# animal and units close to 0
plot(modXpooled.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modXpooled.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modXpooled.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modXpooled.2.Sols <- mcmc.list(list(modXpooled.2.a$Sol, modXpooled.2.b$Sol, modXpooled.2.c$Sol))
plot(modXpooled.2.Sols)
gelman.diag(modXpooled.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modXpooled.2.a$VCV)    # units passed halfwidth
heidel.diag(modXpooled.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modXpooled.2.a)
posterior.mode(modXpooled.2.a$Sol)		# intercept = -0.56
HPDinterval(modXpooled.2.a$Sol)			# Dbreed = -0.90 to -0.26
# I^2 values
between <- posterior.mode(modXpooled.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modXpooled.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modXpooled.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + withinBreeders) * 100		# 0 %
I2.sex <- sex / (between + sex + animal + withinBreeders) * 100				# 1 %
I2.phylo <- animal / (between + sex + animal + withinBreeders) * 100			# 0 %



## SEX-SPECIFIC RESPONSES ##

# metafor #
modXsex.1 <- rma.mv(Dbreed ~ sex-1, V=var.Dbreed, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCorB), data=nLongB)
summary(modXsex.1)   	# female = -0.99 < -0.63 < -0.26; male = -0.98 < -0.62 < -0.25 (Q = 27.3, p < 0.61)
I2.between <- 0.0000 / (0.0000 + 0.0000 + 0.0615 + withinBreeders) * 100		# 0 %
I2.sex <- 0.0000 / (0.0000 + 0.0000 + 0.0615 + withinBreeders) * 100			# 0 %
I2.phylo <- 0.0615 / (0.0000 + 0.0000 + 0.0615 + withinBreeders) * 100		# 20 %

MEV <- nLongB$var.Dbreed
INtree <- inverseA(myTreesB[[1]], nodes="TIPS")
model.start <- MCMCglmm(Dbreed ~ sex-1, random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTreesB[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(Dbreed ~ sex-1, random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modXsex.2.a <- model.start
modXsex.2.b <- model.start
modXsex.2.c <- model.start
save(modXsex.2.a, file=".../modXsex.2.a")
save(modXsex.2.b, file=".../modXsex.2.b")
save(modXsex.2.c, file=".../modXsex.2.c")
# chain convergence
hist(modXsex.2.a$Liab)
plot(modXsex.2.a$VCV)     		# animal and units close to 0
plot(modXsex.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modXsex.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modXsex.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modXsex.2.Sols <- mcmc.list(list(modXsex.2.a$Sol, modXsex.2.b$Sol, modXsex.2.c$Sol))
plot(modXsex.2.Sols)
gelman.diag(modXsex.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modXsex.2.a$VCV)    # units passed halfwidth
heidel.diag(modXsex.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modXsex.2.a)
posterior.mode(modXsex.2.a$Sol)		# F = -0.62, M = -0.61
HPDinterval(modXsex.2.a$Sol)			# F = -1.00 to -0.21, M = -0.96 to -0.24
# I^2 values
between <- posterior.mode(modXsex.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modXsex.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modXsex.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + withinBreeders) * 100		# 0 %
I2.sex <- sex / (between + sex + animal + withinBreeders) * 100				# 1 %
I2.phylo <- animal / (between + sex + animal + withinBreeders) * 100			# 0 %


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

			### QUESTION B. DO INCREASES IN BREEDER SURVIVAL DEPEND ON REDUCED INVESTMENT IN CARE? ###


## POOLED BREEDER RESPONSE ##

# metafor #
modBpooled.1 <- rma.mv(lnRR ~ Dbreed + log(n.Dbreed), V=var.lnRR, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCorB), data=nLongB)
summary(modBpooled.1)   				# slope = -0.10, lwr = -0.18, upr = -0.03 (Q = 91.2, p < 0.001)
I2.between <- 0.0026 / (0.0026 + 0.0084 + 0.0011 + within) * 100		# 16 %
I2.sex <- 0.0084 / (0.0026 + 0.0084 + 0.0011 + within) * 100			# 50 %
I2.phylo <- 0.0011 / (0.0026 + 0.0084 + 0.0011 + within) * 100		# 7 %


# MCMCglmm #
MEV <- nLongB$var.lnRR
INtree <- inverseA(myTreesB[[1]], nodes="TIPS")
model.start <- MCMCglmm(lnRR ~ Dbreed + log(n.Dbreed), random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTreesB[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(lnRR ~ Dbreed + log(n.Dbreed), random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modBpooled.2.a <- model.start
modBpooled.2.b <- model.start
modBpooled.2.c <- model.start
save(modBpooled.2.a, file=".../modBpooled.2.a")
save(modBpooled.2.b, file=".../modBpooled.2.b")
save(modBpooled.2.c, file=".../modBpooled.2.c")
# chain convergence
hist(modBpooled.2.a$Liab)
plot(modBpooled.2.a$VCV)     		# animal and units close to 0
plot(modBpooled.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modBpooled.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modBpooled.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modBpooled.2.Sols <- mcmc.list(list(modBpooled.2.a$Sol, modBpooled.2.b$Sol, modBpooled.2.c$Sol))
plot(modBpooled.2.Sols)
gelman.diag(modBpooled.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modBpooled.2.a$VCV)    # units passed halfwidth
heidel.diag(modBpooled.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modBpooled.2.a)
posterior.mode(modBpooled.2.a$Sol)		# intercept = 0.04, Dbreed = -0.10, n = 0.01
HPDinterval(modBpooled.2.a$Sol)			# Dbreed = -0.18 to -0.02
# I^2 values
between <- posterior.mode(modBpooled.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modBpooled.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modBpooled.2.a$VCV[,1])				# phylogenetic variance
I2.between <- between / (between + sex + animal + within) * 100		# 10 %
I2.sex <- sex / (between + sex + animal + within) * 100				# 36 %
I2.phylo <- animal / (between + sex + animal + within) * 100			# 13 %



## SEX-SPECIFIC RESPONSES ##

# metafor #
modBsex.1 <- rma.mv(lnRR ~ sex+sex:Dbreed-1 + log(n.Dbreed), V=var.lnRR, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCorB), data=nLongB)
summary(modBsex.1)   		# F slope = -0.27 < -0.14 < -0.02, M slope = -0.20 < -0.08 < 0.04 (Q = 82.4, p < 0.001)
I2.between <- 0.0024 / (0.0024 + 0.0087 + 0.0026 + within) * 100		# 13 %
I2.sex <- 0.0087 / (0.0024 + 0.0087 + 0.0026 + within) * 100			# 47 %
I2.phylo <- 0.0026 / (0.0024 + 0.0087 + 0.0026 + within) * 100		# 14 %

MEV <- nLongB$var.lnRR
INtree <- inverseA(myTreesB[[1]], nodes="TIPS")
model.start <- MCMCglmm(lnRR ~ sex+sex:Dbreed-1 + log(n.Dbreed), random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTreesB[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(lnRR ~ sex+sex:Dbreed-1 + log(n.Dbreed), random = ~ animal+commonName, data=nLongB, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modBsex.2.a <- model.start
modBsex.2.b <- model.start
modBsex.2.c <- model.start
save(modBsex.2.a, file=".../modBsex.2.a")
save(modBsex.2.b, file=".../modBsex.2.b")
save(modBsex.2.c, file=".../modBsex.2.c")
# chain convergence
hist(modBsex.2.a$Liab)
plot(modBsex.2.a$VCV)     		# animal and units close to 0
plot(modBsex.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modBsex.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modBsex.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modBsex.2.Sols <- mcmc.list(list(modBsex.2.a$Sol, modBsex.2.b$Sol, modBsex.2.c$Sol))
plot(modBsex.2.Sols)
gelman.diag(modBsex.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modBsex.2.a$VCV)    # units passed halfwidth
heidel.diag(modBsex.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modBsex.2.a)
posterior.mode(modBsex.2.a$Sol)		# F = -0.12, M = -0.07
HPDinterval(modBsex.2.a$Sol)			# F = -0.29 to -0.01, M = -0.20 to 0.06
# I^2 values
between <- posterior.mode(modBsex.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modBsex.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modBsex.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + within) * 100		# 35 %
I2.sex <- sex / (between + sex + animal + within) * 100				# 10 %
I2.phylo <- animal / (between + sex + animal + within) * 100			# 9 %


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

		### QUESTION C. DOES BREEDER INVESTMENT IN PARENTAL CARE DEPEND ON HOW MUCH CARE HELPERS PROVIDE? ###


## POOLED BREEDER RESPONSE ##

# metafor #
modCpooled.1 <- rma.mv(Dbreed ~ Dhelp + log(n.Dhelp), V=var.Dbreed, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCorC), data=nLongC)
summary(modCpooled.1)   				# slope = -0.48, lwr = -0.80, upr = -0.15 (Q = 4.0, p = 0.999)
I2.between <- 0.0000 / (0.0000 + 0.0000 + 0.0000 + withinBreeders) * 100		# 0 %
I2.sex <- 0.0000 / (0.0000 + 0.0000 + 0.0000 + withinBreeders) * 100			# 0 %
I2.phylo <- 0.0000 / (0.0000 + 0.0000 + 0.0000 + withinBreeders) * 100		# 0 %


# MCMCglmm #
MEV <- nLongC$var.lnRR
INtree <- inverseA(myTreesC[[1]], nodes="TIPS")
model.start <- MCMCglmm(Dbreed ~ Dhelp + log(n.Dhelp), random = ~ animal+commonName, data=nLongC, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTreesC[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(Dbreed ~ Dhelp + log(n.Dhelp), random = ~ animal+commonName, data=nLongC, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modCpooled.2.a <- model.start
modCpooled.2.b <- model.start
modCpooled.2.c <- model.start
save(modCpooled.2.a, file=".../modCpooled.2.a")
save(modCpooled.2.b, file=".../modCpooled.2.b")
save(modCpooled.2.c, file=".../modCpooled.2.c")
# chain convergence
hist(modCpooled.2.a$Liab)
plot(modCpooled.2.a$VCV)     		# animal and units close to 0
plot(modCpooled.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modCpooled.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modCpooled.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modCpooled.2.Sols <- mcmc.list(list(modCpooled.2.a$Sol, modCpooled.2.b$Sol, modCpooled.2.c$Sol))
plot(modCpooled.2.Sols)
gelman.diag(modCpooled.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modCpooled.2.a$VCV)    # units passed halfwidth
heidel.diag(modCpooled.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modCpooled.2.a)
posterior.mode(modCpooled.2.a$Sol)		# intercept = -0.25, Dhelp = -0.54, n = -0.14
HPDinterval(modCpooled.2.a$Sol)			# Dhelp = -0.73 to -0.31
# I^2 values
between <- posterior.mode(modCpooled.2.a$VCV[,4])		# between-study variance
sex <- posterior.mode(modCpooled.2.a$VCV[,2])			# between-sex variance
animal <- posterior.mode(modCpooled.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + withinBreeders) * 100		# 26 %
I2.sex <- sex / (between + sex + animal + withinBreeders) * 100				# 0 %
I2.phylo <- animal / (between + sex + animal + withinBreeders) * 100			# 0 %



## SEX-SPECIFIC RESPONSES ##

# metafor #
modCsex.1 <- rma.mv(Dbreed ~ sex+sex:Dhelp-1 + log(n.Dbreed), V=var.Dbreed, random = list(~ 1 | ID, ~ 1 | commonName, ~ 1 | animal), R = list(animal=birdCorC), data=nLongC)
summary(modCsex.1)   		# F slope = -1.00 < -0.54 < -0.08, M slope = -0.91 < -0.44 < 0.02 (Q = 4.4, p = 0.996)
I2.between <- 0.0000 / (0.0000 + 0.0000 + 0.0000 + withinBreeders) * 100		# 0 %
I2.sex <- 0.0000 / (0.0000 + 0.0000 + 0.0000 + withinBreeders) * 100			# 0 %
I2.phylo <- 0.0000 / (0.0000 + 0.0000 + 0.0000 + withinBreeders) * 100		# 0 %

MEV <- nLongC$var.Dbreed
INtree <- inverseA(myTreesC[[1]], nodes="TIPS")
model.start <- MCMCglmm(Dbreed ~ sex+sex:Dhelp-1 + log(n.Dbreed), random = ~ animal+commonName, data=nLongC, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=0, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), verbose=FALSE)
model <- model.start
for(i in 1:1300){
  INtree <- inverseA(myTreesC[[i]], nodes="TIPS")
  start <- list(Liab=model$Liab[1,], R=model$VCV[1,4], G=list(G1=model$VCV[1,1], G1=model$VCV[1,2]))
  model <- MCMCglmm(Dbreed ~ sex+sex:Dhelp-1 + log(n.Dbreed), random = ~ animal+commonName, data=nLongC, pl=TRUE, slice=TRUE, nitt=1000, thin=1, burnin=999, mev=MEV, prior=prior, ginverse=list(animal=INtree$Ainv), start=start, verbose=FALSE)
  if(i>300){
    model.start$VCV[i-300,] <- model$VCV[1,]
    model.start$Sol[i-300,] <- model$Sol[1,]
    model.start$Liab[i-300,] <- model$Liab[1,]
  }
}
modCsex.2.a <- model.start
modCsex.2.b <- model.start
modCsex.2.c <- model.start
save(modCsex.2.a, file=".../modCsex.2.a")
save(modCsex.2.b, file=".../modCsex.2.b")
save(modCsex.2.c, file=".../modCsex.2.c")
# chain convergence
hist(modCsex.2.a$Liab)
plot(modCsex.2.a$VCV)     		# animal and units close to 0
plot(modCsex.2.a$Sol)		   	# intercept estimate well mixed
autocorr(modCsex.2.a$VCV)   	# correlation between successive samples < 0.1 for all components
autocorr(modCsex.2.a$Sol)   	# correlation between successive samples < 0.1 for all components
modCsex.2.Sols <- mcmc.list(list(modCsex.2.a$Sol, modCsex.2.b$Sol, modCsex.2.c$Sol))
plot(modCsex.2.Sols)
gelman.diag(modCsex.2.Sols)     # upper CI = 1.00 for intercept suggesting convergence
heidel.diag(modCsex.2.a$VCV)    # units passed halfwidth
heidel.diag(modCsex.2.a$Sol)    # intercept passed halfwidth
# model parameters
summary(modCsex.2.a)
posterior.mode(modCsex.2.a$Sol)		# F = -0.68, M = -0.47
HPDinterval(modCsex.2.a$Sol)			# F = -1.02 to -0.04, M = -0.93 to 0.08
# I^2 values
between <- posterior.mode(modCsex.2.a$VCV[,4])			# between-study variance
sex <- posterior.mode(modCsex.2.a$VCV[,2])				# between-sex variance
animal <- posterior.mode(modCsex.2.a$VCV[,1])			# phylogenetic variance
I2.between <- between / (between + sex + animal + withinBreeders) * 100		# 1 %
I2.sex <- sex / (between + sex + animal + withinBreeders) * 100				# 1 %
I2.phylo <- animal / (between + sex + animal + withinBreeders) * 100			# 1 %


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

											### THE END ###

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #