###### Script Information ---------------------------------------------------------------------------------------------------------------------------
#
# Title: 03-Analysis.R
# Author: Mandy Moore
# 
# Date: September 2018
# Dataset: Wheat North
#
# Description: 
#   - MM altering the NVT MET code to suit rolling MET
#   - Extract data and models for frosted trials
#   - Conducting analysis of Wheat North MET including 
#     frosted trials to assess differences
#
# Input: Data from 02-METfix.R
# Children: 04-METout.R and 04METheatmaps

###### Data importation and package initialisation ----------------------------------------------------------------------------------------------------
rm(list = ls())

# Initialise packages
library(asreml)
library(ggplot2)

# Required for mac:
#dyn.load("/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/lib/server/libjvm.dylib")

# Required for model.fit:
source("1Functions/NVTfns5.6.R")  
setup.fn(computer='linux')

# Set-up MET details
met.name <- "WheatMainNorth"
date <- "2018-09-01"
date.out <- "2018-09-01"

# Import data from 02-METfixR
load("1RData/WheatMainNorth-2018-09-01-data4RMETcleanedready.RData")

###### Data Checks and Set-up -------------------------------------------------------------------------------------------------------------------------

# Set-up the mt for use in the models
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)

head(Rnvtdata)
str(Rnvtdata)
str(Rnvtmodels)

# Setup VarietyKeep and VarietyDrop if not done in 02-METfix.R, but should have been.
# Generate VarietyKeep and VarietyDrop, put VarietyDrop in sparse

# Rnvtdata$VarietyDrop <- Rnvtdata$VarietyKeep <- Rnvtdata$Variety
# Rnvtdata$VarietyKeep[substr(tolower(Rnvtdata$Variety), 1, 4)=="fill"] <- NA
# Rnvtdata$VarietyDrop[substr(tolower(Rnvtdata$Variety), 1, 4)!="fill"] <- NA
# Rnvtdata$VarietyKeep <- factor(Rnvtdata$VarietyKeep)
# Rnvtdata$VarietyDrop <- factor(Rnvtdata$VarietyDrop)
 
unique(Rnvtdata[order(Rnvtdata$Variety), c("Variety", "VarietyKeep", "VarietyDrop")])

#### Calculation of variety presence in environments -------------------------------------------------
enam <- as.character(levels(Rnvtdata$Experiment));enam # 134  experiments (3 removed by SAGI)
gnam <- as.character(levels(Rnvtdata$VarietyKeep));gnam 
length(gnam) # 129 varieties with no fillers
pres <- with(Rnvtdata[!is.na(Rnvtdata$yield),],table(VarietyKeep,Experiment))
dim(pres) # 129 x 134

ne <- length(unique(Rnvtdata$Experiment)); ne #134
ng <- length(unique(Rnvtdata$VarietyKeep[!is.na(Rnvtdata$VarietyKeep)])); ng # 129

# Degree of balance
nrow(unique(Rnvtdata[, c("Experiment", "VarietyKeep")]))/(ng*ne) #0.3074743

#### Model 1: DIAG ------------------------------------------------------------------------------------ 

# You will have to set up the model to have the correct covariates and lincol/linrow, 
# random and residual terms. 
# Use mt to help you !!!!!!!!!!!!!!!!!!

models.asr <- list()

mt
# ROUND 1: ----
asr.diag <- asreml(yield ~ Experiment + 
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs + 
                     at(Experiment, mt$lcol):lin(Range) + 
                     at(Experiment, mt$lrow):lin(Row),
                   sparse = ~ VarietyDrop,
                   random = ~ diag(Experiment):VarietyKeep +
                     at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   #ROUND1 residuals
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                     dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                     dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                     dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)

asr.diag <- update(asr.diag)
asr.diag$converge
asr.diag$loglik  #ROUND1: 12351.65
models.asr[["diagAR1"]] <- asr.diag #from Round 1 # loglik = 12351.65, 28 iters and 1.0 sec
models.asr[["diagAR1"]]$modelstats <- cbind.data.frame(niter = 28, mniter = 1.0)
models.asr[["diagAR1"]]$"%vaf" <- models.asr[["diagAR1"]]$"site%vaf" <- NA

# Run DIAG model with id:id residuals.
# ROUND 2: ----
asr.diag <- asreml(yield ~ Experiment + 
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs + 
                     at(Experiment, mt$lcol):lin(Range) + 
                     at(Experiment, mt$lrow):lin(Row),
                   sparse = ~ VarietyDrop,
                   random = ~ diag(Experiment):VarietyKeep +
                     at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   #ROUND 2 residuals
                   residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)

asr.diag <- update(asr.diag)
asr.diag$converge
asr.diag$loglik  #ROUND1: 12351.61, ROUND2: 10865.33
models.asr[["diagID"]] <- asr.diag #from Round 2 # loglik = 10865.33, 9 iterations @ 0.6 secs
models.asr[["diagID"]]$modelstats <- cbind.data.frame(niter = 9, mniter = 0.6)
models.asr[["diagID"]]$"%vaf" <- models.asr[["diagID"]]$"site%vaf" <- NA

#### Check Genetic Variances are all positive ---------------------------------------------------------------------------------------------------------

vc1 <- summary(asr.diag)$varcomp
unique(vc1[grep("Experiment:Variety", row.names(vc1)), "bound"]) #All good. none bound
subset(vc1, vc1$"%ch" > 5) # none

##### Check for any outliers ---------------------------------------------------------------------------
asr.diag <- models.asr[["diagAR1"]]
asr.diag <- update(asr.diag, aom = TRUE)
out1 <- cbind(Rnvtdata, stdCondRes = asr.diag$aom$R[,2])
subset(out1, abs(stdCondRes) > 4) # there is one, but marginal, let's leave it
subset(out1, abs(stdCondRes) > 4)[,c('Experiment','Range','Row','Rep','ColRep','RowRep','Variety','yield','stdCondRes')]
# Experiment Range Row Rep ColRep RowRep Variety    yield stdCondRes
# WMaA16BELL2     3  24   2      2      3  Qalbis 5.103571  -4.166661

# Investigate outlier
(tmp <- subset(out1, out1$Experiment=='WMaA16BELL2' & out1$Range=='3'))
hist(tmp$yield); tmp$yield
(tmp <- subset(out1, out1$Experiment=='WMaA16BELL2' & out1$Row=='24'))
hist(tmp$yield); tmp$yield
(tmp <- subset(out1, out1$Experiment=='WMaA16BELL2' & out1$Variety=='Qalbis'))
hist(tmp$yield); tmp$yield
tmp <- subset(out1, out1$Variety=='Qalbis')
hist(tmp$yield); tmp$yield

#### Model 2: rr1 + diag -----------------------------------------------------------------------------------------------

asreml.options(workspace = "1000mb", pworkspace="1000mb")

asr.rr1 <- asreml(yield ~ Experiment + 
                    at(Experiment, mt$animaldmg):animaldmg + 
                    at(Experiment, mt$est):est + 
                    at(Experiment, mt$earlygs):earlygs + 
                    at(Experiment, mt$lcol):lin(Range) + 
                    at(Experiment, mt$lrow):lin(Row),
       sparse =~ Experiment:VarietyDrop, 
       random = ~ rr(Experiment, 1):VarietyKeep + diag(Experiment):VarietyKeep + 
                  at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                  at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
       residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
       na.action = na.method(y='include', x='include'), data=Rnvtdata)
# Do the updates as a for loop so that you don't lose work if an update crashes mid-way. 
# Always save the asr object before running a 2nd loop!

save.image()
for(i in 1:10){
  print(i)
  # if(asr.rr1$converge) stop("converged - yaay")
  asr.rr1 <- update(asr.rr1)
} 

# Comment out 'if(asr.rr1$converge) stop("converged - yaay")' in order 
# to keep updating to remove the varcomp %ch warning

save.image()
asr.rr1$converge
asr.rr1$loglik #11800.08
subset(summary(asr.rr1)$varcomp, summary(asr.rr1)$varcomp$"%ch">1)
subset(summary(asr.rr1)$varcomp, bound=="B")
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 4.666087e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15DUAR4 2.761365e-08        NA      NA     B  NA
models.asr[["rr1ID"]] <- asr.rr1
models.asr[["rr1ID"]]$modelstats <- cbind.data.frame(niter = 38, mniter = 2.8)

#### Calculate %total vaf Model 2 ----------------------------------------------------------------------------------

vv <- asr.rr1$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=1)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/diag(Gmat)) # [1] 0.5135921     *** sum of the means
sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) # 0.5286612 *** mean of the sums
hist(diag(lam%*%t(lam))/diag(Gmat))
models.asr[["rr1ID"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(diag(Gmat))
models.asr[["rr1ID"]]$"site%vaf" <- diag(lam%*%t(lam))/diag(Gmat)

sort(models.asr[["rr1ID"]]$"site%vaf")

# MM experiment
diff <- list()
diff[[1]]<- c('FA1', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

#### Rolling MET sequence Model 3--------------------------------------------------------------------------------------- 
# The following code is for Rolling MET analysis
# vv <- as.data.frame(vv)
# vv$Component <- rownames(vv)
# names(vv)[1] <- "scratch" 
# write.csv(vv[,c(2,1)], "asr.rr1-varcomp-scratch.csv", row.names = F)
# 
# scratchvaf <- models.asr[["rr1ID"]]$"site%vaf"
# scratchvaf <- as.data.frame(scratchvaf)
# scratchvaf$Experiment <- rownames(scratchvaf)
# names(scratchvaf)[1] <- "scratch"
# write.csv(scratchvaf[,c(2,1)], "asr.rr1-vaf-scratch.csv", row.names = F)


#### Set-up Model 3: rr2 + diag -------------------------------------------------------------------------------
rr2.sv <- asreml(yield ~ Experiment + 
                   at(Experiment, mt$animaldmg):animaldmg + 
                   at(Experiment, mt$est):est + 
                   at(Experiment, mt$earlygs):earlygs + 
                   at(Experiment, mt$lcol):lin(Range) + 
                   at(Experiment, mt$lrow):lin(Row),
                 sparse =~ Experiment:VarietyDrop, 
                     random = ~ rr(Experiment, 2):VarietyKeep + diag(Experiment):VarietyKeep + 
                   at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                   at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                 residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                     na.action = na.method(y='include', x='include'), data=Rnvtdata, start.values = T)

#Get starting values set up
rr2.temp <- rr2.sv$vparameters.table

#Get values from previous model
rr1.gam <- matrix(summary(models.asr[['rr1ID']], vparameters=T)$vparameters[['rr(Experiment, 1):VarietyKeep']],ncol=2)
psi.gam <- matrix(summary(models.asr[['rr1ID']], vparameters=T)$vparameters[['Experiment:VarietyKeep']],ncol=1)
rr1.gam <- cbind(psi.gam,rr1.gam)
dimnames(rr1.gam) <- list(levels(Rnvtdata$Experiment),c('psi','rr.psi','lam1'))

#Put values from previous model in the starting values format.
rr2.temp$Value[grep('.*var', rr2.temp$Component)] <- 0 # set rr psi's to zero # 0.2*rr1.gam for fa1
rr2.temp$Value[grep('fa1', rr2.temp$Component)] <- rr1.gam[,"lam1"]
rr2.temp$Value[grep('fa2', rr2.temp$Component)] <- c(0,rep(0.02,ne-1))
rr2.temp$Value[grep('Experiment:VarietyKeep',rr2.temp$Component)] <- 0.8*rr1.gam[,"psi"] #for rr + rr2 # doesnt exist for fa1

tmp <- cbind.data.frame(Component = names(models.asr[["rr1ID"]]$vparameters), Value = models.asr[["rr1ID"]]$vparameters)
tmp$Value[abs(tmp$Value)<0.00001] <- 1e-4
rr2.SV <- merge(rr2.temp, tmp, by = "Component", all=T)

rr2.SV$Value <- rr2.SV$Value.y
unique(is.na(rr2.SV$Value.y)); rr2.SV[is.na(rr2.SV$Value.y)==TRUE,]
rr2.SV$Value[is.na(rr2.SV$Value.y)] <- rr2.SV$Value.x[is.na(rr2.SV$Value.y)]
rr2.SV <- rr2.SV[, c("Component", "Value", "Constraint")]
rr2.SV <- rr2.SV[order(rr2.SV$Component),]
rr2.SV$Component <- as.character(rr2.SV$Component) 
rownames(rr2.SV) <- NULL
rownames(rr2.SV) <- rr2.SV$Component
unique(rr2.SV$Component %in% rr2.temp$Component)
rr2.SV <- rr2.SV[rr2.temp$Component,]
rownames(rr2.SV) <- NULL

#Set any small psi's, probably bound terms to something not so small
rr2.SV$Value[abs(rr2.SV$Value)<0.00001 & rr2.SV$Constraint != "F"] <- 1e-4
save.image()

#### Run Model 3: rr2 + diag -------------------------------------------------------------------------------
## asr.rr1$loglik was #11800.08
## first iteration of rr2+diag was 12274.70
## so all OK
## 6 sec per iteration, 141 iterations before converged
## this is fine

save(list = ls(), file = "1RData/RMETanalysis.RData")
asr.rr2 <- asreml(yield ~ Experiment + 
                    at(Experiment, mt$animaldmg):animaldmg + 
                    at(Experiment, mt$est):est + 
                    at(Experiment, mt$earlygs):earlygs + 
                    at(Experiment, mt$lcol):lin(Range) + 
                    at(Experiment, mt$lrow):lin(Row),
                  sparse =~ Experiment:VarietyDrop, 
                  random = ~ rr(Experiment, 2):VarietyKeep + diag(Experiment):VarietyKeep + 
                    at(Experiment, mt$crep):Rep + at(Experiment, mt$rrep):RowRep + 
                    at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                  residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                  na.action = na.method(y='include', x='include'), data=Rnvtdata, G.param = rr2.SV, R.param = rr2.SV)
#Do this as a for loop so that you don't lose work if an update crashes mid-way. Always save the asr object before running a 2nd loop!
save.image()

for(i in 1:10){#x2
  print(i)
  if(asr.rr2$converge) stop("converged - yaay")
  asr.rr2 <- update(asr.rr2)
}

#1 update, finished on 10th iteration, but then lots of updates to make %ch < 1
asr.rr2$loglik #12274.70  (after updating to remove %ch.)
subset(summary(asr.rr2)$varcomp, summary(asr.rr2)$varcomp$"%ch">1)
subset(summary(asr.rr2)$varcomp, bound == "B")
# Experiment:VarietyKeep!Experiment_WMaA13BILO4 1.382391e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14EMER4 4.974823e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15DUAR4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15ROMA4 2.669558e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15YELA4 3.555585e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA17LUND4 5.192497e-07        NA      NA     B   0
asr.rr2 <- update(asr.rr2) #+30 updates to get % ch smaller than 1% was 

models.asr[["rr2ID"]] <- asr.rr2
models.asr[["rr2ID"]]$modelstats <- cbind.data.frame(niter = 141, mniter = 6.0)

#### Calculate %total vaf Model 3 ----------------------------------------------------------------------------------

vv <- asr.rr2$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=2)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
gvar <- diag(Gmat) 
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/gvar) # [1] 0.6404558    *** sum of the means
sum(diag(lam %*% t(lam)))/sum(gvar) # 0.6362182  *** mean of the sums
hist(diag(lam%*%t(lam))/gvar)

diff[[2]]<- c('FA2', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

models.asr[["rr2ID"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(gvar)
models.asr[["rr2ID"]]$"site%vaf" <- diag(lam%*%t(lam))/gvar

sort(models.asr[["rr2ID"]]$"site%vaf")

#### Rolling MET sequence Model 3 --------------------------------------------------------------------------------------- 
# The following code is for Rolling MET analysis

# vv <- as.data.frame(vv)
# vv$Component <- rownames(vv)
# names(vv)[1] <- "scratch" 
# write.csv(vv[,c(2,1)], "asr.rr2-varcomp-scratch.csv", row.names = F)
# 
# scratchvaf <- models.asr[["rr2ID"]]$"site%vaf"
# scratchvaf <- as.data.frame(scratchvaf)
# scratchvaf$Experiment <- rownames(scratchvaf)
# names(scratchvaf)[1] <- "scratch"
# write.csv(scratchvaf[,c(2,1)], "asr.rr2-vaf-scratch.csv", row.names = F)

#### Set-up Model 4: rr3 + diag -------------------------------------------------------------------------------
rr3.sv <- asreml(yield ~ Experiment + 
                   at(Experiment, mt$animaldmg):animaldmg + 
                   at(Experiment, mt$est):est + 
                   at(Experiment, mt$earlygs):earlygs + 
                   at(Experiment, mt$lcol):lin(Range) + 
                   at(Experiment, mt$lrow):lin(Row),
                 random = ~ rr(Experiment, 3):VarietyKeep + diag(Experiment):VarietyKeep + 
                   at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                   at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                 residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                 na.action = na.method(y='include', x='include'), data=Rnvtdata, start.values = T)
#Get starting values set up
rr3.temp <- rr3.sv$vparameters.table

#Get values from previous model
rr2.gam <- matrix(summary(models.asr[['rr2ID']], vparameters=T)$vparameters[['rr(Experiment, 2):VarietyKeep']],ncol=3)
psi.gam <- matrix(summary(models.asr[['rr2ID']], vparameters=T)$vparameters[['Experiment:VarietyKeep']],ncol=1)
rr2.gam <- cbind(psi.gam,rr2.gam)
dimnames(rr2.gam) <- list(levels(Rnvtdata$Experiment),c('psi','rr.psi','lam1', 'lam2'))


#Put values from previous model in the starting values format.
rr3.temp$Value[grep('.*var', rr3.temp$Component)] <- 0 # set rr psi's to zero # 0.2*rr2.gam for fa1
rr3.temp$Value[grep('fa1', rr3.temp$Component)] <- rr2.gam[,"lam1"]
rr3.temp$Value[grep('fa2', rr3.temp$Component)] <- rr2.gam[,"lam2"]
rr3.temp$Value[grep('fa3', rr3.temp$Component)] <- c(0, 0, rep(0.02,ne-2))
rr3.temp$Value[grep('Experiment:VarietyKeep',rr3.temp$Component)] <- 0.8*rr2.gam[,"psi"] #for rr + rr3 # doesnt exist for fa1

tmp <- cbind.data.frame(Component = names(models.asr[["rr2ID"]]$vparameters), Value = models.asr[["rr2ID"]]$vparameters)
tmp$Value[abs(tmp$Value)<0.00001] <- 1e-4
rr3.SV <- merge(rr3.temp, tmp, by = "Component", all=T)

rr3.SV$Value <- rr3.SV$Value.y
unique(is.na(rr3.SV$Value.y)); rr3.SV[is.na(rr3.SV$Value.y)==TRUE,]
rr3.SV$Value[is.na(rr3.SV$Value.y)] <- rr3.SV$Value.x[is.na(rr3.SV$Value.y)]
rr3.SV <- rr3.SV[, c("Component", "Value", "Constraint")]
rr3.SV <- rr3.SV[order(rr3.SV$Component),]
rr3.SV$Component <- as.character(rr3.SV$Component) 
rownames(rr3.SV) <- NULL
rownames(rr3.SV) <- rr3.SV$Component
unique(rr3.SV$Component %in% rr3.temp$Component)
rr3.SV <- rr3.SV[rr3.temp$Component,]
rownames(rr3.SV) <- NULL

#Set any small psi's, probably bound terms to something not so small
rr3.SV$Value[abs(rr3.SV$Value)<0.00001 & rr3.SV$Constraint != "F"] <- 1e-4
save.image()

#### Run Model 4: rr3 + diag -------------------------------------------------------------------------------
## asr.rr2$loglik was 12274.70
## first iteration of rr3+diag was 12614.17
## so all OK
## 11.5 sec per iteration, 59 iterations before converged
## this is fine

save(list = ls(), file = "1RData/RMETanalysis.RData")
asr.rr3 <-asreml(yield ~ Experiment + 
                   at(Experiment, mt$animaldmg):animaldmg + 
                   at(Experiment, mt$est):est + 
                   at(Experiment, mt$earlygs):earlygs + 
                   at(Experiment, mt$lcol):lin(Range) + 
                   at(Experiment, mt$lrow):lin(Row),
                  sparse =~ Experiment:VarietyDrop, 
                  random = ~ rr(Experiment, 3):VarietyKeep + diag(Experiment):VarietyKeep + 
                   at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                   at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                 residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                  na.action = na.method(y='include', x='include'), data=Rnvtdata, G.param = rr3.SV, R.param = rr3.SV)
#Do this as a for loop so that you don't lose work if an update crashes mid-way. 
#Always save the asr object before running a 2nd loop!
save.image()
for(i in 1:10){#x1
  print(i)
  #if(asr.rr3$converge) stop("converged - yaay")
  asr.rr3 <- update(asr.rr3)
}
#3 updates, finished on 7th iteration up date, but then lots of updates to make %ch < 1
asr.rr3$loglik #12614.17
subset(summary(asr.rr3)$varcomp, summary(asr.rr3)$varcomp$"%ch">1)
subset(summary(asr.rr3)$varcomp, bound == "B")
# Experiment:VarietyKeep!Experiment_WMaA13BILO4 1.529928e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14EMER4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15DUAR4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15ROMA4 6.777597e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15YELA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17LUND4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17MUNG4 5.241470e-07        NA      NA     B   0
asr.rr3 <- update(asr.rr3) 

models.asr[["rr3ID"]] <- asr.rr3
models.asr[["rr3ID"]]$modelstats <- cbind.data.frame(niter = 59, mniter = 11.5)

#### Calculate %total vaf Model 4 ----------------------------------------------------------------------------------

vv <- asr.rr3$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=3)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
gvar <- diag(Gmat) 
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/diag(Gmat)) # [1] 0.7226324    *** sum of the means
sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) #0.7224958  *** mean of the sums

diff[[3]]<- c('FA3', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

hist(diag(lam%*%t(lam))/diag(Gmat))
models.asr[["rr3ID"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(diag(Gmat))
models.asr[["rr3ID"]]$"site%vaf" <- diag(lam%*%t(lam))/diag(Gmat)

sort(models.asr[["rr3ID"]]$"site%vaf")

#### Rolling MET sequence Model 3 --------------------------------------------------------------------------------------- 
# The following code is for Rolling MET analysis

# vv <- as.data.frame(vv)
# vv$Component <- rownames(vv)
# names(vv)[1] <- "scratch" 
# write.csv(vv[,c(2,1)], "asr.rr3-varcomp-scratch.csv", row.names = F)
# 
# scratchvaf <- models.asr[["rr3ID"]]$"site%vaf"
# scratchvaf <- as.data.frame(scratchvaf)
# scratchvaf$Experiment <- rownames(scratchvaf)
# names(scratchvaf)[1] <- "scratch"
# write.csv(scratchvaf[,c(2,1)], "asr.rr3-vaf-scratch.csv", row.names = F)

#### Set-up Model 5: rr4 + diag -------------------------------------------------------------------------------
save(list = ls(), file = "1RData/RMETanalysis.RData")
rr4.sv <- asreml(yield ~ Experiment + 
                   at(Experiment, mt$animaldmg):animaldmg + 
                   at(Experiment, mt$est):est + 
                   at(Experiment, mt$earlygs):earlygs + 
                   at(Experiment, mt$lcol):lin(Range) + 
                   at(Experiment, mt$lrow):lin(Row),
                 sparse =~ Experiment:VarietyDrop, 
                 random = ~ rr(Experiment, 4):VarietyKeep + diag(Experiment):VarietyKeep + 
                   at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                   at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                 residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                 na.action = na.method(y='include', x='include'), data=Rnvtdata, start.values = T)
#Get starting values set up
rr4.temp <- rr4.sv$vparameters.table

#Get values from previous model
rr3.gam <- matrix(summary(models.asr[['rr3ID']], vparameters=T)$vparameters[['rr(Experiment, 3):VarietyKeep']],ncol=4)
psi.gam <- matrix(summary(models.asr[['rr3ID']], vparameters=T)$vparameters[['Experiment:VarietyKeep']],ncol=1)
rr3.gam <- cbind(psi.gam,rr3.gam)
dimnames(rr3.gam) <- list(levels(Rnvtdata$Experiment),c('psi','rr.psi','lam1', 'lam2', 'lam3'))

#Put values from previous model in the starting values format.
rr4.temp$Value[grep('.*var', rr4.temp$Component)] <- 0 # set rr psi's to zero # 0.2*rr3.gam for fa1
rr4.temp$Value[grep('fa1', rr4.temp$Component)] <- rr3.gam[,"lam1"]
rr4.temp$Value[grep('fa2', rr4.temp$Component)] <- rr3.gam[,"lam2"]
rr4.temp$Value[grep('fa3', rr4.temp$Component)] <- rr3.gam[,"lam3"]
rr4.temp$Value[grep('fa4', rr4.temp$Component)] <- c(0,0, 0,rep(0.02,ne-3))
rr4.temp$Value[grep('Experiment:VarietyKeep',rr4.temp$Component)] <- 0.8*rr3.gam[,"psi"] #for rr + rr4 # doesnt exist for fa1

tmp <- cbind.data.frame(Component = names(models.asr[["rr3ID"]]$vparameters), Value = models.asr[["rr3ID"]]$vparameters)
tmp$Value[abs(tmp$Value)<0.00001] <- 1e-4
rr4.SV <- merge(rr4.temp, tmp, by = "Component", all=T)

rr4.SV$Value <- rr4.SV$Value.y
unique(is.na(rr4.SV$Value.y)); rr4.SV[is.na(rr4.SV$Value.y)==TRUE,]
rr4.SV$Value[is.na(rr4.SV$Value.y)] <- rr4.SV$Value.x[is.na(rr4.SV$Value.y)]
rr4.SV <- rr4.SV[, c("Component", "Value", "Constraint")]
rr4.SV <- rr4.SV[order(rr4.SV$Component),]
rr4.SV$Component <- as.character(rr4.SV$Component) 
rownames(rr4.SV) <- NULL
rownames(rr4.SV) <- rr4.SV$Component
unique(rr4.SV$Component %in% rr4.temp$Component)
rr4.SV <- rr4.SV[rr4.temp$Component,]
rownames(rr4.SV) <- NULL

#Set any small psi's, probably bound terms to something not so small
rr4.SV$Value[abs(rr4.SV$Value)<0.00001 & rr4.SV$Constraint != "F"] <- 1e-4
save.image()

#### Run Model 5: rr4 + diag -------------------------------------------------------------------------------
## asr.rr3$loglik was 12614.26
## first iteration of rr4+diag was 12614.26 (12822.6 -converged)
## so all OK
## 2.2 sec per iteration, 51 iterations before converged
## this is fine
#######################
## asr.rr4$loglik
##  first iteration of rr4+diag was 12614.26
## so all OK
## 2.2 sec per iteration, 51 iterations before converged
## this is fine
## 
save(list = ls(), file = "1RData/RMETanalysis.RData")
asr.rr4 <- asreml(yield ~ Experiment + 
                    at(Experiment, mt$animaldmg):animaldmg + 
                    at(Experiment, mt$est):est + 
                    at(Experiment, mt$earlygs):earlygs + 
                    at(Experiment, mt$lcol):lin(Range) + 
                    at(Experiment, mt$lrow):lin(Row),
                  sparse =~ Experiment:VarietyDrop, 
                  random = ~ rr(Experiment, 4):VarietyKeep + diag(Experiment):VarietyKeep + 
                    at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                    at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                  residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                  na.action = na.method(y='include', x='include'), data=Rnvtdata, G.param = rr4.SV, R.param = rr4.SV)
#Do this as a for loop so that you don't lose work if an update crashes mid-way. Always save the asr object before running a 2nd loop!
save.image()
for(i in 1:10){#x1
  print(i)
  # if(asr.rr4$converge) stop("converged - yaay")
  asr.rr4 <- update(asr.rr4)
}
#2 loops, finished on 9th rep pf the 2nd loop. (loops have 13 iterations) but then lots of updates to make %ch < 1
asr.rr4$loglik #loop1: 12795.79, loop2: 12822.6(converged), 12822.82(after remove 1% change)
asr.rr4$converge
subset(summary(asr.rr4)$varcomp, summary(asr.rr4)$varcomp$"%ch">5)
subset(summary(asr.rr4)$varcomp, bound == "B")
# Experiment:VarietyKeep!Experiment_WMaA13BILO4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13DUAR4 5.357432e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14EMER4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14ROMA4 1.008727e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14SURA4 4.734106e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15DUAR4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15ROMA4 2.262035e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15YELA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17COOA2 5.145175e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17LUND4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17MUNG4 2.529822e-08        NA      NA     B  NA
asr.rr4 <- update(asr.rr4) #+40 updates to get % ch smaller than 1% was 

models.asr[["rr4ID"]] <- asr.rr4
models.asr[["rr4ID"]]$modelstats <- cbind.data.frame(niter = 35 , mniter = 20)

#-----------------
## %vaf
#-----------------
vv <- asr.rr4$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=4)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
gvar <- diag(Gmat) 
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/diag(Gmat)) # [1] 0.7697992    *** sum of the means
sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) # 0.7650199  *** mean of the sums

diff[[4]]<- c('FA4', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

hist(diag(lam%*%t(lam))/diag(Gmat))
models.asr[["rr4ID"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(diag(Gmat))
models.asr[["rr4ID"]]$"site%vaf" <- diag(lam%*%t(lam))/diag(Gmat)

# vv <- as.data.frame(vv)
# vv$Component <- rownames(vv)
# names(vv)[1] <- "scratch" 
# write.csv(vv[,c(2,1)], "asr.rr4-varcomp-scratch.csv", row.names = F)
# 
# scratchvaf <- models.asr[["rr4ID"]]$"site%vaf"
# scratchvaf <- as.data.frame(scratchvaf)
# scratchvaf$Experiment <- rownames(scratchvaf)
# names(scratchvaf)[1] <- "scratch"
# write.csv(scratchvaf[,c(2,1)], "asr.rr4-vaf-scratch.csv", row.names = F)


save(list = ls(), file = "1RData/RMETanalysis.RData")
####################################################################################################
## Model 6: rr5 + diag
#####################################################################################################
rr5.sv <- asreml(yield ~ Experiment + 
                   at(Experiment, mt$animaldmg):animaldmg + 
                   at(Experiment, mt$est):est + 
                   at(Experiment, mt$earlygs):earlygs + 
                   at(Experiment, mt$lcol):lin(Range) + 
                   at(Experiment, mt$lrow):lin(Row),
                 sparse =~ Experiment:VarietyDrop, 
                 random = ~ rr(Experiment, 5):VarietyKeep + diag(Experiment):VarietyKeep + 
                   at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                   at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                 residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                 na.action = na.method(y='include', x='include'), data=Rnvtdata, start.values = T)
#Get starting values set up
rr5.temp <- rr5.sv$vparameters.table

#Get values from previous model
rr4.gam <- matrix(summary(models.asr[['rr4ID']], vparameters=T)$vparameters[['rr(Experiment, 4):VarietyKeep']],ncol=5)
psi.gam <- matrix(summary(models.asr[['rr4ID']], vparameters=T)$vparameters[['Experiment:VarietyKeep']],ncol=1)
rr4.gam <- cbind(psi.gam,rr4.gam)
dimnames(rr4.gam) <- list(levels(Rnvtdata$Experiment),c('psi','rr.psi','lam1', 'lam2', 'lam3', 'lam4'))


#Put values from previous model in the starting values format.
# the structure is 1:ne = rr.var, (ne+1):2*ne = rr.fa2, (2*ne + 1):3*ne = rr.fa2, (3*ne+1):4*ne = psi's, and then other terms
#Put values from previous model in the starting values format.
rr5.temp$Value[grep('.*var', rr5.temp$Component)] <- 0 # set rr psi's to zero # 0.2*rr4.gam for fa1
rr5.temp$Value[grep('fa1', rr5.temp$Component)] <- rr4.gam[,"lam1"]
rr5.temp$Value[grep('fa2', rr5.temp$Component)] <- rr4.gam[,"lam2"]
rr5.temp$Value[grep('fa3', rr5.temp$Component)] <- rr4.gam[,"lam3"]
rr5.temp$Value[grep('fa4', rr5.temp$Component)] <- rr4.gam[,"lam4"]
rr5.temp$Value[grep('fa5', rr5.temp$Component)] <- c(0,0,0,0,rep(0.02,ne-4))
rr5.temp$Value[grep('Experiment:VarietyKeep',rr5.temp$Component)] <- 0.8*rr4.gam[,"psi"] #for rr + rr5 # doesnt exist for fa1

tmp <- cbind.data.frame(Component = names(models.asr[["rr4ID"]]$vparameters), Value = models.asr[["rr4ID"]]$vparameters)
tmp$Value[abs(tmp$Value)<0.00001] <- 1e-4
rr5.SV <- merge(rr5.temp, tmp, by = "Component", all=T)

rr5.SV$Value <- rr5.SV$Value.y
unique(is.na(rr5.SV$Value.y)); rr5.SV[is.na(rr5.SV$Value.y)==TRUE,]
rr5.SV$Value[is.na(rr5.SV$Value.y)] <- rr5.SV$Value.x[is.na(rr5.SV$Value.y)]
rr5.SV <- rr5.SV[, c("Component", "Value", "Constraint")]
rr5.SV <- rr5.SV[order(rr5.SV$Component),]
rr5.SV$Component <- as.character(rr5.SV$Component) 
rownames(rr5.SV) <- NULL
rownames(rr5.SV) <- rr5.SV$Component
unique(rr5.SV$Component %in% rr5.temp$Component)
rr5.SV <- rr5.SV[rr5.temp$Component,]
rownames(rr5.SV) <- NULL

#Set any small psi's, probably bound terms to something not so small
rr5.SV$Value[abs(rr5.SV$Value)<0.00001 & rr5.SV$Constraint != "F"] <- 1e-4
save.image()

#######################
## > asr.rr5$loglik
## [1] 3929.821 first iteration of rr5+diag was 12970.48
## so all OK
## 3.3 sec per iteration, 67iterations before converged
## this is fine
## 
save(list = ls(), file = "1RData/RMETanalysis.RData")
asr.rr5 <- asreml(yield ~ Experiment + 
                    at(Experiment, mt$animaldmg):animaldmg + 
                    at(Experiment, mt$est):est + 
                    at(Experiment, mt$earlygs):earlygs + 
                    at(Experiment, mt$lcol):lin(Range) + 
                    at(Experiment, mt$lrow):lin(Row),
                  sparse =~ Experiment:VarietyDrop, 
                  random = ~ rr(Experiment, 5):VarietyKeep + diag(Experiment):VarietyKeep + 
                    at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                    at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                  residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                  na.action = na.method(y='include', x='include'), data=Rnvtdata, G.param = rr5.SV, R.param = rr5.SV)
#Do this as a for loop so that you don't lose work if an update crashes mid-way. Always save the asr object before running a 2nd loop!
save.image()
for(i in 1:100){#x1
  print(i)
  #if(asr.rr5$converge) stop("converged - yaay")
  asr.rr5 <- update(asr.rr5)
}
#1 loop , finished on 6th of the 2nd loop (13 iterations per loop set), but then lots of updates to make %ch < 1
asr.rr5$loglik # 12970.48, 12997.11 (converged), 12997.56
asr.rr5$converge
subset(summary(asr.rr5)$varcomp, summary(asr.rr5)$varcomp$"%ch">5)
subset(summary(asr.rr5)$varcomp, bound == "B")
# Experiment:VarietyKeep!Experiment_WMaA13BILO4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13BROO4 2.866762e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA13COOL2 4.620360e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13DUAR4 1.087070e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA13MERR2 2.530310e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14EMER4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GILG2 1.408766e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14NYNG2 1.198422e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14ROMA4 3.544423e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14SURA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15DUAR4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15MUNG4 3.829349e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15NIND4 5.005403e-09        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15ROMA4 2.678163e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15YELA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17COOA2 4.128400e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA17LUND4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17MUNG4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17WALG2 6.236078e-08        NA      NA     B  NA
asr.rr5 <- update(asr.rr5) #+10 updates to get % ch smaller than 1% was 

models.asr[["rr5ID"]] <- asr.rr5
models.asr[["rr5ID"]]$modelstats <- cbind.data.frame(niter = 20, mniter = 33)

#--------------------
## %vaf
#--------------------
vv <- asr.rr5$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=5)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
gvar <- diag(Gmat) 
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/diag(Gmat)) # [1] 0.8167156   *** sum of the means
sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) # 0.8185117  *** mean of the sums

diff[[5]]<- c('FA5', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

hist(diag(lam%*%t(lam))/diag(Gmat))
models.asr[["rr5ID"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(diag(Gmat))
models.asr[["rr5ID"]]$"site%vaf" <- diag(lam%*%t(lam))/diag(Gmat)

sort(models.asr[["rr5ID"]]$"site%vaf")

save(list = ls(), file = "1RData/frostid.RData")
#load("nofrostid.RData")

# vv <- as.data.frame(vv)
# vv$Component <- rownames(vv)
# names(vv)[1] <- "scratch" 
# write.csv(vv[,c(2,1)], "asr.rr5-varcomp-scratch.csv", row.names = F)
# 
# scratchvaf <- models.asr[["rr5ID"]]$"site%vaf"
# scratchvaf <- as.data.frame(scratchvaf)
# scratchvaf$Experiment <- rownames(scratchvaf)
# names(scratchvaf)[1] <- "scratch"
# write.csv(scratchvaf[,c(2,1)], "asr.rr5-vaf-scratch.csv", row.names = F)

####################################################################################################
## Model 6: rr6 + diag
#####################################################################################################
rr6.sv <- asreml(yield ~ Experiment + 
                   at(Experiment, mt$animaldmg):animaldmg + 
                   at(Experiment, mt$est):est + 
                   at(Experiment, mt$earlygs):earlygs + 
                   at(Experiment, mt$lcol):lin(Range) + 
                   at(Experiment, mt$lrow):lin(Row),
                 sparse =~ Experiment:VarietyDrop, 
                 random = ~ rr(Experiment, 6):VarietyKeep + diag(Experiment):VarietyKeep + 
                   at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                   at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                 residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                 na.action = na.method(y='include', x='include'), data=Rnvtdata, start.values = T)
#Get starting values set up
rr6.temp <- rr6.sv$vparameters.table

#Get values from previous model
rr5.gam <- matrix(summary(models.asr[['rr5ID']], vparameters=T)$vparameters[['rr(Experiment, 5):VarietyKeep']],ncol=6)
psi.gam <- matrix(summary(models.asr[['rr5ID']], vparameters=T)$vparameters[['Experiment:VarietyKeep']],ncol=1)
rr5.gam <- cbind(psi.gam,rr5.gam)
dimnames(rr5.gam) <- list(levels(Rnvtdata$Experiment),c('psi','rr.psi','lam1', 'lam2', 'lam3', 'lam4','lam5'))


#Put values from previous model in the starting values format.
# the structure is 1:ne = rr.var, (ne+1):2*ne = rr.fa2, (2*ne + 1):3*ne = rr.fa2, (3*ne+1):4*ne = psi's, and then other terms
#Put values from previous model in the starting values format.
rr6.temp$Value[grep('.*var', rr6.temp$Component)] <- 0 # set rr psi's to zero # 0.2*rr5.gam for fa1
rr6.temp$Value[grep('fa1', rr6.temp$Component)] <- rr5.gam[,"lam1"]
rr6.temp$Value[grep('fa2', rr6.temp$Component)] <- rr5.gam[,"lam2"]
rr6.temp$Value[grep('fa3', rr6.temp$Component)] <- rr5.gam[,"lam3"]
rr6.temp$Value[grep('fa4', rr6.temp$Component)] <- rr5.gam[,"lam4"]
rr6.temp$Value[grep('fa5', rr6.temp$Component)] <- rr5.gam[,"lam5"]
rr6.temp$Value[grep('fa6', rr6.temp$Component)] <- c(0,0,0,0,0,rep(0.02,ne-5))
rr6.temp$Value[grep('Experiment:VarietyKeep',rr6.temp$Component)] <- 0.8*rr5.gam[,"psi"] #for rr + rr6 # doesnt exist for fa1

tmp <- cbind.data.frame(Component = names(models.asr[["rr5ID"]]$vparameters), Value = models.asr[["rr5ID"]]$vparameters)
tmp$Value[abs(tmp$Value)<0.00001] <- 1e-4
rr6.SV <- merge(rr6.temp, tmp, by = "Component", all=T)

rr6.SV$Value <- rr6.SV$Value.y
unique(is.na(rr6.SV$Value.y)); rr6.SV[is.na(rr6.SV$Value.y)==TRUE,]
rr6.SV$Value[is.na(rr6.SV$Value.y)] <- rr6.SV$Value.x[is.na(rr6.SV$Value.y)]
rr6.SV <- rr6.SV[, c("Component", "Value", "Constraint")]
rr6.SV <- rr6.SV[order(rr6.SV$Component),]
rr6.SV$Component <- as.character(rr6.SV$Component) 
rownames(rr6.SV) <- NULL
rownames(rr6.SV) <- rr6.SV$Component
unique(rr6.SV$Component %in% rr6.temp$Component)
rr6.SV <- rr6.SV[rr6.temp$Component,]
rownames(rr6.SV) <- NULL

#Set any small psi's, probably bound terms to something not so small
rr6.SV$Value[abs(rr6.SV$Value)<0.00001 & rr6.SV$Constraint != "F"] <- 1e-4
save.image()

#######################
## > asr.rr6$loglik
## [1] 3929.821 first iteration of rr6+diag was 12970.48
## so all OK
## 3.3 sec per iteration, 67iterations before converged
## this is fine
## 
save(list = ls(), file = "1RData/RMETanalysis.RData")
asr.rr6 <- asreml(yield ~ Experiment + 
                    at(Experiment, mt$animaldmg):animaldmg + 
                    at(Experiment, mt$est):est + 
                    at(Experiment, mt$earlygs):earlygs + 
                    at(Experiment, mt$lcol):lin(Range) + 
                    at(Experiment, mt$lrow):lin(Row),
                  sparse =~ Experiment:VarietyDrop, 
                  random = ~ rr(Experiment, 6):VarietyKeep + diag(Experiment):VarietyKeep + 
                    at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                    at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                  residual = ~ dsum(~id(Range):id(Row)| Experiment, levels = unlist(mt$resid)),
                  na.action = na.method(y='include', x='include'), data=Rnvtdata, G.param = rr6.SV, R.param = rr6.SV)
#Do this as a for loop so that you don't lose work if an update crashes mid-way. Always save the asr object before running a 2nd loop!
save.image()
for(i in 1:100){#x1
  print(i)
  if(asr.rr6$converge) stop("converged - yaay")
  asr.rr6 <- update(asr.rr6)
}
save.image()
for(i in 1:100){#x1
  print(i)
  #if(asr.rr6$converge) stop("converged - yaay")
  asr.rr6 <- update(asr.rr6)
}
save.image()
#1 loop , finished on 6th of the 2nd loop (13 iterations per loop set), but then lots of updates to make %ch < 1
asr.rr6$loglik # 13135.23, 13190.50 (converged), 13191.45
asr.rr6$converge
subset(summary(asr.rr6)$varcomp, summary(asr.rr6)$varcomp$"%ch">5)
subset(summary(asr.rr6)$varcomp, bound == "B")
# Experiment:VarietyKeep!Experiment_WMaA13BILO4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13COOL2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13DUAR4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13GOON2 7.933422e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13LUND4 5.191984e-09        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13MERR2 6.698633e-09        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13SPRS4 4.904997e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13WEST4 6.161902e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14CAPE4 8.371083e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14COOA2 1.221766e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14EMER4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GILG2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14NORT2 7.300333e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14NYNG2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14ROMA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14SURA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15DUAR4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15MUNG4 3.585868e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15NIND4 2.955653e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15ROMA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15TRAN2 1.759885e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15YELA4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA16WEST4 1.480175e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA17COOA2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17LUND4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17MUNG4 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17WALG2 2.529822e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17WEST4 3.527072e-07        NA      NA     B   0
asr.rr6 <- update(asr.rr6) #+10 updates to get % ch smaller than 1% was 

models.asr[["rr6ID"]] <- asr.rr6
models.asr[["rr6ID"]]$modelstats <- cbind.data.frame(niter = 481 , mniter = 40)

#--------------------
## %vaf
#--------------------
vv <- asr.rr6$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=6)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
gvar <- diag(Gmat) 
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/diag(Gmat)) # [1] 0.8544735   *** sum of the means
sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) # 0.844664  *** mean of the sums

diff[[6]]<- c('FA6', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

hist(diag(lam%*%t(lam))/diag(Gmat))
models.asr[["rr6ID"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(diag(Gmat))
models.asr[["rr6ID"]]$"site%vaf" <- diag(lam%*%t(lam))/diag(Gmat)

sort(models.asr[["rr6ID"]]$"site%vaf")

save(list = ls(), file = "1RData/frostid.RData")
#load("nofrostid.RData")

# vv <- as.data.frame(vv)
# vv$Component <- rownames(vv)
# names(vv)[1] <- "scratch" 
# write.csv(vv[,c(2,1)], "asr.rr6-varcomp-scratch.csv", row.names = F)
# 
# scratchvaf <- models.asr[["rr6ID"]]$"site%vaf"
# scratchvaf <- as.data.frame(scratchvaf)
# scratchvaf$Experiment <- rownames(scratchvaf)
# names(scratchvaf)[1] <- "scratch"
# write.csv(scratchvaf[,c(2,1)], "asr.rr6-vaf-scratch.csv", row.names = F)


#############################
##
## rr6 + diag with AR1 x AR1
##
############################### 
## 
rr6ar.sv <-asreml(yield ~ Experiment + 
                    at(Experiment, mt$animaldmg):animaldmg + 
                    at(Experiment, mt$est):est + 
                    at(Experiment, mt$earlygs):earlygs + 
                    at(Experiment, mt$lcol):lin(Range) + 
                    at(Experiment, mt$lrow):lin(Row),
                       sparse =~ Experiment:VarietyDrop, 
                       random = ~ rr(Experiment, 6):VarietyKeep + diag(Experiment):VarietyKeep + 
                    at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                    at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                  residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                         dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                         dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                         dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                       na.action = na.method(y='include', x='include'), data=Rnvtdata, start.values = T)

rr6.sv <- rr6ar.sv$vparameters.table
rownames(rr6.sv) <- rr6.sv$Component

# Get starting values from rr6ID model and correlations from diag ID
rr6.allID <- models.asr[['rr6ID']]$vparameters #rr6IDwas an ID model.
diag.cor <- models.asr[['diagAR1']]$vparameters #diag was an AR1 model
diag.cor <- diag.cor[grep("cor", names(diag.cor))]
rr6.all <- c(rr6.allID,diag.cor)
# first check terms match up
temp <- names(rr6.all)[!names(rr6.all) %in% rr6.sv$Component]; temp[grep("VarietyKeep", temp, invert=T)] # ok pretty good. Should be 0
temp <- rr6.sv$Component[!rr6.sv$Component %in% names(rr6.all)]; temp[grep("VarietyKeep", temp, invert=T)] # ok pretty good. SHould be 0
# then assign
rr6.all <- data.frame(rr6.all)
rr6.sv[rownames(rr6.all)[rownames(rr6.all)%in%rr6.sv$Component],"Value"] <- 
  rr6.all[rownames(rr6.all)%in%rr6.sv$Component,]

# lambdas/psis
rr6.gam <- rr6.all[grep('Experiment:VarietyKeep', rownames(rr6.all)),]
rr6.sv$Value[grep('fa1', rr6.sv$Component)] <- rr6.all[grep('fa1', rownames(rr6.all)),]
rr6.sv$Value[grep('fa2', rr6.sv$Component)] <- rr6.all[grep('fa2', rownames(rr6.all)),]
rr6.sv$Value[grep('fa3', rr6.sv$Component)] <- rr6.all[grep('fa3', rownames(rr6.all)),]
rr6.sv$Value[grep('fa4', rr6.sv$Component)] <- rr6.all[grep('fa4', rownames(rr6.all)),]
rr6.sv$Value[grep('fa5', rr6.sv$Component)] <- rr6.all[grep('fa5', rownames(rr6.all)),]
rr6.sv$Value[grep('Experiment:VarietyKeep',rr6.sv$Component)] <- 0.8*rr6.gam

#Set any small values to something not so small
rr6.sv$Value[abs(rr6.sv$Value) < 1e-5 & rr6.sv$Constraint != "F"] <- .001
rownames(rr6.sv) <- NULL
save.image()

#######################
## asr.rr6$loglik #3984.813
## rr6 + diag ar1xar1 run
##  8.8 sec/iternation
## starting logl 4087.532
#load("all.RData")
asreml.options(workspace = "2000mb", pworkspace="2000mb")
asr.rr6ar <- asreml(yield ~ Experiment + 
                      at(Experiment, mt$animaldmg):animaldmg + 
                      at(Experiment, mt$est):est + 
                      at(Experiment, mt$earlygs):earlygs + 
                      at(Experiment, mt$lcol):lin(Range) + 
                      at(Experiment, mt$lrow):lin(Row),
                    sparse =~ Experiment:VarietyDrop, 
                    random = ~ rr(Experiment, 6):VarietyKeep + diag(Experiment):VarietyKeep + 
                      at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                      at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                    residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                      dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                      dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                      dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                    na.action = na.method(y='include', x='include'), data=Rnvtdata, 
                    G.param = rr6.sv, R.param = rr6.sv)

save.image()
save(list = ls(), file = "1RData/frostid.RData")
for(i in 1:100){#x1
  print(i)
  if(asr.rr6ar$converge) stop ("Yaay - converged")
  asr.rr6ar <- update(asr.rr6ar)
}

save.image()
save(list = ls(), file = "1RData/frostid.RData")
for(i in 1:100){#x1
  print(i)
  # if(asr.rr6ar$converge) stop ("Yaay - converged")
  asr.rr6ar <- update(asr.rr6ar)
}

#4 updates converging on iteration 11 of update 4
asr.rr6ar$loglik # 14982.96, 15012.46 (converged), 15013.22
subset(summary(asr.rr6ar)$varcomp, summary(asr.rr6ar)$varcomp$"%ch">5)
subset(summary(asr.rr6ar)$varcomp, bound == "B")
# Experiment:VarietyKeep!Experiment_WMaA13BILO4 2.529822e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA13COOL2 5.918973e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13DUAR4 4.647911e-09        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA13GOON2 2.779236e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14CAPE4 1.339717e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14GOON2 4.386305e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14KING4 8.063421e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA14ROMA4 2.529822e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA14YELA4 9.797307e-09        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15NIND4 1.411419e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15ROMA4 2.529822e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15TRAN2 8.361005e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA15TULO2 3.096920e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA16WEST4 3.021283e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17LUND4 2.384787e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA17MEAN4 6.982427e-08        NA      NA     B  NA
# Experiment:VarietyKeep!Experiment_WMaA17NIND4 1.039051e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA17NYNG2 5.250170e-08        NA      NA     B  NA
asr.rr6ar <- update(asr.rr6ar)
subset(summary(asr.rr6ar)$varcomp, summary(asr.rr6ar)$varcomp$"%ch">5)


models.asr[["rr6AR1"]] <- asr.rr6ar
models.asr[["rr6AR1"]]$modelstats <- cbind.data.frame(niter = 167, mniter = 113) # 12*13 +11
asr.rr6ar$loglik #15013.22

save(list = ls(), file = "1RData/frostid.RData")

#--------------------
## %vaf
#--------------------
vv <- asr.rr6ar$vparameters
avar <- vv[grep('VarietyKeep',names(vv))]
lam <- matrix(avar[grep('fa',names(avar))],ncol=6)
psi <- avar[1:ne]
Gmat <- lam%*%t(lam) + diag(psi)
gvar <- diag(Gmat) 
dimnames(Gmat) <- list(enam,enam)
mean(diag(lam%*%t(lam))/diag(Gmat)) # [1] 0.8489626   *** sum of the means
sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) # 0.8426813  *** mean of the sums

diff[[7]]<- c('FA6ar', mean(diag(lam%*%t(lam))/diag(Gmat))-sum(diag(lam %*% t(lam)))/sum(diag(Gmat)))

hist(diag(lam%*%t(lam))/diag(Gmat))
models.asr[["rr6AR1"]]$"%vaf" <- sum(diag(lam %*% t(lam)))/sum(diag(Gmat))
models.asr[["rr6AR1"]]$"site%vaf" <- diag(lam%*%t(lam))/diag(Gmat)
sort(models.asr[["rr6AR1"]]$"site%vaf")

save(list = ls(), file = "1RData/frostid.RData")


###############################-----------------------------------------------------------------
#Output LOGL and %varexpl
save(list = c("models.asr"), file = paste0("1RData/", met.name, "-frostid-allmodels.RData"))

save(list = c("Rnvtdata", "Rnvtmodels", "asr.rr6ar", "met.name"), file = paste0("1RData/", met.name, "-frostid",  "-RMETanalysisOUT.RData"))

save(list = ls(), file = "frostid.RData")


########################

## end of script  #NOW GO TO 04-METout_v2.R

##########################

save.image()
