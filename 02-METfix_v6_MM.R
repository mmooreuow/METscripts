##### Script Information ---------------------------------------------------------------------------------------------------------------------------
#
# Title: 02-METfix_v4.R
# Author: Mandy Moore
# Dataset: Wheat North MET including frosted trials
# Date: September 2018
# Description: 
# - Prepares the MET dataset for analysis
# - Provides some summaries to look at but not for the report.
# - The tables in the report will be based on the analysis dataset.
#
# Modifications:
# 1) Fix of METInclude and METinclude and SAGIInclude and SAGIinclude and SAGIcomments.
# 2) Streamlining covariate tests
#
# Input: Data from 01-DataCheck.R
# Children: 03-METAnalysis-id.R

###### Data importation and package initialisation ----------------------------------------------------------------------------------------------------
rm(list =ls())

# Intialise packages
require(plyr)
require(dplyr)
require(reshape2)
require(asreml)
require(ASExtras4)
require(ggplot2)
require(XLConnect)

# Required for mac:
#dyn.load("/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/lib/server/libjvm.dylib")

# Required for model.fit:
source("1Functions/NVTfns5.6.R")  
setup.fn(computer='linux')

# Set-up MET details
met.name <- "WheatMainNorth"
date <- "2018-09-01"
date.out <- "2018-09-01"

# Import data from 01-DataCheckR
load(paste0("1RData/", met.name, "-", date, "-data4MET.RData"))

# Import trial summary data from csv
trial.info <- read.csv("1RawData/Export_Wheat_Summary_2018-08-09.csv")

###### Data Checks and Set-up -------------------------------------------------------------------------------------------------------------------------

# Make sure the models and data file have the same number of experiments
n_distinct(Rnvtdata$Experiment) # 137
n_distinct(Rnvtmodels$Experiment) # 137

# Generate the factors
str(Rnvtdata)

factor.names <- names(Rnvtdata)[1:(grep("Harvest.Length", names(Rnvtdata))-1)]
for(i in 1:length(factor.names)){
  Rnvtdata[,factor.names[i]] <- factor(Rnvtdata[,factor.names[i]])}

str(Rnvtdata)

(ne.init <- nlevels(Rnvtdata$Experiment)) #137
(nv.init <- nlevels(Rnvtdata$Variety)) #132

# Check if only MET include sites
met.incl <- unique(Rnvtdata[, c("Experiment", "METInclude")])
with(met.incl, table(METInclude)) 
#   METInclude
# FALSE  TRUE 
#     8   129 

# Set up trial.info (Data obtained from trial summary csv)
str(trial.info)
names(trial.info)[1:2] <- c("Experiment", "Year")
trial.info <- droplevels(subset(trial.info, Experiment %in% unique(Rnvtdata$Experiment)))
table(trial.info$SAGIInclude) #109
trial.info$SAGIInclude[is.na(trial.info$SAGIInclude)==TRUE] <- TRUE
trial.info$SAGIcomment[is.na(trial.info$SAGIcomment)==TRUE] <- "NONE"
unique(trial.info$SAGIcomment)

# Check which trials have been excluded by TSP's
with(trial.info, table(METInclude))
# METInclude
# FALSE  TRUE 
#    8   129
trial.info[!trial.info$METInclude, c('Experiment', 'AOVComment')]
# Experiment                                          AOVComment
# WMaA13BUNG4 Severe frost damage - Flagged invalid on 23/12/2016
# WMaA13COOA2                                        Site frosted
# WMaA14BELL2                   Severe frost damage at this site.
# WMaA17BELL2                            Trial affected by frost.
# WMaA17COOL2                            Trial affected by frost.
# WMaA17GOON2                            Trial affected by frost.
# WMaA17WALG2                            Trial affected by frost.
# WMaA17WONG2                            Trial affected by frost.

###### Variety Checks ---------------------------------------------------------------------------------------------------------------------------------

# Check that Fillers are NA's
unique(subset(Rnvtdata, substring(tolower(Variety), 1, 4)=="fill")$yield) 
subset(Rnvtdata, substring(tolower(Variety), 1, 4)=="fill" & is.na(yield)==FALSE)$Experiment
 # If there are some Filler values that are not NA run the following code:
 Rnvtdata$yield[substring(tolower(Rnvtdata$Variety), 1, 4)=="fill"] <- NA
 # Check again
 unique(subset(Rnvtdata, substring(tolower(Variety), 1, 4)=="fill")$yield) 
 
# Variety by Year contingency table
vdata <- unique(Rnvtdata[, c("Variety", "Year")])
vfreq <- with(vdata, table(Variety, Year))

 # Check
 rowSums(vfreq)[rowSums(vfreq) ==0]
 # Format vfreq. This will be part of the output given to the breeder in the prelim data check file!
 vfreq <- as.data.frame(vfreq)
 vfreq <- subset(vfreq, substring(tolower(Variety), 1, 4)!="fill")
 (vfreq <- dcast(vfreq, Variety ~ Year, value.var = "Freq"))

 ###### Check Reps for discrepancies ------------------------------------------------------------------------------------------------------------------
 
# When designs are generated ColRep is equal to Rep. We can check changes in variety allocation
# to trials by checking if Rep==ColRep and flag any trials where this does not.
# Invesitage and consult TSPs/breeders (if required) if Reps==ColReps. 
# Rep IS NOT part of the plot structure factors!
 
table(as.character(Rnvtdata$Rep)==as.character(Rnvtdata$ColRep))
# FALSE  TRUE 
#    14 16378
Rnvtdata[!as.character(Rnvtdata$Rep)==as.character(Rnvtdata$ColRep),c('Experiment','Rep','ColRep','Variety')]
#  Experiment Rep ColRep      Variety
# WMaA17MERR2   2      1       Filler
# WMaA17MERR2   3      1       Filler
# WMaA17MERR2   4      2       Filler
# WMaA17MERR2   5      2       Filler
# WMaA17MERR2   6      2       Filler
# WMaA17MERR2   7      3       Filler
# WMaA17MERR2   8      3       Filler
# WMaA17MERR2   9      3       Filler
# WMaA17SPRS4   2      1   DS Faraday
# WMaA17SPRS4   3      2   DS Faraday
# WMaA17SPRS4   3      2 LRPB Reliant
# WMaA17SPRS4   4      3   DS Faraday
# WMaA17SPRS4   4      3 LRPB Reliant
# WMaA17SPRS4   4      3       Coolah

rep.temp <- reps.check(data=Rnvtdata, BlockID="Rep", outfile1=paste0("1Graphics/", met.name, "-", date,"-layout-RepFixed.pdf"),
                       outfile2=paste0("1Graphics/", met.name, "-", date,"-VarietyFreqperRepFixed.pdf"))
(rep.na <- rep.temp$rep.na)
(rep.bad <- rep.temp$rep.bad) 

# Check plots in output file.
# Comment: Some of the descrepancy is due to fillers being used in the experiment
# and reps of varieties have been swapped in rows.

  # Specific Experiment Investigations
  # Experiments: WMaA17MERR2, WMaA17SPRS4
  # Checks:
  tmp <- droplevels(subset(Rnvtdata, Rnvtdata$Experiment=='WMaA17MERR2'))
  table(tmp$Variety, tmp$Rep) # All good, just filler
  tmp <- droplevels(subset(Rnvtdata, Rnvtdata$Experiment=='WMaA17SPRS4'))
  table(tmp$Variety, tmp$Rep)
  apply(table(tmp$Variety, tmp$Rep), 1, sum) # Not Filler, Coolah 4, DS Faraday 4, LRPB Reliant 4
  tmp[tmp$Variety=='Coolah',c('Experiment','Range', 'Row', 'Rep','ColRep', 'RowRep', 'Variety')]
  tmp[tmp$Variety=='DS Faraday',c('Experiment','Range', 'Row', 'Rep','ColRep', 'RowRep', 'Variety')]
  tmp[tmp$Variety=='LRPB Reliant',c('Experiment','Range', 'Row', 'Rep','ColRep', 'RowRep', 'Variety')]

levels(Rnvtdata$Rep)

###### Check Colreps are aligned with logical blocks --------------------------------------------------------------------------------------------------
# check that ColReps are in same directions for all experiments within co-located trials.

crep.temp <- reps.check(data=Rnvtdata, BlockID="ColRep", outfile1=paste0("1Graphics/", met.name, "-", date,"-layout-ColRepFixed.pdf"),
                       outfile2=paste0("1Graphics/", met.name, "-", date,"-VarietyFreqperColRepFixed.pdf"))
(rep.na <- crep.temp$rep.na)
(rep.bad <- crep.temp$rep.bad) 

# Check plots in output file.
# Comment: NA

  # Specific Experiment Investigations
  # Experiments: NA
  # Checks:
  tmp <- droplevels(subset(Rnvtdata, Rnvtdata$Experiment==''))
  table(tmp$Variety, tmp$ColRep)
  apply(table(tmp$Variety, tmp$ColRep), 1, sum) 

Rnvtdata$ColRep <- factor(Rnvtdata$ColRep)
nlevels(Rnvtdata$ColRep)

###### Check RowReps are aligned with logical blocks --------------------------------------------------------------------------------------------------
# check that RowReps are in same directions for all experiments within for co-located trials.

rrep.temp <- reps.check(data=Rnvtdata, BlockID="RowRep", outfile1=paste0("1Graphics/", met.name, "-", date,"-layout-RowRepOrig.pdf"),
                        outfile2=paste0("1Graphics/", met.name, "-", date,"-VarietyFreqperRowRepOrig.pdf"))
(rrep.na <- rrep.temp$rep.na)
(rrep.bad <- rrep.temp$rep.bad)

# Check plots in output file.
# Comment:Tetris pattern. No further investigation required.

  # Specific Experiment Investigations
  # Experiments: NA
  # Checks:
  tmp <- droplevels(subset(Rnvtdata, Rnvtdata$Experiment==''))
  table(tmp$Variety, tmp$RowRep)
  apply(table(tmp$Variety, tmp$RowRep), 1, sum) 

Rnvtdata$RowRep <- factor(Rnvtdata$RowRep)
nlevels(Rnvtdata$RowRep)

##### Model Set-up (DIAG model to test covariates and spatial terms) ----------------------------------------------------------------------------------
Rnvtdata <- unique(Rnvtdata[order(Rnvtdata$Experiment, Rnvtdata$Range, Rnvtdata$Row),])
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
asreml.options(workspace = "1000mb", pworkspace="1000mb")

###### Test covariates---------------------------------------------------------------------------------------------------------------------------------

# Previous versions of this script have this as a manual step
# It involved using the mt output to obtain the names of the experiments for which there are covariates.  
# This step now sorts itself out. 
not.cov <- c("lrow", "lcol", "crep", "rrep", "rrow", "rcol", "resid")
(mt.cov <- paste0('mt$',names(mt)[!names(mt)%in%not.cov]))
(mt.cov <- mt.cov[-1]) # removes mt$YrCon (leave in for Co-located trials)
(cov2test <- sort(unique(unlist(lapply(mt.cov, function(x) eval(parse(text = x)))))))

temp1 <- cov2test
tmpmodel <- droplevels(Rnvtmodels[Rnvtmodels$Experiment %in% temp1,])
tmpdata <- droplevels(Rnvtdata[Rnvtdata$Experiment %in% temp1,])

mt.tmp <- Colmodel.fit(models = tmpmodel, data = tmpdata)

paste0("at(Experiment, ", mt.cov, "):", unlist(lapply(strsplit(mt.cov, "$", fixed = TRUE), function(x) x[2])), collapse = " + ")
# Copy paste the output from the above command into the second line in the asreml call directly below
# Don't fit any spatial terms in the model for testing covariates
require(asreml)
tmp.diag <- asreml(yield ~ Experiment + 
                      at(Experiment, mt$animaldmg):animaldmg + at(Experiment, mt$est):est + at(Experiment, mt$earlygs):earlygs,
                   sparse=~ at(Experiment):Variety,
                   random = ~ at(Experiment, mt.tmp$crep):ColRep + at(Experiment, mt.tmp$rrep):RowRep, 
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt.tmp$resid$aa),
                   na.action = na.method(y='include', x='include'), data=tmpdata)
tmp.diag <- update(tmp.diag)


# OR simply run this command
# require(asreml)
tmp.diag <- eval(parse(text = paste0("asreml(yield ~ Experiment + ",
                   paste0("at(Experiment, ", mt.cov, "):", unlist(lapply(strsplit(mt.cov, "$", fixed = TRUE), function(x) x[2])), collapse = " + "),",",
                   "sparse=~ at(Experiment):Variety,
                   random = ~ at(Experiment, mt.tmp$crep):ColRep + at(Experiment, mt.tmp$rrep):RowRep,
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt.tmp$resid$aa),
                   na.action = na.method(y='include', x='include'), data=tmpdata)")))
tmp.diag <- update(tmp.diag)

# Note: We don't fit any spatial terms in the model for testing covariates!

# Check wald output for covariates. 
wald(tmp.diag, den.DF="default") 
# Comment: 
# Covariates included:
# at(Experiment, WMaA15DUAR4):animaldmg
# at(Experiment, WMaB17DUAR4):animaldmg
# at(Experiment, WMaA17WEST4):est     
# at(Experiment, WMaA17WONG2):est     
# at(Experiment, WMaA17COOL2):earlygs 
# Covariates removed: NA

# Remove from the Rnvtmodels dataframe if they are not significant.
# e.g est is not significant for experiment WMaA14QUAN2 , remove it.
# Rnvtmodels$est[Rnvtmodels$Experiment=="WMaA14QUAN2"] <- 0

##### Test lin spatial terms --------------------------------------------------------------------------------------------------------------------------

# Identify which trials have < 4 rows or columns
env.sum <- Rnvtdata %>% group_by(Experiment) %>% summarise(ncol = n_distinct(Range),
                                                      nrow = n_distinct(Row))
(exp.col4 <- droplevels(env.sum$Experiment[env.sum$ncol < 4])) 
# 9 exp: WMaA13CAPE4 WMaA13NYNG2 WMaA14CAPE4 WMaA14JAMB4 WMaA14NYNG2 WMaA15THEO4 WMaA17CAPE4 WMaA17JAMB4 WMaA17KILC4
(exp.row4 <- droplevels(env.sum$Experiment[env.sum$nrow < 4])) 
# 0 exp

# Check that if a range has max 3 then set to id for residual and make sure no lincol or rcol terms are fitted
Rnvtmodels$resid <- as.character(Rnvtmodels$resid)

Rnvtmodels$resid[Rnvtmodels$Experiment %in% exp.col4] #all 'ia'
Rnvtmodels$lcol[Rnvtmodels$Experiment %in% exp.col4] #all 0
Rnvtmodels$rcol[Rnvtmodels$Experiment %in% exp.col4] #all 0

# Check that if a range/row has max 3 then set to id for residual and make sure no lincol or rcol terms are fitted
Rnvtmodels$resid[Rnvtmodels$Experiment %in% exp.row4] # NA
Rnvtmodels$lcol[Rnvtmodels$Experiment %in% exp.row4] # NA
Rnvtmodels$rcol[Rnvtmodels$Experiment %in% exp.row4] # NA

# If there is a linrow or lincol, include row (rrow) and range (rcol) in the model 
Rnvtmodels[Rnvtmodels$lrow==1,] 
Rnvtmodels[Rnvtmodels$lcol==1,]
Rnvtmodels$rrow[Rnvtmodels$lrow==1] <- 1
Rnvtmodels$rcol[Rnvtmodels$lcol==1] <- 1

# Manual fix - fix any that need special attention e.g.
#Rnvtmodels$resid[Rnvtmodels$Experiment=="WMaA13WARR5"] <- "ia"

# Subset the data in order to check lincol and linrow
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
lin2test <- sort(unique(c(mt$lcol, mt$lrow)))
lcol <- lrow <- NULL

#### CYCLE BEGIN --------------------------------------------------------------------------------------------------------------------------------------
source("1Functions/MMfns1.R")
# Decide on the number of covariates to test at a time (15 seems to be handled okay - may depend on dataset)
# Needs to be small enough for wald(asr, denDF="default") to work
n.cov <- 15
# Number of cycles
ncycles <- ceiling(length(lin2test)/n.cov)

lcol.Exp <- lrow.Exp <- NULL
for(i in 1:ncycles){
  a <- 1+(i-1)*n.cov
  b <- i*n.cov
  temp1 <- lin2test[a:b] # Will change from 1:15, to 16:30
  tmpmodel <- droplevels(Rnvtmodels[Rnvtmodels$Experiment %in% temp1,])
  tmpdata <- droplevels(Rnvtdata[Rnvtdata$Experiment %in% temp1,])
  
  mt.tmp <- Colmodel.fit(models = tmpmodel, data = tmpdata)
  c <- lintest.call(mt.tmp)
  tmp.diag <- eval(parse(text = c))
  tmp.diag <- update(tmp.diag)
  
  (tmp.wald <- wald(tmp.diag, denDF = "default")$Wald)
  (sig.lin <- row.names(tmp.wald)[tmp.wald$Pr<0.05])
  tmp.lcol <- gsub("):lin(Range)", "", gsub("at(Experiment, ", "", sig.lin[grep("Range", sig.lin)], fixed = TRUE), fixed = TRUE)
  tmp.lrow <- gsub("):lin(Row)", "", gsub("at(Experiment, ", "", sig.lin[grep("Row", sig.lin)], fixed = TRUE), fixed = TRUE)
  
  lcol.Exp <- unique(c(lcol.Exp, unique(c(lcol, tmp.lcol))))
  lrow.Exp <- unique(c(lrow.Exp, unique(c(lrow, tmp.lrow))))
}
lcol.Exp # 2
lrow.Exp # 7

#### CYCLE END ----------------------------------------------------------------------------------------------------------------------------------------

#Correct the models file...
#First the random terms, because they depend on the linear term's
#This allows the identifed random terms to be retained when they weren't associated with lin terms.
rcol.Exp <- Rnvtmodels$Experiment[Rnvtmodels$lcol==Rnvtmodels$rcol]
rrow.Exp <- Rnvtmodels$Experiment[Rnvtmodels$lrow==Rnvtmodels$rrow]

Rnvtmodels$rcol[Rnvtmodels$Experiment %in% rcol.Exp] <- 0
Rnvtmodels$rcol[Rnvtmodels$Experiment %in% lcol.Exp] <- 1 #This is deliberate because want to set rcol to 1 where there is a lcol term.
Rnvtmodels$rrow[Rnvtmodels$Experiment %in% rrow.Exp] <- 0
Rnvtmodels$rrow[Rnvtmodels$Experiment %in% lrow.Exp] <- 1 #This is deliberate because want to set rrow to 1 where there is a lrow term.

#Now the linear terms
Rnvtmodels$lcol <- 0
Rnvtmodels$lrow <- 0
Rnvtmodels$lcol[Rnvtmodels$Experiment %in% lcol.Exp] <- 1
Rnvtmodels$lrow[Rnvtmodels$Experiment %in% lrow.Exp] <- 1

# Once the lin row/random row has been corrected re-run the CYCLE BEGIN-CYCLE END section to check 
# that everything still runs

# Can use this block of code to manually check that the linear terms have a 
# corresponding random term associated with them.
temp1 <- lin2test
tmpmodel <- droplevels(Rnvtmodels[Rnvtmodels$Experiment %in% temp1,])
tmpdata <- droplevels(Rnvtdata[Rnvtdata$Experiment %in% temp1,])
mt.tmp <- Colmodel.fit(models = tmpmodel, data = tmpdata)

sum(Rnvtmodels$lcol); length(lcol) #2
sum(Rnvtmodels$lrow); length(lrow) #7
sum(Rnvtmodels$rcol) #19
sum(Rnvtmodels$rrow) #32

#### Check the residual spatial terms -----------------------------------------------------------------------------------------------------------------
# Summary:
# Any autocorrelation term that has a -ve componenet value then a corresponding random term
# is added to the model, after running the model again any autocorrelation terms 
# still -ve get set to id instead of ar1.

# A spatial model is determined by analysing the global and local spatial variation
# for each site following Gilmour et al. (1997). For co-located trials this involves
# fitting terms at either the trial or the environment level depending on the term as
# described in the vignette. For trials where the number of rows or columns is less
# than or equal to 3 the plots in these directions are assumed to be independent (i.e.
# no autoregressive process (ar1 is not fitted because the number of dimensions is too
# small). In addition, when the autocorrelations are estimated to be less than zero
# then an equivalent random term is also fitted (e.g. if ar1(Row) is negative then Row
# is fitted as a peripheral random effect.). If the autocorrelation estimate is still zero
# after including this random peripheral term then associated residual term is set to be 
# independent, as negative autocorrelations suggests some type of competition between plots 
# which we do not currently model in NVT.


Rnvtmodels0 <- Rnvtmodels #Just to keep what we did above somewhere secure

##### ROUND 1:----
# Initial model to check if autocorrelations are -ve
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
tmp.diag <- asreml(yield ~ Experiment + 
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs +
                     at(Experiment, mt$lcol):lin(Range) + 
                     at(Experiment, mt$lrow):lin(Row),
                   sparse=~ at(Experiment):Variety,
                   random = ~ at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   #ROUND1 & 2
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                                dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia),
                   #ROUND3
                   # residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                   #              dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                   #              dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                   #              dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)
tmp.diag <- update(tmp.diag) 
# ROUND 1:13+13+6 iterations take 0.9,0.8,0.7 sec with 6 (2 iterations each) updates to remove 1% change
# loglik = 9076.00

#### Fix all that have -ve autocorrelation ------------------------------------------------------------------------------------------------------------
# Fix them by fitting the associated random effect and checking if still -ve.
# In ROUND 2 if still -ve then set to id.

# Use ROUND 1 results: assocaite -ve autocorrelations with random row/col as appropriate 
# - look at variograms.
pdf("1Graphics/DIAG-variograms1.pdf")
met.asreml(tmp.diag)
dev.off()

tt <- summary(tmp.diag, vparameters = TRUE)$varcomp
tt[tt$component < 0,]

range2fix <- tt[grep("!Range!cor", row.names(tt)),]
range2fix <- range2fix[range2fix$component <0,]
range2fix <- gsub("!Range!cor", "", gsub("Experiment_", "", rownames(range2fix), fixed = TRUE), fixed = TRUE)
Rnvtmodels$resid <- as.character(Rnvtmodels$resid)
Rnvtmodels$rcol[Rnvtmodels$Experiment %in% range2fix] <- 1 
sum(Rnvtmodels$rcol[Rnvtmodels$Experiment %in% range2fix]) # round 1: 40; round 2: 38

row2fix <- tt[grep("!Row!cor", row.names(tt)),]
row2fix <- row2fix[row2fix$component <0,]
row2fix <- gsub("!Row!cor", "", gsub("Experiment_", "", rownames(row2fix), fixed = TRUE), fixed = TRUE)
Rnvtmodels$rrow[Rnvtmodels$Experiment %in% row2fix] <- 1
sum(Rnvtmodels$rrow[Rnvtmodels$Experiment %in% row2fix])  # round 1: 7; round 2: 8; round 4: 1.

##### ROUND 2:----
# Re-run the model to check if -ve correlations are still -ve
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
tmp.diag <- asreml(yield ~ Experiment + 
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs +
                     at(Experiment, mt$lcol):lin(Range) + 
                     at(Experiment, mt$lrow):lin(Row),
                   sparse=~ at(Experiment):Variety,
                   random = ~ at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   # ROUND 1 & 2
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                     dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia),
                   #ROUND3
                   # residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                   #              dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                   #              dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                   #              dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)
tmp.diag <- update(tmp.diag) 
# ROUND 1:13+13+6 iterations take 0.9,0.8,0.7 sec with 6 (2 iterations each) updates to remove 1% change
# ROUND 2: 13+13+6 iterations take 1.0,0.9,0.8 sec with 6 updates to remove 1% change
# loglik = 9076.00, 9119.021

#### Drop duplicate terms -----------------------------------------------------------------------------------------------------------------------------
# Ok, for some trials Column=ColRep and Row=RowRep and so fitting random term in direction of negative correlation already there
# so simply drop it down to id in that direction straight away

ss <- cbind(with(droplevels(subset(Rnvtdata, Experiment%in%range2fix)), tapply(Range, Experiment, function(x) length(unique(x)))),
             with(droplevels(subset(Rnvtdata, Experiment%in%range2fix)), tapply(ColRep, Experiment, function(x) length(unique(x)))))

crep.rm <- dimnames(ss)[[1]][ss[,1]==ss[,2]] #none

ss <- cbind(with(droplevels(subset(Rnvtdata, Experiment%in%row2fix)), tapply(Row, Experiment, function(x) length(unique(x)))),
             with(droplevels(subset(Rnvtdata, Experiment%in%row2fix)), tapply(RowRep, Experiment, function(x) length(unique(x)))))

rrep.rm <- dimnames(ss)[[1]][ss[,1]==ss[,2]] # none

Rnvtmodels$rcol[Rnvtmodels$Experiment%in%crep.rm] <- 0
Rnvtmodels$resid[Rnvtmodels$Experiment%in%crep.rm] <- "ia"

#### Fix all that still have -ve autocorrelation ------------------------------------------------------------------------------------------------------
# Fix them by setting to id.

# Use ROUND 2 results: setting the autocorrelations that are -ve to id.
tt <- summary(tmp.diag, vparameters = TRUE)$varcomp
cc <- tt[tt$component < 0,]

range2fix <- cc[grep("!Range!cor", row.names(cc)),]
range2fix <- range2fix[range2fix$component <0,]
range2fix <- gsub("!Range!cor", "", gsub("Experiment_", "", rownames(range2fix), fixed = TRUE), fixed = TRUE) # none
 Rnvtmodels$resid <- as.character(Rnvtmodels$resid)
 Rnvtmodels$resid[Rnvtmodels$Experiment %in% range2fix]
 Rnvtmodels$resid[Rnvtmodels$Experiment %in% range2fix & Rnvtmodels$resid=="aa"] <- "ia"
 Rnvtmodels$resid[Rnvtmodels$Experiment %in% range2fix & Rnvtmodels$resid=="ai"] <- "ii"

row2fix <- cc[grep("!Row!cor", row.names(cc)),]
row2fix <- row2fix[row2fix$component <0,]
row2fix <- gsub("!Row!cor", "", gsub("Experiment_", "", rownames(row2fix), fixed = TRUE), fixed = TRUE) # BMaA16WALG2
Rnvtmodels$resid[Rnvtmodels$Experiment %in% row2fix]
Rnvtmodels$resid[Rnvtmodels$Experiment %in% row2fix & Rnvtmodels$resid=="aa"] <- "ai"
Rnvtmodels$resid[Rnvtmodels$Experiment %in% row2fix & Rnvtmodels$resid=="ia"] <- "ii"

##### ROUND 3:----
# Re-run the model to check for bound row and range terms
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
tmp.diag <- asreml(yield ~ Experiment + 
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs +
                     at(Experiment, mt$lcol):lin(Range) + 
                     at(Experiment, mt$lrow):lin(Row),
                   sparse=~ at(Experiment):Variety,
                   random = ~ at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   # ROUND 1 & 2
                   # residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                   #   dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia),
                   #ROUND3
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                                dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                                dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                                dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)
tmp.diag <- update(tmp.diag) 
# Round 1:13+13+6 iterations take 0.9,0.8,0.7 sec with 6 (2 iterations each) updates to remove 1% change
# Round 2: 13+13+6 iterations take 1.0,0.9,0.8 sec with 6 updates to remove 1% change
# Round 3: 13+13+2 iterations take 0.9,0.7,0.7 sec with no updates to remove 1% change
#loglik = 9076.00, 9119.021, 9093.570

#### Check for bound Row and Range terms --------------------------------------------------------------------------------------------------------------
# remove bound random terms unless assocated with lin terms

pdf("1Graphics/DIAG-variograms3.pdf")
met.asreml(tmp.diag)
dev.off()

tt <- summary(tmp.diag, vparameters = TRUE)$varcomp
tt.range <- tt[grep(":Range", rownames(tt)),]
tt.rangeb <- gsub("):Range", "", gsub("at(Experiment, ", "", row.names(tt.range)[tt.range$bound=="B"], fixed = TRUE), fixed = TRUE) # none
rcol20 <- tt.rangeb[!tt.rangeb %in% Rnvtmodels$Experiment[Rnvtmodels$lcol==1]] 
Rnvtmodels$rcol[Rnvtmodels$Experiment %in% rcol20] <- 0

tt.row <- tt[grep(":Row", rownames(tt)),]
tt.rowb <- gsub("):Row", "", gsub("at(Experiment, ", "", row.names(tt.row)[tt.row$bound=="B"], fixed = TRUE), fixed = TRUE)
tt.rowb <- tt.rowb[grep("Rep", tt.rowb, invert = TRUE)] # Just experiments with Row terms not RowRep
rrow20 <- tt.rowb[!tt.rowb %in% Rnvtmodels$Experiment[Rnvtmodels$lrow==1]] # "BMaA16WALG2" "BMaA17GOON2"
Rnvtmodels$rrow[Rnvtmodels$Experiment %in% rrow20] <- 0

#### ROUND 4:----
# Re-run the model without the bound terms as a check
mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
tmp.diag <- asreml(yield ~ Experiment + 
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs +
                     at(Experiment, mt$lcol):lin(Range) + 
                     at(Experiment, mt$lrow):lin(Row),
                   sparse=~ at(Experiment):Variety,
                   random = ~ at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   # ROUND 1 & 2
                   # residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                   #   dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia),
                   # ROUND3
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                     dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                     dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                     dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii),
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)
tmp.diag <- update(tmp.diag) 
# Round 1:13+13+6 iterations take 0.9,0.8,0.7 sec with 6 (2 iterations each) updates to remove 1% change
# Round 2: 13+13+6 iterations take 1.0,0.9,0.8 sec with 6 updates to remove 1% change
# Round 3: 13+13+2 iterations take 0.9,0.7,0.7 sec with no extra updates to remove 1% change
# Round 4: 13+13+2 iterations take 0.8,0.6,0.6 sec with no extra updates to remove 1% change
#loglik = 9076.00, 9119.021, 9093.570, 9093.570 (expect last two to be the same as only removing bound terms (which are already set to zero))

#### Final check for -ve components -------------------------------------------------------------------------------------------------------------------
# Use ROUND 4 results:
tt <- summary(tmp.diag, vparameters = TRUE)$varcomp
tt[tt$component < 0,]
#                                 component std.error    z.ratio bound %ch
# Experiment_WMaA13BILO4!Row!cor -0.070947 0.1710074 -0.4148769     U 0.7

# Investigate -ve component
tt[grep("WMaA13BILO4", rownames(tt)),]
lapply(mt, function(x) grep("WMaA13BILO4",x))
lapply(mt$resid, function(x) grep("WMaA13BILO4",x))

# Add rrow term in model for "WMaA13BILO4"
Rnvtmodels$rrow[Rnvtmodels$Experiment%in%c("WMaA13BILO4")] <- 1

# Run ROUND 4 again:
tt <- summary(tmp.diag, vparameters = TRUE)$varcomp
tt[tt$component < 0,]
#                                 component std.error   z.ratio bound %ch
# Experiment_WMaA13BILO4!Row!cor -0.2917796 0.1518297 -1.921756     U   0
# Adding rrow did not help so will change from ia to ii.

Rnvtmodels$rrow[Rnvtmodels$Experiment%in%c("WMaA13BILO4")] <- 0
Rnvtmodels$resid[Rnvtmodels$Experiment %in% c("WMaA13BILO4")] <- "ii"

# Run ROUND 4 again:
tt <- summary(tmp.diag, vparameters = TRUE)$varcomp
tt[tt$component < 0,] # 0
tmp.diag$loglik # 9094.019
wald(tmp.diag)

xx <- summary(tmp.diag)$varcomp

#### Outlier Check ------------------------------------------------------------------------------------------------------------------------------------

aom.diag <- update(tmp.diag, aom = TRUE)
# Check for any outliers
out1 <- cbind(Rnvtdata, aom.diag$aom$R)
out1$residuals <- resid(aom.diag)
out1$fitted <- fitted(aom.diag)
subset(out1, abs(stdCondRes) > 4) # 2
subset(out1, abs(stdCondRes) > 4)[,c('Experiment','Range','Row','Rep','ColRep','RowRep','Variety','yield','stdCondRes')]

# Outlier Investigations
# Outlier 1 of 2
# Experiment  Range Row Rep ColRep RowRep       Variety    yield stdCondRes
# WMaA14NORT2     3  16   2      2      3 LRPB Spitfire 4.588745   4.070391

(tmp <- subset(out1, out1$Experiment=='WMaA14NORT2' & out1$Range=='3'))
hist(tmp$yield); tmp$yield
(tmp <- subset(out1, out1$Experiment=='WMaA14NORT2' & out1$Row=='16'))
hist(tmp$yield); tmp$yield
(tmp <- subset(out1, out1$Experiment=='WMaA14NORT2' & out1$Variety=='LRPB Spitfire'))
hist(tmp$yield); tmp$yield
tmp <- subset(out1, out1$Variety=='LRPB Spitfire')
hist(tmp$yield); tmp$yield

# Outlier 2 of 2
# Experiment  Range Row Rep ColRep RowRep       Variety    yield stdCondRes
# WMaA16BELL2     3  24   2      2      3        Qalbis 5.103571  -4.111392

(tmp <- subset(out1, out1$Experiment=='WMaA16BELL2' & out1$Range=='3'))
hist(tmp$yield); tmp$yield
(tmp <- subset(out1, out1$Experiment=='WMaA16BELL2' & out1$Row=='24'))
hist(tmp$yield); tmp$yield
(tmp <- subset(out1, out1$Experiment=='WMaA16BELL2' & out1$Variety=='Qalbis'))
hist(tmp$yield); tmp$yield
tmp <- subset(out1, out1$Variety=='Qalbis')
hist(tmp$yield); tmp$yield

# Would check with TSP but not conducting this analysis in real time.
# Things seem to be okay in regards to yield being typical of the variety, row and range
# for the experiment and the variety generally.

# QQNorm plots + scres heatmaps
Colresid.print(data=out1, outfile=paste0("1Graphics/", met.name, "-", date,"-Residplots-tol=3.75.pdf"))


#### Diag Model for Genetic Variance ------------------------------------------------------------------------------------------------------------------
# Now ready to run the diag model to get genetic variances out.

#Generate VarietyKeep and VarietyDrop, put VarietyDrop in sparse
Rnvtdata$VarietyDrop <- Rnvtdata$VarietyKeep <- Rnvtdata$Variety
Rnvtdata$VarietyKeep[substr(tolower(Rnvtdata$Variety), 1, 4)=="fill"] <- NA
Rnvtdata$VarietyDrop[substr(tolower(Rnvtdata$Variety), 1, 4)!="fill"] <- NA
Rnvtdata$VarietyKeep <- factor(Rnvtdata$VarietyKeep)
Rnvtdata$VarietyDrop <- factor(Rnvtdata$VarietyDrop)

unique(Rnvtdata[order(Rnvtdata$Variety), c("Variety", "VarietyKeep", "VarietyDrop")])

mt <- Colmodel.fit(models = Rnvtmodels, data = Rnvtdata)
asr.diag <- asreml(yield ~ Experiment +
                     at(Experiment, mt$animaldmg):animaldmg + 
                     at(Experiment, mt$est):est + 
                     at(Experiment, mt$earlygs):earlygs +
                     at(Experiment, mt$lrow):lin(Row) +
                     at(Experiment, mt$lcol):lin(Range),
                   sparse=~ VarietyDrop,
                   random = ~ diag(Experiment):VarietyKeep +
                     at(Experiment, mt$crep):ColRep + at(Experiment, mt$rrep):RowRep + 
                     at(Experiment, mt$rrow):Row + at(Experiment, mt$rcol):Range,
                   residual = ~ dsum(~ar1(Range):ar1(Row)| Experiment, levels = mt$resid$aa) +
                                dsum(~ar1(Range):id(Row)| Experiment, levels = mt$resid$ai) +
                                dsum(~id(Range):ar1(Row)| Experiment, levels = mt$resid$ia) +
                                dsum(~id(Range):id(Row)| Experiment, levels = mt$resid$ii), 
                   na.action = na.method(y='include', x='include'), data=Rnvtdata)
asr.diag <- update(asr.diag) #Converged in 13+13+2 its and 1.5+1.0+1.0 secs per iteration
asr.diag <- update(asr.diag, aom = TRUE)
asr.diag$loglik # 12567.71
save.image()

#save(list = ls(), file = "metfix.RData")
#load("metfix.RData")

#### Check Genetic Variances are all positive ---------------------------------------------------------------------------------------------------------
vc1 <- summary(asr.diag)$varcomp
unique(vc1[grep("Experiment:VarietyKeep", row.names(vc1)), "bound"]) # Bound
subset(vc1, bound =="B")
# Experiment:VarietyKeep!Experiment_WMaA13NIND4 2.393907e-07        NA      NA     B   0
# Experiment:VarietyKeep!Experiment_WMaA15COOA2 3.128517e-07        NA      NA     B   0
trial.info$SAGIInclude[trial.info$Experiment %in% c("WMaA13NIND4", "WMaA15COOA2")] <- FALSE
trial.info$SAGIcomment <- as.character(trial.info$SAGIcomment)
trial.info$SAGIcomment[trial.info$Experiment %in% c("WMaA13NIND4", "WMaA15COOA2")] <- "Genetic variance bound"

with(trial.info, table(METInclude))
with(trial.info, table(SAGIInclude))
with(trial.info, table(METInclude, SAGIInclude)) 

#### Extract output for breeders use ------------------------------------------------------------------------------------------------------------------
# Extract the trial mean yields, genetic variances, nvarieties and METinclude.
g0 <- vc1[grep("Experiment:VarietyKeep", row.names(vc1)),]
g0$Experiment <- factor(gsub("Experiment:VarietyKeep!Experiment_", "", row.names(g0)))
row.names(g0) <- NULL

g1 <- Rnvtdata %>% group_by(Experiment) %>% summarise(EMY = round(mean(yield, na.rm = TRUE), 2),
                                                     nV = n_distinct(Variety))

#### Check EMY ----------------------------------------------------------------------------------------------------------------------------------------
# Set very low EMY to SAGIinclude = FALSE...
hist(g1$EMY)
subset(g1, EMY < 1) # This cutoff of 1 will be different for each crop
# Don't remove any trials from BarleyNorth for this

# A tibble: 2 x 3
# Experiment    EMY    nV
# <fct>       <dbl> <int>
#WMaA14BULL2  0.38    34
#WMaA14MUNG4  0.99    32

# Comment:
# Are these 'frosted' trials? CHecked report - they are not. 
# BULL2 was removed for low trial mean yield but MUNG4 was not
# Will follow same protocol

 trial.info$SAGIInclude[trial.info$Experiment=="WMaA14BULL2"] <- FALSE
 trial.info$SAGIcomment <- as.character(trial.info$SAGIcomment)
 trial.info$SAGIcomment[trial.info$Experiment=="WMaA14BULL2"] <- "Trial Mean yield v low"
# trial.info$SAGIcomment <- "NONE"

 g2 <- merge(g1, g0[, c("Experiment", "component", "bound")], by = "Experiment")
 g3 <- merge(g2, trial.info[, c("Experiment", "Year", "Region","METInclude", "SAGIInclude")], by = "Experiment")
 names(g3)[4] <- c("GeneticVar")
 g3$GeneticVar <- round(g3$GeneticVar, 4)
 g3 <- g3[order(g3$METInclude, g3$SAGIInclude, g3$EMY),]
 
 
 p1 <- ggplot(g3, aes(x = EMY, y = GeneticVar, color = METInclude, label = substring(Experiment, 5, 11))) + 
   geom_text(size = 3, fontface = "bold") + xlab("Environment Mean Yield (t/ha)") + ylab("Genetic Variance") + 
   scale_colour_manual(values=c("red", "blue"))
 p2 <- ggplot(g3, aes(x = EMY, y = GeneticVar, color = factor(Year), label = Experiment)) + 
   geom_text(size = 2) + xlab("Environment Mean Yield (t/ha)") + ylab("Genetic Variance")
 p3 <- ggplot(g3, aes(x = EMY, y = GeneticVar, color = factor(Region), label = Experiment)) + 
   geom_text(size = 2) + xlab("Environment Mean Yield (t/ha)") + ylab("Genetic Variance")
 library(gridExtra)

pdf(paste0("1Graphics/",  met.name, "-", date.out, "-EMY-GenVar.pdf"))
grid.arrange(p1, p2, p3, nrow = 3); dev.off() #Send this file to the breeders too... ##usually includes , p3

#### Prepare Preliminary Data checking file -----------------------------------------------------------------------------------------------------------
# This file is created to send to breeders, this file must be checked 
# by breeders before MET analysis can continue

# Identify the trials that have METinclude = FALSE and SAGIinclude==FALSe
# Look at the comments online to identify why. Ideally we'd have these comments linked in so that we can see straight away.
(met.xclude <- droplevels(trial.info[trial.info$METInclude=="FALSE" | trial.info$SAGIInclude=="FALSE", 
                                     c("Experiment", "Year", "METInclude", "AOVComment","SAGIInclude", "SAGIcomment")]))
(met.xclude <- droplevels(subset(trial.info, METInclude=="FALSE" | SAGIInclude=="FALSE", 
                                 select = c(MET, Experiment, Year, State, Region, NearestTown, Latitude, Longitude,  
                                            SowingDate, HarvestDate, METInclude, AOVComment, SAGIInclude, SAGIcomment))))

Rnvtdata <- Rnvtdata[, grep("Include", names(Rnvtdata), invert = TRUE)] #Exclude columns that have include in them
Rnvtdata <- merge(Rnvtdata, trial.info[, c("Experiment", "METInclude","SAGIInclude")],
                  by = "Experiment")

# #Set up Info sheet from template and manual
# info.sheet <- read.csv("met-date-prelimDataCheck-InfoSheet.csv",
#                        header = FALSE, as.is = TRUE)
# info.sheet$V2[info.sheet$V1=="Title"] <- paste(met.name, date.out, "prelimDataCheck.xlsx", sep = "-")
# info.sheet$V2[info.sheet$V1=="Author"] <- "Mandy L Moore" #Your name, doh!
# info.sheet$V2[info.sheet$V1=="Date"] <- date.out
# 
# 
# fileXls <- paste0("1Results/", met.name, "-", date.out, "-prelimDataCheck.xlsx")
# unlink(fileXls, recursive = FALSE, force = FALSE)
# exc <- loadWorkbook('fileXls.xlsx', create = TRUE)
# createSheet(exc,'Info')
# writeWorksheet(exc, info.sheet, sheet = "Info", startRow = 1, startCol = 1, header = FALSE)
# createSheet(exc,'VarietybyYear')
# writeWorksheet(exc, vfreq, sheet = "VarietybyYear", startRow = 1, startCol = 1)
# createSheet(exc,'ExperimentCheck')
# writeWorksheet(exc, g3, sheet = "ExperimentCheck", startRow = 1, startCol = 1)
# createSheet(exc,'METexcludeTrials')
# writeWorksheet(exc, met.xclude, sheet = "METexcludeTrials", startRow = 1, startCol = 1)
# saveWorkbook(exc)

# SEND THIS FILE TO THE BREEDERS FOR LAST CHECKS BEFORE STARTING THE MET. OPEN IT AND CHECK THAT ALL IS IN ORDER FIRST!!!!

save(list = c("Rnvtdata", "Rnvtmodels"), file = paste0("1RData/", met.name, "-", date.out, "-data4METcleanedall.RData"))

#### Prepare the MET dataset for analysis ------------------------------------------------------------------------------------------------------------
met.xclude # all are frosted trials

# we're going to include these frosted trials in the Rolling MET
# Rnvtdata <- droplevels(Rnvtdata[Rnvtdata$METInclude==TRUE,])

Rnvtdata <- droplevels(Rnvtdata[Rnvtdata$SAGIInclude==TRUE,])
Rnvtmodels <- Rnvtmodels[is.element(Rnvtmodels$Experiment, Rnvtdata$Experiment),]

#Check to see that all trials are there
n_distinct(Rnvtdata$Experiment); n_distinct(Rnvtmodels$Experiment) #Both the same,keep going!

(ne <- nlevels(Rnvtdata$Experiment)) #134
(nv <- nlevels(Rnvtdata$VarietyKeep)) #129 (without fillers)

#Identify what varieties were lost when trials were exclued.
variety.init <- unique(Rnvtdata$Variety)
variety.final <- unique(Rnvtdata$VarietyKeep)

variety.init[!variety.init%in%variety.final]
#[1] FILLER_1  FILLER_10 Filler 


#### Final Checks--------------------------------------------------------------------------------------------------------------------------------------
#Check the number of experiments and Varieties are as you would expect. Original before subsetting
(met.sum.yr <- Rnvtdata %>% group_by(Year) %>% summarise(nExp = n_distinct(Experiment),
                                                           nVar = n_distinct(VarietyKeep),
                                                           nPlots = length(yield),
                                                           minYield = round(min(yield, na.rm = TRUE),2),
                                                           meanYield = round(mean(yield, na.rm = TRUE),2),
                                                           maxYield = round(max(yield, na.rm = TRUE),2)))
# Year   nExp  nVar nPlots      minYield meanYield maxYield
# <fct> <int> <int>  <int>         <dbl>     <dbl>    <dbl>
# 2013  27(1)    65   2835(90)      0.54      2.89     6.24
# 2014  26(1)    44   2586(102)     0.3       3.18     7.17
# 2015  27(1)    62   3480(144)     1.09      3.58     5.98
# 2016     29    65   4314          2.17      4.74     7.65
# 2017     28    53   3177          0.38      2.5      7.95


Rnvtdata <- merge(Rnvtdata, trial.info[, c("Experiment", "State", "Region")],
                 by.x = "Experiment", by.y = "Experiment", all.x = TRUE)
#Summary by Region
(met.sum.reg <- Rnvtdata %>% group_by(State, Region) %>% summarise(nExp = n_distinct(Experiment),
                                                                    nVar = n_distinct(VarietyKeep),
                                                                    nPlots = length(yield),
                                                                    minYield = round(min(yield, na.rm = TRUE),2),
                                                                    meanYield = round(mean(yield, na.rm = TRUE),2),
                                                                    maxYield = round(max(yield, na.rm = TRUE),2)))
# State Region  nExp  nVar nPlots minYield meanYield maxYield
# <fct> <fct>  <int> <int>  <int>    <dbl>     <dbl>    <dbl>
# NSW   N/E       30   104   4218     0.42      3.92     7.95
# NSW   N/W       35    94   4536     0.3       3.51     7.17
# QLD   CQ        22    69   2130     0.43      3.26     5.98
# QLD   SEQ        9   103   1134     1.04      4.25     7.65
# QLD   SWQ       38    73   4038     0.38      2.97     6.55


#### Variety concurrence by year ----------------------------------------------------------------------------------------------------------------------
# Check variety connectivity over the years
tt <- with(Rnvtdata[!is.na(Rnvtdata$yield),],table(Variety,Year))
tt[tt>1] <- 1
y0 <- t(tt)%*%tt; y0

#        Year
#   Year 2013 2014 2015 2016 2017
#   2013   64   31   34   28   21
#   2014   31   44   35   29   22
#   2015   34   35   62   47   31
#   2016   28   29   47   64   36
#   2017   21   22   31   36   52

#### Save the data and models file -------------------------------------------------------------------------------------------------------------------
save(list = ls(), file = paste0("1RData/",met.name, "-METfix.RData"))
#load(paste0(met.name, "-METsum.RData"))
save(list = c("Rnvtdata", "Rnvtmodels"), file = paste0("1RData/",met.name, "-", date.out, "-data4RMETcleanedready.RData"))  

# For KNVT - includes frosted trials
write.csv(Rnvtmodels, paste0("1ForKNVT/", met.name, "-", date.out, "-modelsFINAL.csv")) 


#### END OF SCRIPT - GO TO 03-METanalysis.R  ---------------------------------------------------------------------------------------------------------
file.edit('Scripts/03-METAnalysis-v2.R')
