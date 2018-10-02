###### Script Information ---------------------------------------------------------------------------------------------------------------------------
#
# Title: 04-METout_v2.R
# Author: Mandy Moore (Ky L Mathews, Daniel Tolhurst (mostly his, just adpated from co-located for wheat and Barley))
# 
# Date: September 2018
# Dataset: Wheat North
#
# Description: 
#   - Getting data from asr.final and obtaining BLUPS etc.
#
# Input: Data from 03-METAnalysis.R
# Siblings: & 04_Heatmaps_v2.R 
# Children: 05_METresults_v1.R  & Report.Rnw
#
# Modifications (made by KM to DT code):
#   - ny = ne cos Experiments, not YrLoc
#   - Removed ChemBlock code
#   - fixed grep of Comp for scores cos I had a variety called Compass....
#   - in allout.df the check.names cos I used Experiment, not YrLoc then I changed code to be more generic by replacing substr() with strsplit()
#
# Modifications (made by MM to KM code):


###### Data importation and package initialisation ----------------------------------------------------------------------------------------------------
rm(list = ls())

met.name <- "WheatMainNorth-frostid"

#Import data (if not run straight after 02-MetAnaysis.R)
load(paste0("1RData/", met.name, "-allmodels.RData"))
# load(paste0("1RData/", met.name, "-allMETinclude-METanalysisOUT.RData"))

#load(paste0("1RData/", met.name, "model.RData"))
load(paste0("1RData/", met.name, "-METanalysisOUT.RData"))

source("1Functions/NVTfns5.6.R")  #need this for model.fit
setup.fn(computer = 'linux') # If you are NOT using a mac the computer argument must be specified!
library(xtable)


#######################
# METout.R
#######################
# final model is rr5 + ar1
# Alter for the FA-k fitted
k <- 6                    #________ALTER HERE__________
final.asr <- models.asr[[paste0("rr",k, "AR1")]] #

(vnam <- sort(as.character(levels(Rnvtdata$VarietyKeep)))) #KLM: just doing this to be safe cos of ordering below...
(nv <- length(levels(Rnvtdata$Variety))) # 132 #Note: this does include fillers
(nvk <- length(levels(Rnvtdata$VarietyKeep))) # 129 #Note: this does not include fillers
(nfill <- length(levels(Rnvtdata$VarietyDrop))) # 3

(enam <- sort(as.character(levels(Rnvtdata$Experiment)))) #KLM: added the sort here cos for some reason this levels(nvtdata$Experiment) is not alphabetical... T
                                                       #This means that ASREML must do a factor() internallly cos the output is all alphabetical
(ne <- length(levels(Rnvtdata$Experiment))) # 134

#######################
# rr output parameters
#######################

lam <- matrix(summary(final.asr, vparameters=T)$vparameters[[paste("rr(Experiment, ",k,"):VarietyKeep",sep="")]],ncol=k+1)[,-1]
dimnames(lam) <- list(enam, paste("lam", 1:k, sep = "_"))
lam

psi <- matrix(summary(final.asr, vparameters=T)$vparameters[[c("Experiment:VarietyKeep")]],ncol=1)
dimnames(psi) <- list(enam,paste("psi"))
summary(as.vector(psi))
psi

Gmat <- lam%*%t(lam)+ diag(as.vector(psi))
summary(diag(Gmat))
dimnames(Gmat) <- list(row.names(lam),row.names(lam)) #KLM: addin in to make transparent
#frostid
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01900 0.05526 0.09253 0.11896 0.15912 0.59305  

CVEmat <- lam%*%t(lam)
dimnames(CVEmat) <- list(row.names(lam),row.names(lam)) #KLM: addin in to make transparent
summary(diag(CVEmat))
#frostid
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.008967 0.046198 0.077964 0.100248 0.137957 0.503708 

Cmat <- cov2cor(Gmat)
temp.cmat <- Cmat
diag(temp.cmat) <- NA
summary(as.vector(temp.cmat))
#frostid
#     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# -0.7511  0.1592  0.3854  0.3472  0.5674  0.9355     134  

### rotated and centred loadings -----------
lamss <- svd(lam)
lamstar <- -lam %*% lamss$v
dd <- 1/sqrt(diag(Gmat))
lamstarc <- diag(dd) %*% lamstar
dimnames(lamstarc) <- dimnames(lamstar) <- dimnames(lam)
# use these in cluster ordered heatmap

#######################
# % vaf
#######################

(vaf <- as.matrix(diag(CVEmat)/(diag(Gmat))*100))
hist(vaf)
sort(vaf)
summary(as.vector(vaf))
# frostid
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.59   77.52   88.15   84.90   94.63  100.00

100 * sum(diag(lam %*% t(lam)))/sum(diag(Gmat)) # 84.26813 overall -- mean of the sums
mean(vaf) #84.89626 -- sum of the means
vaf.df <- data.frame(Env = row.names(vaf), value = vaf)
p2 <- ggplot(vaf.df, aes(value)) + geom_histogram() + 
  xlab("Percentage variance accounted for") + 
  ylab("Frequency")
pdf(paste0("1Graphics/", met.name, "-sitevarexplhist-finalmodel.pdf")); p2;dev.off() #KLM: added this for 1Report


### range in %vaf -----------
summary(vaf)
length(vaf) # 134
sort(vaf[vaf<50,])
# frostid
# WMaA17KILC4 WMaA14SURA4 WMaA16COOL2 
# 34.59114    43.95614    49.33074 

diag(CVEmat)
ggplot(as.data.frame(Gmat), aes(diag(CVEmat), diag(Gmat), label = rownames(Gmat))) + 
  geom_text(size = 3) + labs(x="CVEvar",y="TGvar")

range(psi); hist(psi)
sort(psi[psi>.05,])  #KLM: choose a cutoff based on histogram to see which trials have large psi.
# WMaA16DULA4 WMaA17SPRS4 WMaA16COOA2 WMaA16WALG2 WMaA16ROMA4 WMaA14COOL2 WMaB17DUAR4 WMaA16COOL2 WMaA17KILC4 WMaA13COOA2 
#  0.05427643  0.06601716  0.07431258  0.08557551  0.08934172  0.09023603  0.09110564  0.11149099  0.11579189  0.14976981 

######################
# CVE-BLUPS
#######################

coef.df <- data.frame(summary(final.asr, coef=T)$coef.random)

var.df <- coef.df[grep('VarietyKeep', rownames(coef.df)),]
nrow(var.df) #35346
(ne*2+k)*nvk #35346 regblups, (2 = diag(Experiment) + rr(Experiment, k))

### scores (f) -----------
#scores.df <- var.df[grep('Comp', rownames(var.df)),] #KLM: This didn't work because I had a variety = Compass...
scores.df <- var.df[grep(paste(paste0('Comp', 1:k), collapse="|"), rownames(var.df)),]
nrow(scores.df) #774
nvk*k #774
(scores <- matrix(scores.df[,1],ncol=k))
dimnames(scores) <- list(vnam,paste("f",1:k,sep="_"))
scores

scorestar <- -scores %*% lamss$v
dimnames(scorestar) <- dimnames(scores)
variety.df <- as.data.frame(scorestar)
variety.df$Variety <- rownames(variety.df)
colnames(variety.df)[grep("f",colnames(variety.df))] <- paste("Score",1:k,sep="")
variety.df <- variety.df[,c("Variety", paste("Score",1:k,sep=""))]
variety.df$Variety <- factor(variety.df$Variety)
rownames(variety.df) <- NULL


### SVE-BLUPS (deltas) -----------  
sveBlup.df <- var.df[grep('rr(', rownames(var.df), invert=T, fixed = TRUE),] # i.e deltas
nrow(sveBlup.df) #17286
ne*nvk # 17286

### CVE-BLUPS (regBLUPS) -----------
cveBlup.df <- var.df[grep('rr(', rownames(var.df), fixed = TRUE),]
cveBlup.df <- cveBlup.df[grep(paste(paste0('Comp', 1:k), collapse="|"), rownames(cveBlup.df),invert=T),] 
nrow(cveBlup.df) #17286
ne*nvk # 17286

# check cgblups...
ggplot(cveBlup.df, aes(cveBlup.df$solution, as.vector(as.matrix(variety.df[,c(1+(1:k))])%*%t(lamstar)))) + 
  geom_point() + labs(x="CVE-BLUP",y="lam*f") #KLM: 1:1 ratio, nice check.
  

### Formatting CVE-Blups -----------
allout.df <- cveBlup.df
allout.df$pev <- allout.df$std.error^2 
names(allout.df) <- c("CVE","se","z","PEV")
allout.df <- allout.df[c("CVE","se","PEV")]
head(allout.df)
allout.df$Experiment <- factor(unlist(lapply(strsplit(row.names(allout.df), ":"), function(x) unlist(lapply(strsplit(x[1], "_"), function(x) x[2])))))
#allout.df$Variety <- factor(unlist(lapply(strsplit(row.names(allout.df), ":"), function(x) unlist(lapply(strsplit(x[2], "_"), function(x) x[2])))))
#Had to replace above row of code because a wheat variety had a ":" in it and cannot stsplit on that...
#Now splitting on "_" but this could cause problems in other crops.
allout.df$Variety <- factor(unlist(lapply(strsplit(row.names(allout.df), "_"), function(x) x[3])))#kk <- strsplit(row.names(allout.df), "_")
cve0 <- data.frame(Experiment = names(diag(CVEmat)), CVEvar=diag(CVEmat))
allout.df <- merge(allout.df, cve0, by = "Experiment")
names(allout.df)[names(allout.df)=="se"] <- "seCVE"
allout.df <- allout.df[,c("Experiment","Variety","CVE", "seCVE","PEV","CVEvar")]

save.image()

allout.df$VE <- cveBlup.df$solution + sveBlup.df$solution 
names(allout.df)
head(allout.df)
#######################
# present genos...
#######################

pres.df <- with(nvtdata[!is.na(nvtdata$yield),],table(VarietyKeep,Experiment))
dim(pres.df)
# 129, 126
sum(pres.df>0)/prod(dim(pres.df)) #[1] 0.3044789
pres <- as.vector(pres.df)

### 0, 1 or 2 ----------- THis is for checking dual tolerance varieties and for co-located trials
tt <- with(nvtdata[!is.na(nvtdata$yield),],table(Experiment,VarietyKeep))
tt[tt>1]<-1
tt<-data.frame(tt)
head(tt)
range(tt$Freq)
tt <- tt[order(tt$Experiment,tt$VarietyKeep),]
rownames(tt) <- NULL
head(tt)

#THis is for checking dual tolerance varieties and for co-located trials
ggplot(allout.df, aes(as.numeric(as.vector(pres>0)),as.numeric(as.vector(tt$Freq>0)))) + 
  geom_point() + labs(x="0-1 pres",y="0-1-2 pres") # 0 and 1 only
ggplot(allout.df, aes(as.numeric(as.vector(pres<1)),as.numeric(as.vector(tt$Freq<1)))) + 
  geom_point() + labs(x="0-1 pres",y="0-1-2 pres") # 0 and 1 only
head(tt)

allout.df$pres <- as.vector(tt$Freq)

#######################
# accuracies
#######################

# accuracy ie correlation between true and predicted
allout.df$ACC <- sqrt(1-allout.df$PEV/allout.df$CVEvar) # should be no error message, but if NA then look below.
#Check where ACC = NA
subset(allout.df, is.na(ACC)==TRUE)
subset(allout.df, ACC <0.10)
unique(subset(allout.df, is.na(ACC)==TRUE)$Variety) #only for Filler where no effects anyway - won't happen when using VarietyKeep
#The following were all only in one trial so to estimate the accuracy for some trials is difficult.
#This is partly due to the lam coming from 's+1' iterations of ASREML but the PEV (i.e. COEFF) coming from 's' iterations.
#Set the accuracies to 0. 
kk1  <- subset(allout.df, Experiment=="WMaA14CAPE4") #Not sure why acc for these checks are NA... 
kk1$Variety[is.na(kk1$ACC)==TRUE]
kk1  <- subset(allout.df, Experiment=="WMaA17KILC4") #Not sure why acc for these checks are NA... 
kk1$Variety[is.na(kk1$ACC)==TRUE]
sort(unique(subset(nvtdata, Experiment=="WMaA17KILC4")$Variety)) #None of the varieties with ACC=NA are in this trial

#diag(lamstar[rownames(lamstar)=="WMaA17KILC4",]%*%t(lamstar[rownames(lamstar)=="WMaA17KILC4",]))/diag(Gmat)["WMaA17KILC4"]*100
nvtdata[nvtdata$Variety %in% subset(allout.df, is.na(ACC)==TRUE)$Variety,]
Cmat[c("WMaA13BROO4","WMaA14CAPE4","WMaA17KILC4"),c("WMaA13BROO4","WMaA14CAPE4","WMaA17KILC4")]

#Shows no correlation with the site where these NA accuracies were grown and the two where the NA accuracies occur.

 unique(subset(allout.df, is.na(ACC)==TRUE)$pres)
 allout.df$ACC[is.na(allout.df$ACC)==TRUE] <- 0

# summary of accuracies ...
ACC <- with(allout.df,tapply(ACC,Experiment,mean,na.rm=T))
summary(ACC)

#nofrostid
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7528  0.8262  0.8589  0.8543  0.8852  0.9130 

head(allout.df)
allout.df <- allout.df[,c("Experiment","Variety","VE","CVE","seCVE","pres","ACC","CVEvar")] # leave PEV and CVEvar in for now

save.image()

#######################
# Trial details
#######################

allout.df$"%vaf" <- rep(vaf,each=nvk)

trials.df <- data.frame(Experiment=rownames(vaf),
                        EMY=as.vector(with(nvtdata,tapply(yield,Experiment,mean,na.rm=T))), #KLM: had to make this as.vector() so that when it's read to Excel it is read as numeric, not a character
                        "%vaf"=vaf,
                        RotatedLoads=lamstar) 
rownames(trials.df) <- NULL
names(trials.df)[grep("Load",names(trials.df))] <- paste("Load",1:k,sep="")
head(trials.df)
nrow(trials.df)

# Obtain Trial Level environment details - note, we'll have to update the SAGIinclude in this file and reimport into KNVT
met0 <- read.csv("1FromKNVT/Export_Wheat_Summary_20180103170359.csv")
met1 <- droplevels(met0[is.element(met0$TrialAcronym, enam)==TRUE,])
head(met1)
environment.df <- met1[,c("SiteYear","TrialAcronym","State","Region","Latitude","Longitude")] # also used for connectivity heatmap. This will be required from 30 Dec onwards.
names(environment.df)[grep("TrialAcronym",names(environment.df))] <- "Experiment"
names(environment.df)[grep("Year",names(environment.df))] <- "Year"
head(environment.df)
# check for duplicate entries (occurs when different regions are recorded for expts within a yrloc)
unique(with(environment.df,table(Experiment))) # 1 
# check for duplicate region entries (occurs when different regions are recorded for yrlocs within a location)
tt <- with(environment.df,table(substring(Experiment,4,11),Region))
tt[tt>1]<-1
unique(rowSums(tt)) # 1 - perfect!!

# check no Experiments are missing
trials.df$Experiment[!trials.df$Experiment %in% environment.df$Experiment] #length 0 - good
environment.df$YrLoc[!environment.df$Experiment %in% trials.df$Experiment] #NULL - good

summary.df <- merge(environment.df, unique(nvtdata[, c("Experiment", "METinclude", "SAGIinclude")]))
summary.df <- merge(summary.df,trials.df, by= "Experiment")
nrow(summary.df) 
ne # 126
summary.df <- summary.df[,c("Year","Experiment","State","Region","Latitude","Longitude","EMY", "X.vaf","METinclude",
                             "SAGIinclude", paste0("Load",1:k))]  #replaced paste(, sep="") with paste0 - a neat function!!!
head(summary.df)

#######################
# Results
#######################

allresults.df.temp <- merge(allout.df,summary.df, by="Experiment")
allresults.df <- merge(allresults.df.temp,variety.df, by="Variety")
allresults.df <- allresults.df[order(allresults.df$Experiment,allresults.df$Variety),]
names(allresults.df)[grep("vaf", names(allresults.df))] <- "%vaf"
allresults.df <- allresults.df[,c("Year","Experiment","Variety","VE","CVE","seCVE","pres","ACC", 
                                  "State","Region","Latitude","Longitude", "EMY", "METinclude",
                                  "SAGIinclude", "CVEvar",
                                  "%vaf",paste0("Load",1:k),paste0("Score",1:k))]
head(allresults.df)
nrow(allresults.df)
nvk*ne #16254
rownames(allresults.df) <- NULL

###############################
#Calculate overall and RMSD
###############################
#First check that the loadings from k = 1 are mostly positive.
#Check that the majority of the loadings for factor 1 are positive ,i.e.
# that first loading is a main effect - should be a fab effect with little cross-over at zero for Load1
ggplot(allresults.df, aes(Load1, Load1*Score1)) + geom_point()

allresults.df$OP <- mean(allresults.df$Load1)*allresults.df$Score1

allresults.df$dev <- allresults.df$CVE - allresults.df$Load1*allresults.df$Score1
tt <- sqrt(with(allresults.df,tapply(dev^2,Variety,mean)))
allresults.df$RMSD <- tt[as.character(allresults.df$Variety)]
variety.df <- merge(variety.df, unique(allresults.df[, c("Variety", "OP", "RMSD")]),
                    by = "Variety")
head(variety.df)
ggplot(variety.df, aes(RMSD, OP)) + geom_point()

             
# take what I want for results file  
results.df <- allresults.df[,c("Year","Experiment","Variety","VE","CVE","seCVE","pres","ACC", 
                               "State","Region","EMY","METinclude","CVEvar","%vaf",paste0("Load",1:k),paste0("Score",1:k))]
#round to 4 decimal places
for(i in 1:ncol(results.df)){
  if(is.numeric(results.df[,i])){
    results.df[,i] <- round(results.df[,i],4)
  }
}
head(results.df)

names(results.df) <- gsub("pres", "present", names(results.df))
names(results.df) <- gsub("ACC", "accuracy", names(results.df))

# # switch YrLoc for std trials to Experiment code -- NOT FOR WHEAT or BARLEY
# unique(results.df$YrLoc[results.df$YrLoc %in% std.df]) #17, YES!
# expt.df <- nvtdata[,c("YrLoc","Experiment")]
# expt.df <- expt.df[!duplicated(expt.df),]
# expt.df <- expt.df[expt.df$YrLoc %in% std.df,]
# nrow(expt.df) # 17
# 
# results.df$YrLoc <- as.character(results.df$YrLoc)
# unique(substring(results.df$YrLoc[results.df$YrLoc %in% std.df],2,8) == 
# substring(rep(expt.df$Experiment,each=nvk),5,11)) # TRUE
# results.df$YrLoc[results.df$YrLoc %in% std.df] <- as.character(rep(expt.df$Experiment,each=nvk))
# results.df$YrLoc <- factor(results.df$YrLoc, levels = unique(results.df$YrLoc))
# (ry <- as.character(levels(results.df$YrLoc)))

# some very quick checks
head(results.df)
head(coef.df[grep("rr(",rownames(coef.df), fixed = TRUE),]) #KLM: ensuring that we don't grab a variety with rr
for(i in 1:ne){#i <- 3
print(unique(vnam[!vnam %in% nvtdata$Variety[nvtdata$Experiment %in% enam[i]]] %in%
              results.df$Variety[results.df$present==0 & results.df$Experiment %in% enam[i]])) # KLM: had to change this when including all accuracy values
              #results.df$Variety[is.na(results.df$accuracy) & results.df$Experiment %in% enam[i]])) # check pres with data
  }
ggplot(results.df, aes(results.df$CVE, 
                       coef.df$solution[grep("rr(",rownames(coef.df), fixed = TRUE)][1:c(nvk*ne)])) + 
  geom_point() + labs(x="CVE",y="coef rrBlups")

#######################
# Output CSVs
#######################


csv1 <- results.df[,c("Year","Experiment","Variety", "CVE", "seCVE","present","accuracy")]
head(csv1)
nrow(csv1) #16254
ne*(nvk) # 16254

csv2 <- results.df[,c("Year","Experiment", "State","Region","METinclude","EMY","%vaf",paste0("Load",1:k))]
csv2 <- csv2[!duplicated(csv2),]
head(csv2)
nrow(csv2)
ne # 126

csv3 <- results.df[,c("Variety",paste0("Score",1:k))]
csv3 <- csv3[!duplicated(csv3),]
head(csv3)
nrow(csv3) #129
nv-nfill # 129

save(list = ls(), file = "nfrostidresults.RData")

#!!!!!!!!!!!!!!!!!!!!!!!!
#Need to run 04-METheatmaps.R first before the information is available for outputting
csv4 <- Cmat.ord.all
csv4 <- round(csv4,4)
nrow(csv4)
ne # 126`

csv5 <- Conmat.ord
nrow(csv5)
ne # 126


#########################
#write final data to file
#########################

#------------------------------
#Results for the app
#------------------------------CHANGE THE NAME FOR YOUR MET DATASET, THIS IS IMPORTANT FOR CROPS WITH MULTIPLE METS ---------------------

WheatMainNorthallresults.df <- allresults.df
names(WheatMainNorthallresults.df)[grep("seCVE", names(WheatMainNorthallresults.df))] <- "SE"
WheatMainNorthCSV4 <- Cmat 
diag(WheatMainNorthCSV4) <- NA

WheatMainNorthCSV5 <- Conmat 

WheatMainNorthlamstar <- lamstar
WheatMainNorthGmat <- Gmat 

save(list=c("nvtdata","WheatMainNorthallresults.df","WheatMainNorthCSV4","WheatMainNorthCSV5", "WheatMainNorthGmat", "WheatMainNorthlamstar"), 
     file= "1RData/WheatMainNorthAppResults13-17.RData")

#------------------------------
#Results for NVTdb
#------------------------------
ClusterNumber*	Cultivar*	BLUP*	PEV	Accuracy	VAF	SMY*	TrialCode*	IncludeInReport

NVTdb <- results.df[,c("Variety", "CVE", "seCVE", "accuracy", "%vaf", "EMY","Experiment")]
NVTdb <- cbind.data.frame(ClusterNumber=NA, NVTdb, IncludeInReport="Y")
NVTdb$PEV <- NVTdb$seCVE^2
names(NVTdb) <- c("ClusterNumber", "Cultivar",	"BLUP", "SE", "accuracy", "VAF","SMY", "TrialCode","IncludeInReport", "PEV")
NVTdb <- NVTdb[, c("ClusterNumber", "Cultivar",	"BLUP", "SE", "PEV", "accuracy", "VAF","SMY", "TrialCode","IncludeInReport")]
head(NVTdb)

write.csv(NVTdb, paste0("1results/WheatMainNorth-METresults-NVTdbimport-", Sys.Date(), ".csv"),row.names=F)

write.csv(results.df, paste0("1Report/xlsx/", met.name, "-results2017-", Sys.Date(), ".csv"),row.names=F)

save(list=ls(),file= paste0("1RData/1stageMETresults_", met.name, "13-17.RData"))
save(list = c("csv1","csv2","csv3","csv4","csv5","csv6","results.df","allresults.df","cveBlup.df","summary.df","trials.df","variety.df","environment.df",
              #"models.asr",
              "lam","psi","scores","CVEmat","Gmat","Cmat","lamstar","lamstarc","vaf",
              "nvtdata", "nvtmodels"), #"all.sv", "var.exclude", "mtdiag","drop.trials"), 
     file = paste0("1RData/1stageMETresults_", met.name, "13-17-finalMETresults.RData"))
save.image()


#Read in info sheet
#Set up Info sheet from template and manual
met.name.public <- paste(unlist(strsplit(met.name, "-"))[1], "2013-17")
info.sheet <- read.csv("~/Dropbox/NVT17/MET/0Resources/MET-results-InfoSheet-v1.csv", 
                       header = FALSE, as.is = TRUE)
info.sheet$V2[info.sheet$V1=="MET:"] <- paste(met.name.public,  ".xlsx", sep = "-")
info.sheet$V2[info.sheet$V1=="Author:"] <- "Ky L Mathews" #Your name, doh!
info.sheet$V2[info.sheet$V1=="Date:"] <- "10 Jan 2018"
info.sheet$V2[info.sheet$V1=="Years:"] <- "2013-17"
info.sheet$V3[info.sheet$V2=="Trials"] <- paste0("Summary of Environments in the MET (", ne, " in total).")
info.sheet$V3[info.sheet$V2=="Varieties"] <- paste0("Summary of Varieties in the MET (", nvk, " in total).")
info.sheet$V2[grep("NOTE:", info.sheet$V2)] <- paste0("NOTE: A Glossary of terms is provided with the accompanying report ",
                                                       met.name.public, " MET.pdf")


library(XLConnect)
fileXls <- paste0("1Report/xlsx/", met.name.public, "-METnvtresults.xlsx")
unlink(fileXls, recursive = FALSE, force = FALSE) 
exc <- loadWorkbook(fileXls, create = TRUE)
createSheet(exc,'Info')
writeWorksheet(exc, info.sheet, sheet = "Info", startRow = 1, startCol = 1, header = FALSE)
createSheet(exc,'CVEeffects')
writeWorksheet(exc, csv1, sheet = "CVEeffects", startRow = 1, startCol = 1)
createSheet(exc,'Trials')
writeWorksheet(exc, csv2, sheet = "Trials", startRow = 1, startCol = 1)
createSheet(exc,'Varieties')
writeWorksheet(exc, csv3, sheet = "Varieties", startRow = 1, startCol = 1)
createSheet(exc,'Correlation Matrix')
writeWorksheet(exc, csv4, sheet = "Correlation Matrix", startRow = 1, startCol = 1, rownames="Experiment")
createSheet(exc,'Connectivity Matrix')
writeWorksheet(exc, csv5, sheet = "Connectivity Matrix", startRow = 1, startCol = 1, rownames="Experiment")
saveWorkbook(exc)


#########################################################
#
#                  End of MET out Script
#
#               Move to 1 stage MET Heatmaps (but come back here for csv4 and csv5.)
#
#########################################################
