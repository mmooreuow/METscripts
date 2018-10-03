##################################################################################################
# DT's NVT Functions - VERSION 5.8, 22 Sep 2018
#
# Input from others as noted
#
##################################################################################################
setup.fn <- function(computer="MAC", APP=FALSE){
                     options(width = 150)
                     require(asreml)
                     require(ASExtras4)
                     require(cluster)
                     require(plyr)
                     require(dplyr)
                     require(ggplot2)
                     require(grid)
                     require(gridExtra)
                     require(gtools)
                     require(lattice)
                     require(mclust)
                     require(od)
                     require(reshape2)
                     if(tolower(computer) %in% c("mac","apple","apple mac")){dyn.load("/Library/Internet Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/lib/server/libjvm.dylib")}
                     require(XLConnect)
                     require(xtable)
                     if(APP){require(GGally)
                             require(ggmap)
                             require(leaflet)
                             require(magrittr)
                             require(mapdata)
                             require(markdown)
                             require(randomcoloR)
                             require(RColorBrewer)
                             require(shiny)
                             require(shinyBS)}
                } # end of setup function
########################################

##################################################################################################
# 
##################################################################################################
data.read <- function(data, infile, rows.skip=0, comments=FALSE, export=FALSE){
                      if(length(grep("xls", infile))==1){fileType <- "xls"}
                      if(length(grep("csv", infile))==1){fileType <- "csv"}
                      if(fileType %in% c("csv")){
                                          data.out <- read.csv(infile,na.strings = c("NA",""," ","-"), skip=rows.skip, strip.white=TRUE) 
                      }
                      if(fileType %in% c("xls")){
                                     require(XLConnect)
                                     require(data.table)
                                     data.out.temp <- loadWorkbook(infile)
                                     setMissingValue(data.out.temp, value = c("NA",""," ","-"))
                                     if(export){data.out.temp <- readWorksheet(data.out.temp,3,startRow=2)
                                                data.out.temp <- dcast(data.out.temp, TrialCode+Cultivar~Measurement, value.var="Value")
                                                names(data.out.temp) <- c("Experiment","Variety","ems","reps","predictedYield","se","weight")
                                                data <- droplevels(data[grep("buff|fill", tolower(data$Variety), invert=T),])
                                                # data <- droplevels(data[!is.na(data$yield),])
                                                EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                                                tt <- data.table::melt(with(data, tapply(yield,list(Variety,Experiment,get(paste(EnvID))),mean,na.rm=T)))
                                                colnames(tt) <- c("Variety","Experiment",EnvID,"rawYield")
                                                data.out <- merge(data.out.temp, tt, by= c("Experiment","Variety"), all = TRUE)
                                                data.out <- data.out[!is.element(data.out$rawYield,NA),]
                                                data.out$rawYield[is.na(data.out$rawYield)] <- NA
                                                yy <- levels(data[[EnvID]])
                                                ee <- levels(data$Experiment)
                                                vv <- levels(data$Variety)
                                                data.out[[EnvID]] <- factor(data.out[[EnvID]], levels=yy)
                                                data.out$Experiment <- factor(data.out$Experiment, levels=ee)
                                                data.out$Variety <- factor(data.out$Variety, levels=vv)
                                                data.out <- droplevels(data.out[order(data.out[[EnvID]], data.out$Experiment, data.out$Variety),])
                                                rownames(data.out) <- NULL
                                                data.out <- data.out[,c(EnvID,"Experiment","Variety","reps","rawYield","predictedYield","se","weight","ems")]
                                                tmp <- data.out[is.na(data.out$rawYield) & !is.na(data.out$predictedYield),]
                                                if(length(tmp$rawYield)){
                                                    message("You should NOT have a prediction if the raw yield is NA!")
                                                    rownames(tmp) <- NULL
                                                    print(tmp)
                                                }
                                     }
                                     if(!export & !comments){data.out <- readWorksheet(data.out.temp,1,startRow=6)}
                                     if(comments){data.out <- readWorksheet(data.out.temp,3,startRow=6)}
                      }
                      if(!export & !comments){data.out$Experiment<- as.character(data.out$Experiment)
                                          tt <- droplevels(subset(data.out[is.element(data.out$Harvest.Length, c(NA,0,""," ","-"))|
                                                                           is.element(data.out$Harvest.Width, c(NA,0,""," ","-"))|
                                                                           is.element(data.out$Kg.Plot, c(NA,0,""," ","-")),]))
                                          if(nrow(tt)>0){message("Warning, Harvest Attributes not Recorded for the Following Experiments")
                                                         sort(unique(print(unique(tt$Experiment))))}else{tt<- logical(0)}
                                          # for ASReml v4 we need a character prior to numeric for Environment
                                          data.out$Environment <- paste0(tolower(substring(as.character(data.out$Experiment),1,1)),
                                                                             substring(as.character(data.out$Experiment),5,11))
                                          data.out$Environment[substring(data.out$Experiment,1,1) %in% "C" & substring(data.out$Experiment,4,4) %in% c("B","C")] <- 
                                                  paste0(data.out$Environment[substring(data.out$Experiment,1,1) %in% "C" & substring(data.out$Experiment,4,4) %in% c("B","C")],"-",
                                                         substring(data.out$Experiment[substring(data.out$Experiment,1,1) %in% "C" & substring(data.out$Experiment,4,4) %in% c("B","C")],4,4))
                                          data.out$Year <- paste0('20',substring(as.character(data.out$Experiment),5,6))
                                          data.out$Location <- substring(as.character(data.out$Experiment),7,11)
                                          data.out$Series <- substring(as.character(data.out$Experiment),2,3)
                                          data.out$Series[!substring(data.out$Experiment,1,1) %in% c("C","W")] <- 
                                                   substring(as.character(data.out$Experiment),4,4)[!substring(data.out$Experiment,1,1) %in% c("C","W")]
                                          data.out$Series[substring(data.out$Experiment,1,1) %in% "B"] <- 
                                                   substring(as.character(data.out$Experiment),2,2)[substring(data.out$Experiment,1,1) %in% "B"]
                                          data.out$Series[substring(data.out$Experiment,1,1) %in% c("D","L")] <-
                                                   substring(as.character(data.out$Experiment),3,4)[substring(data.out$Experiment,1,1) %in% c("D","L")]
                                          data.out$State <- substring(as.character(data.out$Experiment),11,11)
                                          data.out <- data.out[order(data.out$Environment, as.numeric(data.out$Range), as.numeric(data.out$Row)),]
                                          data.out$Plots.orig <- paste0("c",sprintf("%03d", data.out$Range),"r", sprintf("%03d", data.out$Row))
                                          data.out$Blocks.orig <- paste0("cr",data.out$ColRep,"rr",data.out$RowRep)
                                          data.out$Rep.orig <- data.out$Rep
                                          
                                          nnames <- c("Environment", "Year", "State", "Location", "Experiment", 
                                                      "Series", "Range", "Row", "Rep", "ColRep", "RowRep",
                                                      "Variety", "Harvest.Length", "Harvest.Width", "Kg.Plot")
                                          covnames <- colnames(data.out)[!is.element(colnames(data.out), nnames)]
                                          # could clean up covariate names here
                                          data.out <- data.out[,c(nnames, covnames)]
                                          data.out <- data.out[order(as.character(data.out$Environment), as.character(data.out$Experiment), 
                                                                     as.numeric(data.out$Range), as.numeric(data.out$Row)),]
                                          rownames(data.out) <- NULL
                                          data.out <- list(data=data.out,nas=tt)
                       }
                       return(data.out)
              } # end of data read function
########################################

##################################################################################################
# 
##################################################################################################
layouts.read <- function(data, infile){
                         require(plyr)
                         EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                         if(length(grep("xls", infile))==1){fileType <- "xls"}
                         if(length(grep("csv", infile))==1){fileType <- "csv"}
                         if(fileType %in% c("csv")){
                                options(warn=-1)
                                temp <- read.csv(infile,na.strings = c("NA"," ","-"), strip.white=TRUE,header=TRUE)
                                options(warn=0)
                         }
                         if(fileType %in% c("xls")){
                                require(XLConnect)
                                temp <- loadWorkbook(infile)
                                setMissingValue(temp, value = c("NA"," ","-"))
                                temp <- readWorksheet(temp,1)
                         }
                         if(nrow(temp)<1){message("Warnings, No Layouts in File")
                                          outdata <- data[order(data[[EnvID]]),]
                                          rownames(outdata) <- NULL
                                          temp2 <- logical(0)
                                          ready <- data.frame(e1=NA,e2=NA,e3=NA,e4=NA)}
                         if(nrow(temp)>0){
                              temp$Environment <- paste0(tolower(substring(temp$e1,1,1)),substring(temp$e1,5,11))
                              temp$Year <- paste('20',substring(temp$e1,5,6),sep='')
                              temp$State <-substring(temp$Environment,8,8)
                              nnames <- c("Environment", "Year", "State", "Location", 
                                          "Series", "rows", "cols", "e1", "e2", "e3",
                                          "e4", "g1", "g2", "g3","g4")
                              temp <- temp[,c(nnames)]
                              temp <- temp[order(temp$Year,temp$Series,temp$State,temp$Location,temp$e1),]
                              otherdata <- droplevels(data[!data[[EnvID]] %in% temp$Environment,])
                              data <- droplevels(data[data[[EnvID]] %in% temp$Environment,])
                              temp2 <- droplevels(temp[temp$Environment %in% data[[EnvID]],])
                              rownames(temp2) <- NULL
                              ready <- temp2[temp2$e1 %in% data$Experiment & is.element(temp2$e2, c(data$Experiment,NA,""," ")) &
                                            is.element(temp2$e3, c(data$Experiment,NA,""," ")) & is.element(temp2$e4, c(data$Experiment,NA,""," ")),]
                              not.ready <- as.vector(as.matrix(temp2[!temp2[[EnvID]] %in% ready[[EnvID]],c("e1","e2","e3","e4"),]))
                              not.ready <- not.ready[!is.element(not.ready, c(NA, "", " "))]
                              if(nrow(temp2)<1){temp2 <- logical(0)}
                              if(length(not.ready)>0){message("Experiments not Ready for Analysis")
                                                      temp.notready <- not.ready[not.ready %in% data$Experiment]
                                                      temp.notready <- temp.notready[order(paste0(substring(temp.notready,1,1),substring(temp.notready,5,11)))]
                                                      print(temp.notready)
                                                      message("(missing these remaining co-located experiments)")
                                                      temp.wait <- not.ready[!not.ready %in% data$Experiment]
                                                      temp.wait <- temp.wait[order(paste0(substring(temp.wait,1,1),substring(temp.wait,5,11)))]
                                                      print(temp.wait)
                                                      data <- droplevels(data[!data$Experiment %in% not.ready,])
                                                      temp2 <- droplevels(temp2[temp2$Environment %in% data$Environment,])}
                              outdata <- rbind.fill(otherdata,data)
                              outdata <- outdata[order(outdata[[EnvID]]),]
                              rownames(outdata) <- NULL
                        }
                        layouts <- list()
                        layouts$Layouts <- temp2    
                        layouts$Experiments <- as.vector(as.matrix(ready[,c("e1","e2","e3","e4")]))
                        layouts$Experiments <- layouts$Experiments[!is.na(layouts$Experiments)]
                        layouts$Data <- outdata
                        layouts$not.ready <- logical(0)
                        layouts$waiting.on <- logical(0)
                        if(length(not.ready)>0){
                              layouts$not.ready <- temp.notready
                              layouts$waiting.on <- temp.wait
                        }
                        return(layouts)
                } # end of layout read function
########################################

##################################################################################################
# 
##################################################################################################
# DT with input from BC
# loop over Year locations
# if gap is null then data frame is as is
# assume Range x Row layout
# Expts are laid out in a super array
# ny x nx
# process each with each ny level
# as we may get ragged arrays
# denoted by NA in the set2.loc.layout matrix

layouts.build <- function(data, layouts){
                          require(plyr)
                          require(dplyr)
                          EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                          if(length(layouts)<1){temp.df11 <- data
                                                                Environment.df <- data[data$Experiment == "DT",]}
                          if(length(layouts)>1){
                             temp.df11 <- data[!data$Experiment %in% unique(c(as.character(layouts$e1),as.character(layouts$e2),
                                               as.character(layouts$e3),as.character(layouts$e4))),]
                             data <- data[data$Experiment %in% c(as.character(layouts$e1),as.character(layouts$e2),
                                          as.character(layouts$e3),as.character(layouts$e4)),]
                             temp <- split(layouts,layouts[[EnvID]])
                             temp <- lapply(temp,function(x){rr <- as.numeric(x[6])
                                                            cc <- as.numeric(x[7])
                                                            if(rr==1){xx <- sapply(x[,8:11],as.character)
                                                            xx <- matrix(xx[!is.na(xx)],byrow=T,nrow=1,ncol=cc)}
                                                            else{xx <- sapply(x[,8:11],as.character)
                                                            xx <- matrix(xx,byrow=T,ncol=cc,nrow=rr)}
                                                            xx})
                             Environment.layout <- temp
                             temp <- split(layouts,layouts[[EnvID]])
                             temp <- lapply(temp,function(x){rr <- as.numeric(x[6])
                                                            cc <- as.numeric(x[7])
                                                            if(rr==1){xx <- sapply(x[,12:15],as.numeric)
                                                            xx <- 2.5*matrix(xx[!is.na(xx)],byrow=T,nrow=1,ncol=cc-1)}
                                                            if(rr>1 & cc>1){xx <- sapply(x[,12:15],as.numeric)
                                                            xx <- 2.5*matrix(xx,byrow=T,ncol=cc-1,nrow=rr)}
                                                            if(rr>1 & cc==1){xx <- sapply(x[,12:15],as.numeric)
                                                            xx <- 2.5*matrix(xx,byrow=T,ncol=cc,nrow=rr)}
                                                            xx})
                             Environment.gaps <- temp
                             data$Variety <- as.character(data$Variety)
                             Environment.lst <- list()
                             yy <- levels(factor(data[[EnvID]]))
                             for(yloc in yy){gap <- Environment.gaps[[yloc]]
                                     if(length(gap)>0){this.df <- data[data[[EnvID]]==yloc,]
                                               this.df$Experiment <- as.character(this.df$Experiment)
                                               Nranges <- with(this.df,tapply(Range,Experiment,max))
                                               Nrows <- with(this.df,tapply(Row,Experiment,max))
                                               this.layout <- Environment.layout[[yloc]]
                                               Environment.rows <- with(this.df,tapply(Row,as.character(Experiment),max))
                                               Environment.rows <- Environment.rows[as.character(this.layout)]
                                               Environment.rows <- matrix(Environment.rows,ncol=ncol(this.layout),nrow=nrow(this.layout))
                                               layout.rows <- apply(Environment.rows,1,sum,na.rm=T) + apply(gap/2.5,1,sum,na.rm=T)
                                               row.lst <- list()
                                               for(i in 1:dim(this.layout)[1]){ne <- this.layout[i,]
                                                        ne <- ne[!is.na(ne)]
                                                        row.df <- subset(this.df,is.element(Experiment,ne))
                                                        row.df$Experiment <- factor(as.character(row.df$Experiment),levels=ne)
                                                        row.df$Variety <- as.character(row.df$Variety)
                                                        if(length(ne)==1){insert.lst <- list()
                                                                  endbit <- sum(Environment.rows[i,],na.rm=T) +
                                                                  sum(round(gap[i,]/2.5),na.rm=T)
                                                                  endbit <- max(layout.rows,na.rm=T) - endbit
                                                                  if(endbit>0){tt <- subset(row.df,Experiment==ne[1] & Row==1)
                                                                             tt$Variety <- 'Buffer'
                                                                             tt$Kg.Plot <- NA
                                                                             tt$Plots.orig <- NA
                                                                             tt$Blocks.orig <- NA
                                                                             tt$Rep.orig <- NA
                                                                             tt <- lapply(tt,function(x,gg){rep(x,each=gg)},gg=endbit)
                                                                             tt$Row <- rep(1:endbit,times=length(unique(row.df$Range)))
                                                                             tt$Row <- tt$Row + sum(Environment.rows[i,1:1],na.rm=T)
                                                                             insert.lst[[1]] <- tt
                                                                             insert.df <- as.data.frame(insert.lst[[1]])
                                                                             row.df <- rbind(row.df,insert.df)}
                                                                  row.lst[[i]] <- row.df}
                                                        else{insert.lst <- list()
                                                        for(j in 1:(length(ne)-1)){tt <- subset(row.df,Experiment==ne[j] & Row==1)
                                                                  tt$Variety <- 'Buffer'
                                                                  tt$Kg.Plot <- NA
                                                                  tt$Plots.orig <- NA
                                                                  tt$Blocks.orig <- NA
                                                                  tt$Rep.orig <- NA
                                                                  tt <- lapply(tt,function(x,gg){rep(x,each=gg)},gg=round(gap[i,j]/2.5))
                                                                  tt$Row <- rep(1:round(gap[i,j]/2.5),times=length(unique(row.df$Range)))
                                                                  if(j==1){tt$Row <- tt$Row + sum(Environment.rows[i,1:j])}
                                                                  else{tt$Row <- tt$Row + sum(Environment.rows[i,1:j]) + sum(round(gap[i,1:(j-1)]/2.5))}
                                                                  insert.lst[[j]] <- tt}
                                                        endbit <- sum(Environment.rows[i,]) + sum(round(gap[i,]/2.5))
                                                        endbit <- max(layout.rows) - endbit
                                                        if(endbit>0){j <- length(ne)
                                                                  tt <- subset(row.df,Experiment==ne[j] & Row==1)
                                                                  tt$Variety <- 'Buffer'
                                                                  tt$Kg.Plot <- NA
                                                                  tt$Plots.orig <- NA
                                                                  tt$Blocks.orig <- NA
                                                                  tt$Rep.orig <- NA
                                                                  tt <- lapply(tt,function(x,gg){rep(x,each=gg)},gg=endbit)
                                                                  tt$Row <- rep(1:endbit,times=length(unique(row.df$Range)))
                                                                  tt$Row <- tt$Row + sum(Environment.rows[i,1:j]) + sum(round(gap[i,1:(j-1)]/2.5))
                                                                  insert.lst[[j]] <- tt}
                                                        insert.df <- as.data.frame(insert.lst[[1]])
                                                        if(length(insert.lst)>1){for(j in 2:length(insert.lst)){
                                                                                          insert.df <- rbind(insert.df,as.data.frame(insert.lst[[j]]))}}
                                                        temp.df <- subset(row.df,Experiment==ne[1])
                                                        for(j in 2:(length(ne))){tt <- subset(row.df,Experiment==ne[j])
                                                                  tt$Row <- tt$Row + sum(Environment.rows[i,1:(j-1)]) + sum(round(gap[i,1:(j-1)]/2.5))
                                                                  temp.df <- rbind(temp.df,tt)}
                                                        row.df <- temp.df
                                                        row.df <- rbind(row.df,insert.df)
                                                        row.lst[[i]] <- row.df}}
                                               Environment.df <- row.lst[[1]]
                                               nrange <- max(Environment.df$Range)
                                               if(dim(this.layout)[[1]]>1){for(i in 2:length(row.lst)){
                                                                                           row.lst[[i]]$Range <- row.lst[[i]]$Range + nrange*(i-1)
                                                                                           Environment.df <- rbind(Environment.df,row.lst[[i]])}}
                                               Environment.lst[[yloc]] <- Environment.df}
                             else{Environment.lst[[yloc]] <- data[data[[EnvID]]==yloc,]}}
                             Environment.df <- Environment.lst[[1]]
                             Environment.df <- Environment.df[order(Environment.df$Range,Environment.df$Row),]
                             if(length(yy)>1){for(k in 2:length(yy)){Environment.df <- rbind(Environment.df,Environment.lst[[k]])}}
  
                             Environment.df <- Environment.df[order(as.character(Environment.df[[EnvID]]), 
                                                                    as.numeric(as.character(Environment.df$Range)),as.numeric(as.character(Environment.df$Row))),]
                          }
                             temp.df11 <- temp.df11[order(as.character(temp.df11[[EnvID]]), as.character(temp.df11$Experiment)),]
                             Environment.df2 <- rbind.fill(Environment.df,temp.df11)
                             Environment.df2[[EnvID]] <- as.character(Environment.df2[[EnvID]])
                             Environment.df2$Experiment <- as.character(Environment.df2$Experiment)
                             Environment.df2 <- Environment.df2[order(Environment.df2[[EnvID]]),]
                             
                             Environment.df2$Experiment <- factor(Environment.df2$Experiment, levels = unique(Environment.df2$Experiment))
                             Environment.df2[[EnvID]] <- factor(Environment.df2[[EnvID]])
                             
                             Environment.df2$Year <- factor(as.numeric(Environment.df2$Year))
                             Environment.df2$YrCon <- factor(Environment.df2$YrCon)
                             Environment.df2$Location <- factor(Environment.df2$Location)
                             Environment.df2$Series <- factor(Environment.df2$Series)
                             Environment.df2$State <- factor(as.numeric(Environment.df2$State))
                             Environment.df2$Range <- factor(as.numeric(Environment.df2$Range))
                             Environment.df2$Row <- factor(as.numeric(Environment.df2$Row))
                             Environment.df2$Rep <- factor(as.numeric(Environment.df2$Rep))
                             Environment.df2$ColRep <- factor(as.numeric(Environment.df2$ColRep))
                             Environment.df2$RowRep <- factor(as.numeric(Environment.df2$RowRep))
                             Environment.df2$Variety <- factor(Environment.df2$Variety)
                             Environment.df2$Plots.orig <- factor(Environment.df2$Plots.orig)
                             Environment.df2$Blocks.orig <- factor(Environment.df2$Blocks.orig)
                             Environment.df2$Rep.orig <- factor(Environment.df2$Rep.orig)
                             
                             # finally, add in yield
                             Environment.df2$yield <- Environment.df2$Kg.Plot/Environment.df2$Harvest.Length/Environment.df2$Harvest.Width*10
                             Environment.df2$yield.orig <- Environment.df2$yield
                             nnames <- c(EnvID, "Year", "State", "Location", "Experiment", 
                                         "Series", "Range", "Row", "Rep", "ColRep", "RowRep",
                                         "Variety", "Harvest.Length", "Harvest.Width", "Kg.Plot", "yield")
                             onames <- c("YrCon","Plots.orig","Blocks.orig","Rep.orig","yield.orig")
                             covnames <- colnames(Environment.df2)[!is.element(colnames(Environment.df2),c(nnames,onames))]
                             # NA Fillers/buffers
                             var.exclude <- sort(unique(as.character(Environment.df2$Variety[grep("fill|buff",tolower(Environment.df2$Variety))])))
                             if(length(covnames)>0){Environment.df2[Environment.df2$Variety %in% var.exclude,c("yield",covnames)] <- NA}
                             Environment.df2$RowRep[Environment.df2$Variety == "Buffer"] <- NA
                             # could clean up covariate names here
                             Environment.df2 <- Environment.df2[,c(nnames, covnames, onames)]
                             rownames(Environment.df2) <- NULL
                             M <- list()
                             M$data <- Environment.df2
                             M$var.exclude <- var.exclude
                             return(M) 
                     } # end of co-located layouts function
########################################


##################################################################################################
# 
##################################################################################################
trials.sum <- function(data){
  EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
  yy <- levels(factor(data[[EnvID]]))
  temp <- data.frame(State=factor(tapply(as.character(data$State),data[[EnvID]],unique)[yy]),
                     Environment=yy,
                     Year=factor(tapply(as.character(data$Year),data[[EnvID]],unique)[yy]),
                     Location=factor(tapply(as.character(data$Location),data[[EnvID]],function(x) unique(as.character(x)))[yy]),
                     Nexpts=as.numeric(tapply(as.character(data$Series),data[[EnvID]],function(x) length(unique(as.character(x))))[yy]),
                     Series=tapply(as.character(data$Series),data[[EnvID]],function(x) unique(as.character(x))),
                     Nrange=as.numeric(with(droplevels(subset(data, Variety != "Buffer")), tapply(as.numeric(as.character(Range)),get(paste(EnvID)),function(x) length(unique(x))))[yy]),
                     Nrow=as.numeric(with(droplevels(subset(data, Variety != "Buffer")), tapply(as.numeric(as.character(Row)),get(paste(EnvID)),function(x) length(unique(x))))[yy]))
  temp <- temp[order(substring(as.character(temp$Environment),1,1),
                     temp$Environment,temp$Location,temp$Series),]
  rownames(temp) <-NULL
  print(temp)
  return(temp)
} # end of trial summary function
########################################

##################################################################################################
# 
##################################################################################################
reps.check <- function(data, BlockID= "Rep.orig", outfile1, outfile2, ExptID = "Trial", VarietyID = "Genotype", ColID = "Col", RowID = "Row"){
                       require(ggplot2)
                       if(length(grep("Variety$",colnames(data)))==0){colnames(data)[grep(paste0(VarietyID,"$"),colnames(data))] <- "Variety"}
                       if(length(grep("Experiment$",colnames(data)))==0){colnames(data)[grep(paste0(ExptID,"$"),colnames(data))] <- "Experiment"}
                       if(length(grep("Range$",colnames(data)))==0){colnames(data)[grep(paste0(ColID,"$"),colnames(data))] <- "Range"}
                       if(length(grep("Row$",colnames(data)))==0){colnames(data)[grep(paste0(RowID,"$"),colnames(data))] <- "Row"}
                       data <- droplevels(data[data$Variety!="Buffer",])
                       rr <- list()
                       ee <- levels(data$Experiment)
                       for(i in ee){data.temp <- droplevels(subset(data,Experiment==i))
                                    rr[[i]] <- table(data.temp$Variety,data.temp[[BlockID]])} 
                       rr.1 <- lapply(rr, function(x) table(apply(x, 2, function(x) length(unique(x)))))
                       rep.bad1 <- names(rr.1)[sapply(rr.1,length) > 1 | lapply(sapply(rr.1,names),max) >1]
                       
                       rro <- list()
                       for(i in ee){data.temp <- droplevels(subset(data,Experiment==i))
                                    rro[[i]] <- tapply(data.temp[[BlockID]], list(data.temp$Row,data.temp$Range), unique)}
                       rro <- lapply(rro, function(x) table(apply(x, 2, function(x) length(unique(x)))))
                       rep.bad2 <- names(rro)[sapply(rro,length)!=1]
                                    
                       rrt <- tapply(data[[BlockID]], data$Experiment, function(x) max(as.numeric(x)))
                       rep.bad3 <- names(rrt[rrt>3 & !is.na(rrt)])
                                    
                       rry <- tapply(data[[BlockID]], data$Experiment, function(x) sum(is.na(x))/length(x))
                       rep.bad4 <- names(rry[rry>0.1 & rry <1]) # more than 10% misssing values.
                       rep.bad5 <- names(rry[rry ==1]) # 100% missing values
                       rep.na <- logical(0)
                       if(length(rep.bad5)>0){rep.na <- sort(unique(rep.bad5))}

                       rep.bad <- sort(unique(c(rep.bad1, rep.bad2, rep.bad3, rep.bad4)))
                       if(length(rep.bad)<1){message("No Figures Produced (no fixes required)")}
                       if(length(rep.bad)>0){
                                   a.0 <- droplevels(data[data$Experiment %in% rep.bad,]) # only spit out pdf of those with issues
                                   a.0$plot.grp <- ceiling(as.numeric(a.0$Experiment)/9)
                                   np <- length(unique(a.0$plot.grp))
                                   pdf(outfile1)
                                   for(i in 1:np){#i <- 1
                                             p.data <- droplevels(subset(a.0, plot.grp==i))
                                             p.data$Row <- factor(p.data$Row, levels = rev(levels(p.data$Row)))
                                             p1 <- ggplot(p.data, aes_string(x = "Range", y = "Row", fill = BlockID)) + geom_tile() + 
                                                          theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0), 
                                                                axis.text.y = element_text(size=6, angle=30), title = element_text(size=9)) + # change here if many rows 
                                                          facet_wrap(~Experiment, scales = "free")
                                             options(warn=-1)       
                                             print(p1)
                                             options(warn=0)}
                                   dev.off()
                                            
                                   # check that reps resolvable wrt varieties...
                                   a.01 <- data.frame(table(a.0$Variety,a.0[[BlockID]],a.0$Experiment))
                                   colnames(a.01) <- c("Variety",BlockID,"Experiment","Freq")
                                   a.01 <- a.01[a.01$Freq > 0,]
                                   a.01$plot.grp <- ceiling(as.numeric(a.01$Experiment)/6)
                                   np2 <- length(unique(a.01$plot.grp))
                                   pdf(outfile2)
                                   for(j in 1:np2){#i <- 1
                                             p.data <- droplevels(subset(a.01, plot.grp==j))
                                             p1 <- ggplot(p.data, aes_string(x =BlockID, y = "Variety", fill = "factor(Freq)")) + geom_tile() + labs(fill="Freq") +
                                                          theme(strip.text.x = element_text(size = 10, colour = "black", angle = 0),
                                                                axis.text.y = element_text(size=6, angle=30), title = element_text(size=9)) + # change here if many varieites 
                                                          facet_wrap(~Experiment, scales = "free", nrow = 2, ncol = 3) 
                                             options(warn=-1)
                                             print(p1)
                                             options(warn=0)}
                                   dev.off()
                                   if(i != np | j != np2){message("Warning, Function Error, Please Re-run Rep Check Function")}
                                   if(i == np & j == np2){message("Output Figures Printed to File")}
                       }
                       M <- list()
                       M$rep.na <- rep.na
                       M$rep.bad <- rep.bad
                       return(M)
              }# end of rep checks function
########################################

##################################################################################################
# 
##################################################################################################
layouts.print <- function(data, outfile){
                          require(ggplot2)
                          EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                          tt <- with(droplevels(data[data$YrCon != as.character(data$Experiment),]),table(get(paste(EnvID)),Experiment))
                          tt[tt>1] <- 1
                          colocated <- rownames(tt)[rowSums(tt)>1]
                          if(length(colocated)<1){message("No Co-located Experiments in Data")}
                          if(length(colocated)>0){
                                     a.02 <- droplevels(data[data[[EnvID]] %in% colocated,])
                                     yy <- levels(a.02[[EnvID]])
                                     np <- nlevels(a.02[[EnvID]])
                                     pdf(outfile)
                                     for(i in 1:np){#i <- 1
                                                p.data <- droplevels(subset(a.02, Environment==yy[[i]]))
                                                p.data$Row <- factor(p.data$Row, levels = rev(levels(p.data$Row)))
                                                p.data$ExperimentFill <- as.character(p.data$Experiment)
                                                p.data$ExperimentFill[p.data$Variety =="Buffer"] <- "Buffer"
                                                p.data$ExperimentFill <- factor(p.data$ExperimentFill)
                                                expts <- cols <- levels(p.data$ExperimentFill)
                                                cols[cols == "Buffer"] <- "grey"
                                                p.data$Crop <- substring(p.data[[EnvID]],1,1)
                                                # colnames(p.data)[grep(EnvID, colnames(p.data))] <- "Environment"
                                                crop <- unique(p.data$Crop)
                                                if(crop %in% c("c")){
                                                          cols[substring(cols,3,3)=="C"] <- "#F8766D"
                                                          cols[substring(cols,3,3)=="I"] <- "#7CAE00"
                                                          cols[substring(cols,3,3)=="R"] <- "#00BFC4"
                                                          cols[substring(cols,3,3)=="T"] <- "#C77CFF"}
                                                if(crop %in% c("d","l")){
                                                          cols[substring(cols,3,4)%in%c("DA","AA")] <- "#F8766D"
                                                          cols[substring(cols,3,4)%in%c("DB","AB")] <- "#7CAE00"
                                                          cols[substring(cols,3,4)%in%c("KA","NA")] <- "#00BFC4"
                                                          cols[substring(cols,3,4)%in%c("KB","NB")] <- "#C77CFF"}
                                                if(!crop %in% c("c","d","l")){
                                                          cols[substring(cols,4,4) == "A"] <- "#F8766D"
                                                          cols[substring(cols,4,4) == "B"] <- "#7CAE00"
                                                          cols[substring(cols,4,4) == "C"] <- "#00BFC4"
                                                          cols[substring(cols,4,4) == "D"] <- "#C77CFF"}
                                                p1 <- ggplot(p.data, aes(x = Range, y = Row, fill = ExperimentFill)) + geom_tile()+ 
                                                             geom_text(aes(label=Variety), size=2.3) + 
                                                             scale_fill_manual("Experiment",breaks = expts[expts!="Buffer"],values=cols) +
                                                             theme(legend.position = "top",
                                                                   strip.text.x = element_text(size = 10, colour = "black", angle = 0), 
                                                                   axis.text.y = element_text(size=6, angle=30), title = element_text(size=9))
                                                options(warn=-1)
                                                print(p1)
                                                options(warn=0)}
                                     dev.off()
                                     if(i != np){message("Warning, Function Error, Please Re-run Colocated Graphics Function")}
                                     if(i==np){message("Output Figures Printed to File")}
                           }
                 }# end of ColocatedGraphics function
########################################

##################################################################################################
# 
##################################################################################################
covariates.print <- function(data, infile, outfile, ExptID = "Trial", VarietyID = "Genotype", ColID = "Col"){
                             analID <- "NVT"
                             if(length(grep("Variety$",colnames(data)))==0){colnames(data)[grep(paste0(VarietyID,"$"),colnames(data))] <- "Variety";analID <- "PBA"}
                             if(length(grep("Experiment$",colnames(data)))==0){colnames(data)[grep(paste0(ExptID,"$"),colnames(data))] <- "Experiment";analID <- "PBA"}
                             if(length(grep("Range$",colnames(data)))==0){colnames(data)[grep(paste0(ColID,"$"),colnames(data))] <- "Range";analID <- "PBA"}
                             yield.ord <- grep("yield$",colnames(data))
                             yrcon.ord <- min(grep("YrCon|Plots.orig",colnames(data)))
                             cov.df <- logical(0)
                             if(!yrcon.ord>yield.ord+1){message("No Covariates in Data")}
                             if(yrcon.ord>yield.ord+1){
                                      if(length(grep("xls", infile))==1){fileType <- "xls"}
                                      if(length(grep("csv", infile))==1){fileType <- "csv"}
                                      if(fileType %in% c("csv")){
                                              options(warn=-1)
                                              cov.data <- read.csv(infile,na.strings = c("NA"," ","-"), strip.white=TRUE)[,1:2]
                                              options(warn=0)
                                      }
                                      if(fileType %in% c("xls")){
                                              require(XLConnect)
                                              cov.data <- loadWorkbook(infile)
                                              setMissingValue(cov.data, value = c("NA"," ","-"))
                                              cov.data <- readWorksheet(cov.data,1)
                                      }
                                      cov.names <- colnames(data)[(yield.ord+1):(yrcon.ord-1)]
                                      index1 <- match(cov.names, cov.data[,1])
                                      if(sum(is.na(index1))>1){message("Covariates not in File")
                                              print(cov.names[is.na(index1)])}
                                      index2 <- match(cov.names, colnames(data))[!is.na(index1)]
                                      index1 <- index1[!is.na(index1)]
                                      colnames(data)[index2] <- as.character(cov.data[index1,2])
                                      cov.names <- colnames(data)[(yield.ord+1):(yrcon.ord-1)]

                                      # list of YrLoc, Experiment x Co-variate
                                      cov.df1 <- list()
                                      EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                                      for(i in 1:length(cov.names)){#i<-1
                                              # only print out those with > 1 unique value
                                              # one unique value- NA
                                              tt <- tapply(data[[cov.names[i]]], data[[EnvID]], function(x) length(unique(x[!is.na(x)])))
                                              tt2 <- unique(as.character(data$Experiment[!is.na(data[cov.names[i]])]))
                                              exp.df <- unique(as.character(data$Experiment[data[[EnvID]] %in% names(tt[tt>1]) & data$Experiment %in% tt2]))
                                              if(length(exp.df)>0){cov.df1[[i]] <- data.frame(cov=cov.names[i],Experiment=exp.df)}
                                      }
                                      cov.df <- do.call(rbind,cov.df1)
                                      if(length(cov.df)==0){message("No Suitable Covariates in Data (only one unique value)")
                                                            cov.df <- logical(0)
                                        }else{cov.df$Experiment <- factor(cov.df$Experiment, levels=levels(data$Experiment))
                                              cov.df <- cov.df[order(cov.df$cov,cov.df$Experiment),]
                                              cov.df$Experiment <- as.character(cov.df$Experiment)
                                              cov.df$cov <- as.character(cov.df$cov)
                                              cov.names <- sort(unique(cov.df$cov))

                                              # print to pdf
                                              if(length(cov.names)==0){message("No Suitable Covariates in Data (only one unique value)")
                                                               cov.df <- logical(0)}
                                              if(length(cov.names)>0){
                                                    pdf(outfile,width=6,height=6)
                                                    np <- cov.names
                                                    for(i in np){#i <- "animaldmg"
                                                          cov.exp <- cov.df$Experiment[cov.df$cov==i]
                                                          p.data <- droplevels(subset(data, Experiment%in%cov.exp))
                                                          p.data <- p.data[order(p.data[[EnvID]], p.data$Experiment),]
                                                          cov.env <- unique(as.character(p.data[[EnvID]]))
                                                      
                                                          for(j in cov.env){# j<-"w17TULO2"
                                                               p.data.temp2 <- p.data.temp <- droplevels(data[data[[EnvID]] %in% j,])
                                                               p.data.temp <- droplevels(p.data.temp[!is.na(p.data.temp[[i]]),])
                                                               col.df <- 1 < length(unique(p.data.temp$Experiment))
                                                               con.df <- 1 < length(unique(p.data.temp$YrCon == as.character(p.data.temp$Experiment)))
                                                               expts <- cols <- levels(p.data.temp$Experiment)
                                                               p.data.temp <- p.data.temp[order(p.data.temp$Experiment),]
                                                               p.data.temp$Crop <- substring(p.data.temp[[EnvID]],1,1)
                                                               crop <- unique(p.data.temp$Crop)
                                                               if(crop %in% c("c")){
                                                                   cols[substring(cols,3,3)=="C"] <- "#F8766D"
                                                                   cols[substring(cols,3,3)=="I"] <- "#7CAE00"
                                                                   cols[substring(cols,3,3)=="R"] <- "#00BFC4"
                                                                   cols[substring(cols,3,3)=="T"] <- "#C77CFF"}
                                                               if(crop %in% c("d","l")){
                                                                   cols[substring(cols,3,4)%in%c("DA","AA")] <- "#F8766D"
                                                                   cols[substring(cols,3,4)%in%c("DB","AB")] <- "#7CAE00"
                                                                   cols[substring(cols,3,4)%in%c("KA","NA")] <- "#00BFC4"
                                                                   cols[substring(cols,3,4)%in%c("KB","NB")] <- "#C77CFF"}
                                                               if(analID != "NVT"){
                                                                 cols[substring(cols,2,2)==1] <- "#F8766D"
                                                                 cols[substring(cols,2,2)==2] <- "#7CAE00"
                                                                 cols[substring(cols,2,2)==3] <- "#00BFC4"
                                                                 cols[substring(cols,2,2)==4] <- "#C77CFF"}
                                                               if(!col.df){
                                                                   pp <- ggplot(p.data.temp,aes_string(x=i,y="yield")) + 
                                                                                ggtitle(paste(unique(as.character(p.data.temp$Experiment))))
                                                               }
                                                               if(!col.df & crop %in% c("c","d","l") & analID=="NVT"){
                                                                   pp <- ggplot(p.data.temp,aes_string(x=i,y="yield", colour="Experiment")) + 
                                                                                ggtitle(paste(unique(as.character(p.data.temp$Experiment))))+
                                                                                scale_colour_manual("Experiment",breaks = expts,values=cols)
                                                               }
                                                               if(col.df){
                                                                   pp <- ggplot(p.data.temp,aes_string(x=i, y="yield", colour="Experiment")) +
                                                                                ggtitle(paste(unique(as.character(p.data.temp[[EnvID]]))))+
                                                                                scale_colour_manual("Experiment",breaks = expts,values=cols)
                                                               }
                                                      
                                                               pp <- pp + geom_point() + xlab(i) +
                                                                     theme(legend.position="top", strip.text.x = element_text(size = 10, colour = "black", angle = 0),
                                                                     axis.text.y = element_text(size=6), title = element_text(size=9))
                                                               options(warn=-1)
                                                               print(pp)
                                                               options(warn=0)
                                                      
                                                               # covariate heatmap plot
                                                               p.data.temp2$Row <- factor(p.data.temp2$Row, levels= rev(sort(unique(as.numeric(as.character(p.data.temp2$Row))))))
                                                               p.data.temp$Row <- factor(p.data.temp$Row, levels= rev(sort(unique(as.numeric(as.character(p.data.temp$Row))))))
                                                               if(!col.df&!con.df){
                                                                  pp2 <- ggplot(p.data.temp,aes_string(x="Range",y="Row", fill=i)) + 
                                                                         ggtitle(paste(unique(as.character(p.data.temp$Experiment)))) + geom_tile() +
                                                                         scale_fill_gradient(low = "red", high = "green")
                                                                  options(warn=-1)
                                                                  print(pp2)
                                                                  options(warn=0)
                                                               }
                                                               if(col.df&!con.df){
                                                                  pp2 <- ggplot(p.data.temp2,aes_string(x="Range",y="Row", fill=i)) +
                                                                         ggtitle(paste(unique(as.character(p.data.temp2[[EnvID]])))) + geom_tile() +
                                                                         scale_fill_gradient(low = "red", high = "green")
                                                                  options(warn=-1)
                                                                  print(pp2)
                                                                  options(warn=0)
                                                               }
                                                               if(con.df){
                                                                  for(l in unique(p.data.temp2$Experiment)){
                                                                        pp2 <- ggplot(p.data[p.data.temp2$Experiment == l,],aes_string(x="Range",y="Row", fill=i)) +
                                                                               ggtitle(l) + geom_tile() +
                                                                               scale_fill_gradient(low = "red", high = "green")
                                                                        options(warn=-1)
                                                                        print(pp2)
                                                                        options(warn=0)
                                                                  }
                                                               }
                                                          }
                                                    }
                                                    if(i == last(np)){message("Output Figures Printed to File")}
                                                    if(i != last(np)){message("Warning, GGplot Error, Please Re-run Covariates Function")}
                                                    dev.off() 
                                              }
                                          }        
                             }
                             M <- list()
                             M$cov.df <- cov.df
                             if(length(cov.df>0)){print(cov.df)}
                             if(analID == "PBA"){colnames(data)[grep("Variety$",colnames(data))] <- VarietyID
                                                 colnames(data)[grep("Experiment$",colnames(data))] <- ExptID
                                                 colnames(data)[grep("Range$",colnames(data))] <- ColID
                             }
                             M$data <- data
                             return(M)
                     }# end of covariates graphic function
########################################

######################################
equal.phiv3 <- function(obj=BLasrA.sv, data=rice.df, ConstraintID = "YrCon", 
                        constraint = "vcc", ExptID = "Trial", RowID='Row!cor', ColID='Range!cor', ResID = "R$") {
                        if(length(grep("Experiment$",colnames(data)))==0){colnames(data)[grep(paste0(ExptID,"$"),colnames(data))] <- "Experiment"}
                        EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                        envs <- levels(data[[EnvID]])
                        temp <- obj$vparameters.table
                        temp <- temp[grep(ConstraintID,temp$Component),]
                        order <- temp$Component
                        sigmas <- temp[grep(ResID,temp$Component),]
                        colphis <- temp[grep(ColID,temp$Component),]
                        rowphis <- temp[grep(RowID,temp$Component),]
                        temp.ss <- rbind(sigmas,colphis,rowphis)
                        temp.ss$YrLoc <- sapply(strsplit(as.character(temp.ss$Component),split='_|!'), function(x) x[2])
                        temp.ss$phi <- sapply(strsplit(as.character(temp.ss$Component),split='!'), function(x) x[2])
                        yrloc.temp <- temp.ss$YrLoc[temp.ss$YrLoc %in% data$Experiment]
                        temp.df <- data[data$Experiment %in% yrloc.temp, c(EnvID,"Experiment")]
                        temp.df <- temp.df[!duplicated(temp.df),]
                        temp.ss <- merge(temp.ss, temp.df, by.x="YrLoc", by.y="Experiment", all.x = T)
                        rownames(temp.ss) <- temp.ss$Component
                        temp.ss <- temp.ss[order,]
                        temp.ss$YrLoc <- as.character(temp.ss$YrLoc)
                        temp.ss$YrLoc[!is.na(temp.ss$Environment)] <- trimws(temp.ss$Environment[!is.na(temp.ss$Environment)])
                        temp.ss$Env.Phi <- paste(temp.ss$YrLoc,temp.ss$phi, sep="!")
                        temp.ss <- temp.ss[,grep("YrLoc|phi|Environment",colnames(temp.ss),invert=T),]  
                        temp.ss$Value <- with(temp.ss,tapply(Value,Env.Phi,mean))[temp.ss$Env.Phi]
                        rownames(temp.ss) <- NULL
  
                        M <- vcm.lm( ~ Env.Phi, data=temp.ss)
                        attr(M,'assign') <- NULL; attr(M,'contrasts') <- NULL
                        if(constraint == "vcm"){return(list(gammas=temp.ss,M=M))}
                        if(constraint == "vcc"){
                               M2 <- data.frame(Constraint= temp.ss[["Env.Phi"]], Scale=1)
                               rownames(M2) <- temp.ss[["Component"]]
                               M2$Constraint <- as.numeric(factor(M2$Constraint, levels=unique(M2$Constraint)))
                               M2 <- as.matrix(M2)
                               return(list(gammas=temp.ss,M=M2))
                        }
                        }# end of equal.phi function
########################################

##################################################################################################
#model.fit()
#model function: to determine which terms (design or covariate) will be at which Experiments.
#Note the residual is broken up into as many levels as there are in models$resid
#Uses require(reshape2) written by DT for co-located trials with input from KM for original fn

# fit linear term but not random term.... Hey Hey Hey
# fit linear term to 3 ranges/row.... Hey Hey Hey
Colmodel.fit <- function(models, data){
                         residID <- names(models)[grep("res",tolower(names(models)))]
                         EnvID <- colnames(models)[grep("Env|YrLoc",colnames(models))]
                         models <- models[!duplicated(models),]
                         tt <- table(models[[EnvID]], models$Experiment)
                         tt[tt>1] <- 1
                         tt <- rowSums(tt)
                         if(length(grep("YrCon",colnames(models)))>0){standard <- droplevels(models[models[[EnvID]] %in% names(tt[tt==1]) & models$Experiment != as.character(models$YrCon),c(EnvID,"Experiment")])
                              constrained <- droplevels(models[models[[EnvID]] %in% names(tt[tt>1]) & models$Experiment == as.character(models$YrCon),c(EnvID,"Experiment")])
                              colocated <- droplevels(models[models[[EnvID]] %in% names(tt[tt>1]) & models$Experiment != as.character(models$YrCon),c(EnvID,"Experiment")])
                         }else{constrained <- data.frame()
                               colocated <- data.frame()
                               standard <- levels(models[[EnvID]])}
                         cnames <- colnames(models)
                         if("Fmodel" %in% cnames){models <- models[,-grep("%ch|X.ch",cnames)];models <- models[,-grep("Fmodel",cnames)];models <- models[,-grep("Best",cnames)]}
                         require(reshape2)
                         indicator <- "YrCon"
                         if(length(grep("Q|B|W|O", substring(data$Experiment,1,1))) > 0 &
                            nrow(rbind(colocated,constrained))<1){indicator <- "Experiment"}
                         standard <- models[,c(EnvID, "Experiment")]
                         options(warn=-1)
                         model.long <- melt(models, id.vars = c(unique(c(indicator,EnvID,"Experiment"))), 
                                            measure.vars = names(models)[!is.element(names(models), c(EnvID,"Experiment", residID))])
                         options(warn=0)
                         model.long$YrConVar <- paste(model.long[[indicator]],model.long$variable,sep="_")
                         variables.table <- as.data.frame(tapply(model.long[[indicator]],with(model.long,list(variable,value)),unique))
                         
                         tt <- with(models, table(Experiment))
                         names(tt[tt>1])
                         if(length(names(tt[tt>1]))>0){message("Warning, Experiments that have more than one row in models.csv, Please Fix!")
                                                       print(names(tt[tt>1]))}
                         # costrained trials below
  #                           if(errors >0){model.short <- model.resolve[rowSums(model.resolve,na.rm=T)>0,]
  #                                         model.short$YrCon <- rownames(model.short)
  #                           YrCon.err <- as.character(rownames(model.short))
  #   model.melt <- subset(melt(model.short,id.vars = "YrCon"),value>0)
  #   model.melt$value <- NULL
  #   model.melt$YrConVar <- paste(model.melt$YrCon,model.melt$variable,sep="_")
  #   model.long2 <- droplevels(subset(model.long, YrConVar %in% model.melt$YrConVar))
  #   # what variable should now be modelled at the experiment level?
  #   # - only those that have values for one experiment
  #   v1.df <- data.frame(YrConVar = model.long2$YrConVar[model.long2$value==1], vExpt = model.long2$Experiment[model.long2$value==1]) # variable score Expts
  #   v2.df <- v1.df[duplicated(v1.df$YrConVar),]
  #   v1.df <- v1.df[!duplicated(v1.df$YrConVar),]
  #   v3.df <- v2.df[duplicated(v2.df$YrConVar),]
  #   v2.df <- v2.df[!duplicated(v2.df$YrConVar),]
  #   all.v1 <- merge(v1.df,v2.df,by = "YrConVar",all=T)
  #   all.v <- merge(all.v1,v3.df,by = "YrConVar",all=T)
  #   #### and for the expts without scores
  #   n1.df <- data.frame(YrConVar = model.long2$YrConVar[model.long2$value==0], nExpt = model.long2$Experiment[model.long2$value==0]) # variable score Expts
  #   n2.df <- n1.df[duplicated(n1.df$YrConVar),]
  #   n1.df <- n1.df[!duplicated(n1.df$YrConVar),]
  #   n3.df <- n2.df[duplicated(n2.df$YrConVar),]
  #   n2.df <- n2.df[!duplicated(n2.df$YrConVar),]
  #   all.n1 <- merge(n1.df,n2.df,by = "YrConVar",all=T)
  #   all.n <- merge(all.n1,n3.df,by = "YrConVar",all=T)
  #   all.df <- merge(all.v,all.n,by = "YrConVar",all=T)
  #   colnames(all.df)[2:ncol(all.df)] <- c(paste("vExpt",1:3,sep=""),paste("nExpt",1:3,sep=""))
  #   all.df$YrCon <- strsplit(as.character(all.df$YrConVar),"_")
  #   all.df$variable <- unlist(lapply(all.df$YrCon, function (x) x[[2]]))
  #   all.df$YrCon <- unlist(lapply(all.df$YrCon, function (x) x[[1]]))
  #   all.df <- all.df[,c(c(ncol(all.df)-1),ncol(all.df),1:c(ncol(all.df)-2))]
  #   # those to go via at(Expt)
  #   all.df$ExptLevel <- TRUE
  #   all.df$ExptLevel[!is.na(all.df$vExpt2)] <- FALSE # those to go via at(YrCon), or
  #   all.df$ExptLevel[all.df$variable %in% c("rep","crep","rrep","lrow","lcol","rrow","rcol")] <- FALSE # design terms not to go to Expt level.
  #   # print out which EnvID and Experiments dont have values for variable
  #   message("Warning, Experiments that do NOT have Variable Scores at a YrCon (i.e. nExpt), while others do (vExpt)")
  #   print(all.df[,c(1:2,4:ncol(all.df))]) 
  #   model.long <- model.long[,-grep("Experiment",colnames(model.long))]
  #   model.long <- model.long[!duplicated(model.long),]  
  #   model.long$value[model.long$YrConVar %in% all.df$YrConVar[all.df$ExptLevel == TRUE]] <- 0 #Expt level
  #   model.long$value[model.long$YrConVar %in% all.df$YrConVar[all.df$ExptLevel == FALSE]] <- 1 #YrCon level
  #   model.long <- model.long[!duplicated(model.long),] 
  #   #######
  #   model.long3 <- model.long2[model.long2$YrConVar %in% all.df$YrConVar[all.df$ExptLevel == TRUE],]
  #   model.ls2 <- split(model.long3, f = model.long3$variable)
  #   model.terms2 <- lapply(model.ls2, function(x) as.character(x$Experiment[x$value==1]))
  #   #Identify mterms that are to be used in this model
  #   mterms.keep2 <- names(model.terms2)[unlist(lapply(model.terms2, function(x) length(x)>0))==TRUE]
  #   mt2 <- model.terms2[mterms.keep2]
  #   #Set up the residual terms -- done at YrCon
  # }
                         # if(errors <1){model.long <- model.long[,-grep("Experiment",colnames(model.long))]
                                       # model.long <- model.long[!duplicated(model.long),] 
                                       # mt2<-logical(0)}
                         #############################
                         # colocated + standard trials
                         if(nrow(colocated)>0|nrow(standard)>0){
                              col.df <- model.long[model.long$Experiment %in% colocated$Experiment|model.long$Experiment %in% standard$Experiment,]
                              
                              # check for errors- only co-located trials
                              if(nrow(colocated)>0){
                                    model.resolve <- as.data.frame(with(model.long,tapply(value,list(YrCon,variable),function(x) length(unique(x)))))
                                    model.resolve[model.resolve<2] <- 0
                                    model.resolve[model.resolve>1] <- 1
                                    yield.ord <- grep("yield$",colnames(data))
                                    yrcon.ord <- grep("YrCon",colnames(data))
                                    cov.names <- c("rep","rrep")
                                    colnames(model.resolve)[grep("rep|rrep|lrow|lcol",colnames(model.resolve))] <- c("Rep","RowRep","Row","Range")
                                    if(yrcon.ord>yield.ord+1){cov.names <- c("Rep","RowRep","Row","Range",colnames(data)[(yield.ord+1):(yrcon.ord-1)])
                                            for(i in cov.names){# i <- "est"
                                                    temp.data <- list()
                                                    env.temp <- rownames(model.resolve)[model.resolve[[i]]==1]
                                                    for(j in env.temp){# j <- env.temp
                                                            temp.data <- droplevels(data[data[[EnvID]] == j,])
                                                            tt <- table(temp.data$Experiment,temp.data[[i]])
                                                            if(length(tt[rowSums(tt)<1])>0){
                                                                    model.resolve[rownames(model.resolve) == j, i] <- 0}}
                                    }}
                                    colnames(model.resolve)[grep("Rep|RowRep|Row|Range",colnames(model.resolve))] <- c("rep","rrep","lrow","lcol")
                                    errors <- length(model.resolve[model.resolve>0])
                                    if(errors>0){# print out which EnvID and Experiments dont have values for variable
                                                 message("Warning, Environments that have different models according to Experiment")
                                                 print(model.resolve[rowSums(model.resolve)>0,])}
                              }
                              #Set up the fixed/random terms
                              ls <- split(model.long, f = model.long$variable)
                              terms <- lapply(ls, function(x) as.character(x[[EnvID]][x$value==1]))
                              if(indicator == "Experiment"){terms <- lapply(ls, function(x) as.character(x$Experiment[x$value==1]))
                                                            if(length(terms)>6){terms <- terms[c(7:length(terms),3:4,1:2,5:6)]
                                                                }else{terms <- terms[c(3:4,1:2,5:6)]}}
                              if(indicator != "Experiment"){if(length(terms)>6){covnames <- c(sort(names(terms)[c(8:length(terms))]),names(terms)[1:7])
                                                                }else{covnames <- c(names(terms)[1:7])}
                                                            terms <- terms[covnames]}
                              terms.keep <- names(terms)[unlist(lapply(terms, function(x) length(x)>0))==TRUE]
                              tt <- terms[terms.keep]
                              tt <- lapply(tt, function(x) x[!duplicated(x)])
                              mt2 <- logical(0)
                         }
                         #Set up the residual terms
                         mt.resid <- split(models[,c(indicator, residID)], f = factor(models[[residID]]))
                         mt.resid$aa <- mt.resid$aa[!duplicated(mt.resid$aa),]
                         mt.resid$ai <- mt.resid$ai[!duplicated(mt.resid$ai),]
                         mt.resid$ia <- mt.resid$ia[!duplicated(mt.resid$ia),]
                         mt.resid$ii <- mt.resid$ii[!duplicated(mt.resid$ii),]
                         yrcon.df <- list()
                         yrcon.df$resid <- lapply(mt.resid, function(x) as.character(x[[indicator]]))
                         mt <- list(YrLoc = tt, YrCon = yrcon.df)
                         names(mt)[1] <- EnvID
                         if(indicator=="Experiment"){mt <- list(Experiment = c(tt,yrcon.df)); mt <- mt$Experiment}
                         print(mt)
                         return(mt)  
                } # end of co-located model.fit
########################################

##################################################################################################
# 
##################################################################################################
Colvario.print <- function(obj, data, outfile, EnvID = "Environment", ExptID = "Experiment", ColID='Range', RowID='Row', tol = 10){
  # my function for compositing variograms for co-located trials...
  require(ggplot2)
  EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
  env.df <- unique(alldata[,c(EnvID,ExptID)])
  rownames(env.df) <- NULL
  vario.data <- varioGram(obj)
  vario.data$gamma2 <- vario.data$gamma*vario.data$np
  pdf(outfile)
  for(i in levels(env.df$Environment)){#i<- "AN15"
    temp.vario.data <- droplevels(vario.data[vario.data$vv.groups %in% trimws(env.df[env.df[[EnvID]]==i,ExptID]),])
    cc <- sort(unique(temp.vario.data[[ColID]]))
    rr <- unique(temp.vario.data[[RowID]])
    temp.vario.data <- data.frame(Range= rep(cc,each=length(rr)),
                                  Row=rr,
                                  gamma=as.vector(t(tapply(temp.vario.data$gamma2, list(temp.vario.data[[ColID]],temp.vario.data[[RowID]]), function(x) sum(x)))),
                                  np= as.vector(t(tapply(temp.vario.data$np, list(temp.vario.data[[ColID]],temp.vario.data[[RowID]]), function(x) sum(x)))))
    temp.vario.data$gamma <- temp.vario.data$gamma/temp.vario.data$np
    temp.vario.data <- temp.vario.data[order(temp.vario.data$Row,temp.vario.data$Range),]
    rownames(temp.vario.data) <- NULL
    pp <- wireframe(gamma ~ Row * Range, data = droplevels(temp.vario.data[temp.vario.data$np>tol,]), 
                    drape = T, colorkey = F, zoom = 0.8, 
                    par.settings = list(axis.line = list(col = 'transparent'),
                                        layout.widths=list(left.padding=2,right.padding=-7),
                                        layout.heights=list(top.padding=-12,bottom.padding=-12)),
                    xlab = list(label = "Row displacement", cex = 1, rot=20), 
                    ylab = list(label = "Column displacement", cex = 1, rot=310),
                    zlab = list(label = "Semi variance", cex = 1, rot=94),
                    screen = list(z = 30, x = -60, y = 0), aspect = c(1, 0.66), 
                    scales = list(distance = c(1, 1, 1.2), tck = c(1.4, 1, 1), arrows = F, cex = 0.8, col = "black"))
    print(pp)
    grid::grid.text(i, x=unit(0.5, "npc"), y=unit(0.92, "npc"), gp = gpar(fontface="bold", fontsize=20))
  }
  dev.off()
  message("Output Figures Printed to File")
} # end of co-located variogram.print
########################################

##################################################################################################
# 
##################################################################################################
Colresid.print <- function(data, outfile, ExptID = "Trial", resID='residuals', scresID='stdCondRes', tol=3.5, ColID='Range', RowID='Row'){
                           require(ggplot2)
                           require(gridExtra)
                           # creates qqnorm plots for residuals
                           # Requires data-frame with indexing factors, residuals, studentised conditional residuals;
                           # tol=tolerance for "large" studentised conditional residuals
                           # loop over YrCons (i.e. residual)
                           if(length(grep("Experiment$",colnames(data)))==0){colnames(data)[grep(paste0(ExptID,"$"),colnames(data))] <- "Experiment"}
                           EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                           tt <- with(droplevels(data[data$YrCon == as.character(data$Experiment),]),table(get(paste(EnvID)),Experiment))
                           tt[tt>1] <- 1
                           colocated <- rownames(tt)[rowSums(tt)>1]
                           resname <- "Experiment"
                           if(length(colocated)>0){resname <- "YrCon"}
                           yc <- levels(data[[resname]])
                           pdf(outfile,width=6,height=6)
                           for(i in 1:length(yc)){#i<-1
                                   res.df <- droplevels(data[data[[resname]] == yc[i],])
                                   res.df[[resID]][abs(res.df[[resID]]) < 1e-9] <- NA
                                   res.df <- res.df[order(res.df[[ColID]],res.df[[RowID]]),]
                                   qq <- qqnorm(res.df[[resID]],plot.it=F)
                                   res.df$sampquant <- qq$y
                                   res.df$theoquant <- qq$x
                                   y <- quantile(res.df[!is.na(res.df[[resID]]),resID], c(0.25, 0.75))
                                   x <- qnorm(c(0.25, 0.75))
                                   res.df$slope <- diff(y)/diff(x)
                                   res.df$int <- y[1L] - res.df$slope * x[1L]
                                   res.df$big <- paste('<',tol,sep=' ')
                                   res.df$big[abs(res.df[[scresID]])>tol] <- paste('>',tol,sep=' ')
                                   if(!is.null(RowID)){res.df$lab <- paste(res.df[[ColID]],res.df[[RowID]],sep=":")}
                                   if(is.null(RowID)){res.df$lab <- res.df[[ColID]]}
                                   plot.title <- yc[i]
                                   if(!is.null(RowID)){plot.subtitle <- paste0('Residual effects indexed as ',ColID,':',RowID)}
                                   if(is.null(RowID)){plot.subtitle <- paste0(ColID,' effects')}
                                   par(mar=c(7, 6, 6, 4) + 0.1)
                                   p1 <- ggplot(data=res.df,aes(y=sampquant,x=theoquant,label=factor(lab),colour=big)) + 
                                                geom_abline(aes(slope = slope, intercept = int), col='red') +
                                                geom_text(size=3) + theme(legend.position=c(0.85,0.25)) + 
                                                scale_colour_manual(values=c('blue','red')) +
                                                theme(plot.title = element_text(hjust = 0.5)) +
                                                labs(title=bquote(atop(.(plot.title), atop(italic(.(plot.subtitle)), ""))), 
                                                     x='Theoretical quantiles', y='Sample quantiles', colour='Absolute scres')
                                   options(warn=-1)
                                   print(p1)
                                   options(warn=0)
                                   
                                   # Scres heatmap plot
                                   res.df[[RowID]] <- factor(res.df[[RowID]],levels=rev(sort(unique(res.df[[RowID]]))))
                                   res.df <- res.df[order(res.df[[ColID]],res.df[[RowID]]),]
                                   res.df[!is.na(res.df[[scresID]]) & abs(res.df[[scresID]])>ceiling(tol),scresID] <- 
                                          sign(res.df[!is.na(res.df[[scresID]]) & abs(res.df[[scresID]])>ceiling(tol),scresID])*ceiling(tol)
                                   res.df <- as.data.frame(res.df)
                                   p2 <- ggplot(res.df, aes_string(ColID, RowID, fill=scresID, color="as.factor(big)")) + 
                                                geom_tile(size=1.5) + 
                                                scale_fill_gradient(low = "red", high = "green", limits=c(floor(-tol),ceiling(tol))) +
                                                scale_color_manual(values=c("transparent","yellow")) +
                                                theme(plot.title = element_text(hjust = 0.5)) +
                                                labs(title=plot.title, x=ColID, y=RowID, fill="scres", colour="abs scres")
                                   par(mar=c(6, 5, 5, 3) + 0.1)
                                   options(warn=-1)
                                   print(p2)
                                   options(warn=0)
                           }
                           dev.off()
                           if(i != length(yc)){message("Warning, GGplot Error, Please Re-run Residual Function")}
                           if(i == length(yc)){message("Output Figures Printed to File")}
                  } # end of resid plots function
########################################

##################################################################################################
# 
##################################################################################################
Colavsed.fn <- function(data, object, var.exclude=NULL, ...){
                        EnvID <- colnames(data)[grep("Env|YrLoc",colnames(data))]
                        var.exclude <- sort(unique(as.character(data$Variety[grep("fill|buff",tolower(data$Variety))])))
                        varlev.list <- with(droplevels(subset(data,!Variety %in% var.exclude & !is.na(yield))),tapply(Variety,Experiment,unique))
                        varlev.list <- lapply(varlev.list, function(x) sort(x))
                        varlev.list <- lapply(varlev.list, function(x) as.character(x))
                        enam <- levels(data[["Experiment"]])
                        ns <- length(enam)
                        ## loop
                        model.avsed <- list()
                        for(i in 1:ns){#i<- 1
                                message(paste0("Experiment ",i,"/",ns))
                                tmp <- data[data[["Experiment"]]==enam[i],]
                                tmp[["Variety"]] <- factor(tmp[["Variety"]])
                                pred.list <- list(enam[i],varlev.list[[enam[i]]])
                                names(pred.list) <- c("Experiment","Variety")
                                print(pred.list)
                                options(warn=-1)
                                temp.names <- rownames(object$coefficients$fixed)
                                temp.names <- temp.names[grep(":",temp.names, invert = T)]
                                if(length(temp.names[grep(EnvID,temp.names)])>0){
                                        pred <- predict(object, classify="Experiment:Variety",levels=pred.list, sed=T, maxit=1,
                                                        present= c("Experiment","Variety"),
                                                        average= c("Experiment","Variety"),
                                                        associate = ~Environment:Experiment)
                                }else{pred <- predict(object, classify="Experiment:Variety",levels=pred.list, sed=T, maxit=1)}
                                options(warn=0)
                                # means + sed
                                wh <- !is.na(pred$pvals$predicted.value)
                                model.avsed[[i]] <- pred$sed[wh,wh]
                       } # end of loop
                       if(i != ns){message("Warning, Function Error, Please Re-run Colavsed Function")}
                       if(i == ns){message("Avsed Data Saved to Object")}
                       names(model.avsed) <- enam
                       enam2 <- data[!duplicated(data$Experiment),c(EnvID,"Experiment")]
                       rownames(enam2) <- enam2$Experiment
                       enam2$Environment <- as.character(enam2$Environment)
                       enam2$Experiment <- as.character(enam2$Experiment)
                       all.df <- data.frame(Environment= enam2[enam,EnvID],
                                            Experiment = enam,
                                            Sum= sapply(model.avsed, function(x) sum(as.vector(x), na.rm=T)),
                                            Length= sapply(model.avsed, function(x) length(as.vector(x)[!is.na(as.vector(x))])))
                       expt.df <- data.frame(Experiment = factor(enam),
                                             avsed= round(all.df$Sum/all.df$Length,4))
                       rownames(expt.df) <- NULL
                       env.df <- data.frame(Environment = factor(unique(all.df$Environment)),
                                             avsed= round(with(all.df, tapply(Sum, Environment, sum))/with(all.df, tapply(Length, Environment, sum)),4))
                       rownames(env.df) <- NULL
                       colnames(env.df)[1] <- EnvID
                       if(length(grep("C|D|F", substring(data$Experiment,1,1))) > 0){
                         print(env.df);return(env.df)
                       }else{print(expt.df);return(expt.df)}
               }# end of avsed.fn
########################################

##################################################################################################
# This function builds the final model.csv based on avsed according to either 
# Experiment or Environment
Colfmodel.build <- function(avsed2, avsed3, models2, models3, tol=0.01){
                            name1 <- paste0("asr-",strsplit(deparse(substitute(avsed2)), split = "avsed")[[1]][2])
                            name2 <- paste0("asr-",strsplit(deparse(substitute(avsed3)), split = "avsed")[[1]][2])
                            EnvID <- names(models2)[grep("Env|YrLoc",colnames(models2))]
                            compare <- merge(avsed2,avsed3, by = 1)
                            colnames(compare)[2:3] <- c(name1,name2)
                            compare$Best <- compare[[name1]]
                            names(compare$Best) <- "Best"
                            compare$Fmodel <- name1
                            compare$Fmodel[compare[[name2]] < (compare[[name1]] - compare[[name1]]*tol)] <- name2
                            compare$Best[compare[[name2]] < (compare[[name1]] - compare[[name1]]*tol)] <- 
                                            compare[[name2]][compare[[name2]] < (compare[[name1]] - compare[[name1]]*tol)]
                            compare$"%ch" <-  round(((compare[[name1]] - compare$Best)/compare[[name1]]*100),2)
                            names(compare$"%ch") <- "%ch"
                            rownames(compare) <- NULL
                            print(compare)
                            compare2 <- compare
                            if(colnames(compare2)[1] == EnvID){compare2 <- merge(compare2,models2, by=EnvID)
                              }else{compare2 <- merge(compare2,models2, by= "Experiment")}
                            compare2$Experiment <- factor(compare2$Experiment, levels= unique(compare2$Experiment))
                            compare2 <- compare2[!duplicated(compare2),]
                            rownames(compare2) <- NULL
                            nnames <- colnames(compare2)[! colnames(compare2) %in% c(EnvID,"Experiment")]
                            nnames <- nnames[!nnames %in% c("asr-1","asr-2","asr-3","Best","Fmodel","%ch")]
                            Fmodels <- compare2[,c(EnvID,"Experiment",nnames,"Best","Fmodel","%ch")]
                            
                            # loop over avsed term
                            for(i in 1:nrow(Fmodels)){temp.model <- Fmodels$Fmodel[i] # i<- 1
                                      temp.Expt <- as.character(Fmodels$Experiment[i])
                                      if(temp.model == name2){temp.model2 <- models3
                                                          temp.model2$resid <- as.character(temp.model2$resid)
                                                          Fmodels[i,3:(ncol(Fmodels)-4)] <- temp.model2[temp.model2$Experiment == temp.Expt,3:(ncol(temp.model2)-1)]
                            }}
                            Fmodels <- Fmodels[order(Fmodels[[EnvID]],Fmodels$Experiment),]
                            rownames(Fmodels) <- NULL
                            print(Fmodels)
                            Fmodels.ls <- list()
                            Fmodels.ls$compare <- compare
                            Fmodels.ls$Fmodels <- Fmodels
                            return(Fmodels.ls)
                }# end of final models build function
########################################

##################################################################################################
# 
##################################################################################################
nvt <- function(outfile,data,model,ExptID='Experiment',VarietyID='Variety',large.no.exp=FALSE,non.standard.trial=FALSE, ...){
  # Author J Taylor with input from CL + BC
  # DT to do- I need to alter some things here
  #ASReml-R and NVT spreadsheet function
  #This function exports the results from a ASReml-R job of the analysis of NVT trials
  #to a format that it easily imported back into the NVT database
  # DT- updated 24 October 2016 for use in ASReml ver4.
  #DT + KM updated to fix errors
  message("Version October 2017")
  message("Running Predict on Experiments...")
  require(asreml)
  require(ASExtras4)
  require(XLConnect)
  trace <- "trace.txt"
  ftrace <- file(trace, "w")
  sink(trace, type = "output", append = FALSE)
  on.exit(sink(type = "output"))
  on.exit(close(ftrace), add = TRUE)
  ### start
  trait <- deparse(model$call$fixed[[2]])
  enam <- levels(data[[ExptID]])
  EnvID <- names(data)[grep("Env|YrLoc",colnames(data))]
  nv <- length(levels(data[[VarietyID]]))
  ns <- length(enam)
  pn <- paste(ExptID,VarietyID,sep=':')
  ExpVar <- c(ExptID,VarietyID)
  blues <- list(); sed <- c()
  for(i in 1:ns) {#i<-28
    evnam <- unique(trimws(data$Environment[data$Experiment==enam[i]]))
    message(paste0("Experiment ",i,"/",ns))
    tmp <- data[data[[ExptID]]==enam[i],]
    tmp <- tmp[tmp[[VarietyID]]!="Buffer",]
    tmp[[VarietyID]] <- factor(tmp[[VarietyID]])
    tt <- tapply(tmp[[trait]],tmp[[VarietyID]],mean,na.rm=T)
    pred.list <- list(enam[i],levels(tmp[[VarietyID]]))
    names(pred.list) <- ExpVar
    pred.list[[VarietyID]] <- pred.list[[VarietyID]][!pred.list[[VarietyID]] %in% names(tt[is.na(tt)])]
    options(warn=-1)
    if(ns > 1){temp.names <- rownames(model$coefficients$fixed)
               temp.names <- temp.names[grep(":",temp.names, invert = T)]
               if(length(temp.names[grep(EnvID,temp.names)])>0){
                    ignore.lst <- rownames(model$coefficients$fixed)[grep(paste0("at.*",evnam,".*:"),rownames(model$coefficients$fixed))]
                    if(length(ignore.lst)>0){cov <- unlist(lapply(strsplit(ignore.lst,":"),function(x) x[2]))
                          for(j in length(cov))if(length(unique(data[data$Experiment==enam[i],cov[j]]))>1|
                                                  length(unique(data[data$Experiment==enam[i],cov[j]]))==1 & !NA %in% unique(data[data$Experiment==enam[i],cov[j]]))ignore.lst[[j]] <- NA
                                ignore.lst <- ignore.lst[!is.na(ignore.lst)]
                          }
                    if(length(ignore.lst)==0)ignore.lst <- logical(0)
                    pvs <- predict(model, classify=pn,levels=pred.list, vcov=T, maxit=1, 
                                       present= c(ExptID,VarietyID),
                                       average= c(ExptID,VarietyID), ignore= ignore.lst,
                                       associate = ~Environment:Experiment)
               }else{pvs <- predict(model, classify=pn, levels=pred.list, vcov=T, maxit=1)}
    }else {pvs <- predict(model, classify=VarietyID, vcov=T,maxit=1)#, ...)
      pvs$pvals <- cbind.data.frame(rep(enam[i], nrow(pvs$pvals)), pvs$pvals)
      names(pvs$pvals)[1] <- ExptID
    }
    options(warn=0)
    pred <- pvs$pvals
    wh <- !is.na(pred$predicted.value)
    pred <- pred[wh,]
    pred$weights <- diag(solve(pvs$vcov[wh,wh]))
    tmp <- droplevels(tmp[tmp[[VarietyID]] %in% pred.list[[VarietyID]],])
    treps <- tapply(tmp[[trait]],tmp[[VarietyID]],function(x) length(x[!is.na(x)]))
    pred$truereps <- treps[names(treps) %in% pred[[VarietyID]]]
    pred$ems <- mean(pred$truereps/pred$weights, na.rm = T)
    blues[[i]] <- pred
    sed[i] <- pvs$avsed[2]
  }
  blues <- do.call("rbind", blues)
  print(str(blues))
  pred.mean <- tapply(blues$predicted.value, blues[[ExptID]],
                      mean, na.rm = TRUE)
  ems.mean <- tapply(blues$ems, blues[[ExptID]], mean, na.rm = TRUE)
  
  #################################################################################
  ##blues data.frame is in a column format, but NVT wants the data in row format
  ## i.e. different values stacked on top of each other.
  ##the rest of the script will put everything in the correct place.
  
  ##need 4 tables;
  #1. ETDI Trial - basic info on trials
  #2. ETDI Trial Cultivar Measurement - Analysed means, weights,... etc
  #3. ETDI Trial Measurement - analaysed trial summaries - CV, TMY, ... etc
  #4. ETDI Plot Measurement - plot data.. don't need to enter anything
  #################################################################################
  
  #############1. ETDI Trial#################################
  ###This contains the trial info...but will be blank for this purpose. but must be in file.
  #TrialCode*  Year*  Crop*	TrialTypeName*	Locality*
  
  nvt.table1 <- data.frame(TrialCode = NA,Year=NA,Crop=NA,TrialTypeName=NA, Locality=NA)
  message("Predict on Experiments Complete")
  message("Table 1 Complete")
  ##############2. ETDI Trial Cultivar Measurement########
  ##Need EMS, treps, blues, weight
  ##Separate values need to be stacked ontop of each other
  #TrialCode*  Cultivar*  Measurement*	MeasurementDate*	Value*
  
  tcm <- c("Predicted Yield","Weighting","Number of Successful Replicates in Trial",
           "Error Mean Squared","Standard Error")
  ltcm <- length(tcm)
  nvt.table2 <- data.frame(TrialCode = rep(blues[[ExptID]], ltcm), Cultivar = rep(blues[[VarietyID]],ltcm))
  nvt.table2$Measurement <- rep(tcm, each = nrow(blues))
  nvt.table2$MeasurementDate <- format(Sys.time(), "%d/%m/%Y")
  nrow(nvt.table2)
  head(blues)
  nvt.table2$Value <- c(blues$predicted.value, blues$weights, blues$truereps, blues$ems, blues$std.error)
  message("Table 2 Complete")
  ##################################################################################
  ########################3. ETDI Trial Measurement######################
  ###Need CV, varietal probability, SED, TMY, model line
  #TrialCode*  Measurement*	MeasurementDate*	Value*
  
  t3cm <- c("CV","Probability","Standard Error Difference","Yield Mean","Model Line")
  lt3cm <- length(t3cm)
  if(!large.no.exp & ns > 1) {
    wld <- wald.asreml(model, denDF='algebraic', ssType='conditional') 
    wldprob <- wld$Wald$Pr[grep(VarietyID,row.names(wld$Wald))]
    dendf <- wld$Wald$denDF[grep(VarietyID,row.names(wld$Wald))] ##doesn't always work??
    if(any(is.na(dendf))){
      wld <- wald.asreml(model)
      wldprob.i <- wld[grep(VarietyID,row.names(wld)),"Pr(Chisq)"]
      wldprob[is.na(dendf)] <- wldprob.i[is.na(dendf)]
    }
  }else {
    wld <- wald.asreml(model)
    wldprob <- wld[grep(VarietyID,row.names(wld)),'Pr(Chisq)']
  }
  wldprob[wldprob < 0.001] <- '<0.001'
  if(!non.standard.trial){
    formn <- c("fixed","random","sparse","residual")
    forml <- list()
    for(i in 1:length(formn)){# i<-4
      form <- model$call[[formn[i]]]
      wform <- formn[i]
      fixed <- FALSE
      if(!is.null(form)){
        if(wform == "fixed"){
          fixed <- TRUE
          resp <- deparse(form[[2]])
          form <- deparse(form[[3]])
        } else form <- deparse(form[[2]])
        str <- strsplit(paste(form, collapse = ""), "\\+")[[1]]
        str <- gsub(" ","", str)
        imat <- NULL
        strm <- str
        if(ns != 1){
          wh <- grep(":", str)
          if(length(wh)){
            strm <- strm[-wh]
            strm <- strm[!(strm %in% ExptID)]
            stri <- str[wh]
            ilevs <- gsub(".*?t,(.*?)\\):.*", "\\1", stri)
            wh <- sapply(gregexpr("\\)", stri), function(x) x[1])
            iv <- substring(stri, wh + 2, nchar(stri))
            ch <- sapply(gregexpr(":", iv), "[", 1) %in% 1
            iv[ch] <- substring(iv[ch], 2, nchar(iv[ch]))
            imat <- matrix(nrow = length(enam), ncol = length(stri))
            for(j in 1:length(stri)){#j<-1
              if(!length(grep(ExptID, ilevs[j]))){
                ilevs[j] <- gsub("c\\(|\\)","", ilevs[j])
                tlev <- strsplit(ilevs[j], "\\,")[[1]]
                if(length(grep("\"",tlev))){
                  wh <- pmatch(gsub("\"","",tlev), enam)
                  imat[wh,j] <- iv[j]
                } else {
                  if(length(grep(":",tlev)))
                    nums <- unlist(lapply(strsplit(tlev, ":"), function(x){
                      x <- as.numeric(factor(x))
                      if(length(x) > 1)
                        x[1]:x[2]
                      else x  }))
                  else nums <- as.numeric(tlev)
                  imat[nums,j] <- iv[j]
                }
              } else  imat[,j] <- iv[j]
            }
          }
        }
        if(length(strm))
          for(j in length(strm):1) imat <- cbind(strm[j],imat)
        forml[[wform]] <- apply(imat, 1, function(x){
          x <- x[!is.na(x)]
          if(!length(x)) NULL
          else paste("~",paste(x, collapse = "+"), collapse = "")
        })
        if(fixed) forml[[wform]] <- paste(resp, forml[[wform]], sep = " ")
        forml[[wform]] <- paste(paste(wform, "=", collapse = ""), forml[[wform]], sep = " ")
      }
    }
    forms <- do.call("cbind", forml)
    forms <- apply(forms, 1, function(x) {
      x <- x[grep("~",x)]
      paste(x, collapse = " , ")
    })
    res <- residuals(model)
    res[is.na(res)] <- 0
    nvt.table4 <- data.frame(TrialCode=data[[ExptID]],Row=as.numeric(as.character(data$Row)),
                             Range=as.numeric(as.character(data$Range)),
                             Measurement='Residual',MeasurementDate=format(Sys.time(), "%d/%m/%Y"),
                             Value=res)
  }else {
    model.line <- "Yield ~ VarietyID, random=~ Rep, residual=~ ar1(Range):ar1(Row)"
    forms <- rep(model.line, ns)
    nvt.table4 <- data.frame(TrialCode=NA,Row=NA,Range=NA,
                             Measurement=NA,MeasurementDate=NA,Value=NA)
  }
  nvt.table3 <- data.frame(TrialCode = rep(enam,lt3cm), Measurement= rep(t3cm,each = ns))
  nvt.table3$MeasurementDate <- format(Sys.time(), "%d/%m/%Y")
  nvt.table3$Value <- c(sqrt(ems.mean)/pred.mean*100, wldprob, sed, pred.mean, forms)
  nvt.table3 <- nvt.table3[order(nvt.table3$TrialCode,rep(1:lt3cm, each = ns)),]
  
  message("Table 3 Complete")
  ########################4. ETDI Plot Measurement######################
  ###Don't need to enter anything in here!
  #TrialCode*  Row*	Range*	Measurement*	MeasurementDate*	Value*
  
  message("Table 4 Complete")
  
  ##############export all tables to xls file, with different tables in different spreadsheets
  
  filenam <- paste(outfile,'xls',sep='.') #KLM: changes this from enam (which is already in use to filenam)
  analysedmeans <- loadWorkbook(filenam, create = TRUE) #KLM: changes this from enam (which is already in use to filenam)
  
  head <- rep(NA, dim(nvt.table1)[2])
  head[c(1,length(head))] <- c("HeaderStart","HeaderEnd")
  createSheet(analysedmeans, name = "ETDI Trial")
  writeWorksheet(analysedmeans, t(head), sheet = "ETDI Trial",rownames=NULL,header=FALSE)
  setColumnWidth(analysedmeans, sheet = "ETDI Trial",1:5,width=-1)
  
  head <- rep(NA, dim(nvt.table4)[2])
  head[c(1,length(head))] <- c("HeaderStart","HeaderEnd")
  createSheet(analysedmeans, name = "ETDI Plot Measurement")
  writeWorksheet(analysedmeans, t(head), sheet = "ETDI Plot Measurement",rownames=NULL,header=FALSE)
  writeWorksheet(analysedmeans, nvt.table4, sheet = "ETDI Plot Measurement",startRow=2,rownames=NULL)
  setColumnWidth(analysedmeans, sheet = "ETDI Plot Measurement",1:6,width=-1)
  
  # nvt.table2[nvt.table2$Measurement=="Standard Error" & nvt.table2$TrialCode== "FMaA16WALG2",]
  head <- rep(NA, dim(nvt.table2)[2])
  head[c(1,length(head))] <- c("HeaderStart","HeaderEnd")
  createSheet(analysedmeans, name = "ETDI Trial Cultivar Measurement")
  writeWorksheet(analysedmeans, t(head), sheet = "ETDI Trial Cultivar Measurement",rownames=NULL,header=FALSE)
  writeWorksheet(analysedmeans, nvt.table2, sheet = "ETDI Trial Cultivar Measurement",startRow=2,rownames=NULL)
  setColumnWidth(analysedmeans, sheet = "ETDI Trial Cultivar Measurement",1:5,width=-1)
  
  head <- rep(NA, dim(nvt.table3)[2])
  head[c(1,length(head))] <- c("HeaderStart","HeaderEnd")
  createSheet(analysedmeans, name="ETDI Trial Measurement")
  writeWorksheet(analysedmeans, t(head), sheet = "ETDI Trial Measurement",rownames=NULL,header=FALSE)
  writeWorksheet(analysedmeans, nvt.table3, sheet = "ETDI Trial Measurement",startRow=2,rownames=NULL)
  setColumnWidth(analysedmeans, sheet = "ETDI Trial Measurement",1:4,width=-1)
  
  saveWorkbook(analysedmeans)
  
  message("Data Results Exported to File")
} #END of NVT function
########################################

##################################################################################################
# 
##################################################################################################
results.print <- function(data, outfile){
                          require(ggplot2)
                          data$plot.grp <- ceiling(as.numeric(factor(data$Experiment))/9)
                          pdf(outfile)
                          np <- length(unique(data$plot.grp))
                          for(i in 1:np){#i <- 1
                                    p.data <- droplevels(subset(data, plot.grp==i))
                                    p1 <- ggplot(data=p.data,mapping=aes(x=rawYield,y=predictedYield)) + geom_point() + 
                                                 geom_abline(col='blue') + facet_wrap(~Experiment, scales='free') +
                                                 labs(x="Raw Yield (t/ha)", y="Predicted Yield (t/ha)")
                                    options(warn=-1)
                                    print(p1)
                                    options(warn=0)}
                          dev.off()
                          if(i == np){message("Output Figures Printed to File")}
                          if(i != np){message("Warning, GGplot Error, Please Re-run Results Function")}
                  } #end of results print function
########################################

