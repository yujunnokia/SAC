# Join checklists with site covariates
# 
#0 SAMPLING_EVENT_ID
#1 LATITUDE
#2 LONGITUDE
#3 YEAR
#4 MONTH
#5 DAY
#6 TIME
#7 COUNTRY
#8 STATE_PROVINCE
#9 COUNT_TYPE
#10 EFFORT_HRS
#11 EFFORT_DISTANCE_KM
#12 EFFORT_AREA_HA
#13 OBSERVER_ID
#14 NUMBER_OBSERVERS
#15 GROUP_ID
#16 PRIMARY_CHECKLIST_FLAG
#
# Author: Jun Yu
# Version: Oct 2013
###############################################################################

rm(list=ls())

setwd("/nfs/guille/tgd/wonglab/yuju/SAC")
        
# load ebird data
args <- commandArgs(trailingOnly = TRUE)
state <- args[1]
years <- c("2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")

for (year in years) {
    cat(state,year,"\n")
    
    # load checklist file
    dataFile <- paste("../data/SoCS/",year,"/",state,"_checklists.csv",sep="")
    data <- read.csv(dataFile,head=TRUE,stringsAsFactors=FALSE)
    nCol <- ncol(data)

    # load covariate file
    covFile <- paste("../data/SoCS/",year,"/core-covariates.csv",sep="")
    cov <- read.csv(covFile,head=TRUE,stringsAsFactors=FALSE)
    
    # merge
    data <- merge(data,cov,by="SAMPLING_EVENT_ID")
    
    # aggregate count
    observations <- data[,18:nCol]
    observations[observations > 1] <- 1
    
    # put species to the end of a checklist
    #counts <- rowSums(observations)
    #data <- cbind(data[,1:17],counts,data[,(nCol+1):ncol(data)],observations)
    data <- cbind(data[,1:17],data[,(nCol+1):ncol(data)],observations)
    
    outputFile <- paste("../data/SoCS/",year,"/",state,"_cov_checklists.csv",sep="")
    write.csv(data,outputFile,row.names=FALSE,quote=FALSE)
    
    # save data
    #outputFile <- paste("../data/SoCS/",year,"/",state,"_cov_checklists.RData",sep="")
    #save(state,year,data,file=outputFile)
} 
