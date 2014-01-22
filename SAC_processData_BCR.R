# Process BCR data:
# - for each BCR, find the sampling_IDs from core-covariates.csv file
# - extract the checklist for each BCR from checklists.csv file
# - fit mixture of SACs 
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
#args <- commandArgs(trailingOnly = TRUE)
BCRs <- c("30","31","32","37")
years <- c("2011")

# convert ERD to bird binary data
for (year in years) {    
    # load checklists
    dataFile <- paste("../data/SoCS/",year,"/all_checklists.RData",sep="")
    load(dataFile)
    
    # load covs
    covFile <- paste("../data/SoCS/",year,"/core-covariates.csv",sep="")
    cov <- read.csv(covFile,head=TRUE,stringsAsFactors=FALSE)
    
    for (BCR in BCRs) {
        cat(BCR,year,"\n")
        
        samplingIDs <- cov$SAMPLING_EVENT_ID[which(cov$BCR == BCR)]
        BCRData <- data[which(data$SAMPLING_EVENT_ID %in% samplingIDs),]
        
        # save data
        outputFile <- paste("../data/SoCS/",year,"/BCR",BCR,"_checklists.RData",sep="")
        save(BCR,year,BCRData,file=outputFile)
    }
} 
