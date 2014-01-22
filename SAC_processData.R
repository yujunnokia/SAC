# 
# 
#0 SAMPLING_EVENT_ID,
#1 LATITUDE,
#2 LONGITUDE,
#3 YEAR,
#4 MONTH,
#5 DAY,
#6 TIME,
#7 COUNTRY,
#8 STATE_PROVINCE,
#9 COUNT_TYPE,
#10 EFFORT_HRS,
#11 EFFORT_DISTANCE_KM,
#12 EFFORT_AREA_HA,
#13 OBSERVER_ID,
#14 NUMBER_OBSERVERS,
#15 GROUP_ID,
#16 PRIMARY_CHECKLIST_FLAG,
# Author: Jun Yu
###############################################################################

rm(list=ls())

setwd("/nfs/guille/tgd/wonglab/yuju/SAC")
        
# load ebird data
args <- commandArgs(trailingOnly = TRUE)
state <- args[1]
years <- c("2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012")

# convert ERD to bird binary data
for (year in years) {
    cat(state,year,"\n")
    
    file <- paste("../data/SoCS/",year,"/",state,"_checklists.csv",sep="")
    data <- read.csv(file,head=TRUE,stringsAsFactors=FALSE)
    
    # aggregate count
    observations <- data[,18:ncol(data)]
    observations[observations > 1] <- 1
    data <- cbind(data[,1:17],observations)
    
    # save data
    outputFile <- paste("../data/SoCS/",year,"/",state,"_allchecklists.RData",sep="")
    save(state,year,data,file=outputFile)
} 
