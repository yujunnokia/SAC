# SAC experiment
# 
# Author: Jun Yu
# Version: Oct 2013
###############################################################################

rm(list=ls())

# set working directory
#setwd("/Users/yujunnokia/workspace/SAC")
setwd("/nfs/guille/tgd/wonglab/yuju/SAC")
library(ggplot2)

source("SAC.R")

# turn off warning messages
options(warn=-1)

EFFORT_DURATION <- 120
CL_THRESHOLD    <- 20

states <- c("New_York") #c( "New_York", "Florida", "Texas", "California")
year <- "2012"  

Ks <- c(3) #c(1:6)
sigma <- NULL 
lambda <- 0.1
nRandomRestarts <- 10

for (state in states) {
	# load birder data
	birderData <- read.csv("../data/SoCS/birder/birders_NY.csv",stringsAsFactors=FALSE)
	
	# load eBird data
	dataFile <- paste("../data/SoCS/",year,"/",state,"_checklists.RData",sep="")
	load(dataFile)
	
	# remove checklists whose effort time is 0
	data <- data[data$EFFORT_HRS != "?",]
	data$EFFORT_HRS <- as.numeric(data$EFFORT_HRS)
	data <- data[data$EFFORT_HRS > 0,]
	# observation time window
	data$EFFORT_HRS <- data$EFFORT_HRS*60
	data <- data[data$EFFORT_HRS <= EFFORT_DURATION,]
	
	X <- matrix(sqrt(data$EFFORT_HRS),ncol=1)
	Y <- data$counts

	B <- list()
	uniqueBirders <- unique(data$OBSERVER_ID)
	birderIDs <- NULL
	
	count <- 1
	for (i in 1:length(uniqueBirders)) {
		birder <- uniqueBirders[i]
		idx <- which(data$OBSERVER_ID == birder)
		if (length(idx) > CL_THRESHOLD) {
			birderIDs <- c(birderIDs,birder)
			B[[count]] <- idx
			count <- count + 1
		}
	}
	cat("Number of birders:",length(B),"\n")
	
	# learn MSAC
	for (K in Ks) {
		# train SACs
		model <- RandomRestartEM.SAC(K,X,Y,B,lambda,nRandomRestarts,sigma) 
		
		# predict
		predictions <- NULL
		for (k in 1:K) {
			prediction <- data.frame(EFFORT_HRS=seq(0, EFFORT_DURATION, length.out=100),
					X=sqrt(seq(0, EFFORT_DURATION, length.out=100)))
			prediction <- cbind(prediction,predict(model$models[[k]], prediction, type = "link", se.fit=TRUE))
			prediction$group <- paste("G",k,"(",round(model$pi[k], digits=2)*100,"%)",sep="")
			
			predictions <- rbind(predictions, prediction)
		}		
		predictions <- within(predictions, {group <- as.factor(group)
					LL <- fit - 1.96 * se.fit
					UL <- fit + 1.96 * se.fit
				})
		
		# plot the SACs
		p <- ggplot(predictions, aes(x=EFFORT_HRS, y=fit)) +	
				scale_colour_brewer(palette="Set1") +
				coord_cartesian(xlim=c(-1,EFFORT_DURATION*1.1), ylim=c(0, 50)) +
				ggtitle(gsub("_", " ", state)) +
				xlab("Duration in Mins") + 
				ylab("Number of species on a checklist") +
				theme(plot.title=element_text(family="Times", face="bold", size=20),legend.position="right") +
				theme(legend.title=element_text(face="bold",size=15), legend.text=element_text(face="bold",size=15)) +
				theme(axis.title.x=element_text(face="bold",size=15), axis.text.x=element_text(face="bold",size=15)) +
				theme(axis.title.y=element_text(face="bold",size=15), axis.text.y=element_text(face="bold",size=15)) +
				geom_ribbon(aes(ymin = LL, ymax = UL, fill = group), alpha = .3) +
				geom_line(aes(colour=group), size = 1) 
		print(p)
		#file <- paste("plots/",state,"_",year,"_",K,"_SAC.jpeg",sep="")
		#ggsave(filename=file, plot=p, width=5,height=4)
		
		# assign each birder to a group
		partition <- data.frame(obsID=birderIDs,name=array(0,length(birderIDs)),numCLs=array(0,length(birderIDs)),
				groupID=array(0,length(birderIDs)),expertise=array(0,length(birderIDs)))
		for (i in 1:length(birderIDs)) {
			birderID <- birderIDs[i]
            partition$numCLs[i] <- length(B[[i]])
			partition$groupID[i] <- order(model$Z[i,],decreasing=TRUE)[1]
			idx <- which(birderData$OBSERVER_ID == birderID)
			if (length(idx) == 0) {
				partition$name[i] <- "Not Found"
				partition$expertise[i] <- "Not Known"
			} else {
				partition$name[i] <- birderData$OBSERVER[idx]
				if (birderData$expertise[idx] == "e" || birderData$expertise[idx] == "n") {
					partition$expertise[i] <- birderData$expertise[idx]
				} else {
					partition$expertise[i] <- "Not Known"
				}
			}
		}
		
		# print first class
		partition <- partition[order(partition[,4]),]
		print(partition[partition[,4]==1,])
		#write.csv(partition,file=paste("./results/",state,"_cluster.csv",sep=""),row.names=FALSE,quote=FALSE)
		
		# save result
		#resultFile <- paste("./results/",state,"_",year,"_",K,"_SAC.RData",sep="")
		#save(state,year,K,data,X,Y,B,sigma,model,partition,file=resultFile)
	}
}


