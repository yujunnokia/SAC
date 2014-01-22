# report the most differential species between two groups
# 
# Author: Jun Yu
# Version: Oct 2013
###############################################################################

rm(list=ls())

library(ggplot2)

# set working directory
#setwd("/Users/yujunnokia/workspace/SAC")
setwd("/nfs/guille/tgd/wonglab/yuju/SAC")

source("./SAC_util.R")

args <- commandArgs(trailingOnly = TRUE)
state <- args[1]
diff <- args[2] # most least
g1 <- as.integer(args[3])
g2 <- as.integer(args[4])
year <- "2012"
reportSpecies <- 20
TOP_SPECIES <- 200

#if (!(state %in% c("New_York","Florida","Texas","California"))) { stop("state is invalid...\n") }
if (!(diff %in% c("most","least"))) { stop("diff is invalid...\n") }

cat("Report in",state,"...\n")

# load bird names
birds <- read.csv("../data/SoCS/bird/taxonomy.csv",stringsAsFactors=FALSE) 
uniqueSpecies <- unique(birds$SPECIESNAME)
uniqueSpecies <- uniqueSpecies[uniqueSpecies != ""]

# load birder data
birders <- read.csv(paste("../data/SoCS/birder/",state,"_cluster.csv",sep=""),stringsAsFactors=FALSE)
K <- length(unique(birders$groupID))
nBirders <- nrow(birders)

# load checklist data
dataFile <- paste("../data/SoCS/",year,"/",state,"_allchecklists.RData",sep="")
load(dataFile)
data <- BCRData

# merge subspecies into species
observations <- matrix(0,nrow=nrow(data),ncol=length(uniqueSpecies))
commonNames <- array(0,length(uniqueSpecies))
for (s in 1:length(uniqueSpecies)) {
    species <- uniqueSpecies[s]
    indices <- birds$IDX[which(birds$SPECIESNAME == species)]
    if (length(indices) == 1) {
        observations[,s] <- as.integer(data[,indices])
        commonNames[s]   <- birds$COMMONNAME[which(birds$SPECIESNAME == species)]
    } else {
        observations[,s] <- rowSums(data.matrix(data[,indices]))
        commonNames[s]   <- birds$COMMONNAME[which(birds$SPECIESNAME == species)[1]]
    }
}
observations[observations > 0] <- 1
colnames(observations) <- commonNames

# keep the top 200 most frequently reported species
speciesIdx <- order(colMeans(observations),decreasing=TRUE)[1:TOP_SPECIES]
frequentBirds <- commonNames[speciesIdx]
data <- cbind(data[,1:17],observations[,speciesIdx])
print(frequentBirds[1:10])

# calculate detection ratio
df <- data.frame(species=array(0,nBirders*TOP_SPECIES),
                 group  =array(0,nBirders*TOP_SPECIES),
                 metric =array(0,nBirders*TOP_SPECIES))
means <- matrix(0,nrow=K,ncol=TOP_SPECIES)
ses <- matrix(0,nrow=K,ncol=TOP_SPECIES)
count <- 1
for (k in 1:K) {
    group <- birders$obsID[birders$groupID == k]
    
    totalCLs <- 0
    groupStats <- matrix(0,nrow=length(group),ncol=TOP_SPECIES)
    for (i in 1:length(group)) {
        member <- group[i]
        memberData <- data[which(data$OBSERVER_ID == member),]
        totalCLs <- totalCLs + nrow(memberData)
        observations <- memberData[,18:ncol(memberData)]
        
        groupStats[i,] <- colMeans(observations)
        
        for (s in 1:TOP_SPECIES) {
            df$species[count] <- frequentBirds[s]
            df$group[count] <- k
            df$metric[count] <- mean(observations[,s])
            count <- count + 1
        }
    }
    
    means[k,] <- colMeans(groupStats)
    ses[k,] <- apply(groupStats, 2, sd) / sqrt(length(group))
    
    cat("State:",state,"Group:",k,"Birders:",length(group),"Total CLs:",totalCLs,"\n")
}
colnames(means) <- frequentBirds
colnames(ses) <- frequentBirds
output <- rbind(means,ses)

if (diff == "most") {
    topSpecies <- frequentBirds[order(means[g1,] / means[g2,], decreasing=TRUE)[1:reportSpecies]]
} else {
    topSpecies <- frequentBirds[order(means[g1,] / means[g2,], decreasing=FALSE)[1:reportSpecies]]
}
print(topSpecies)

df <- df[which(df$species %in% topSpecies),]
df$group <- factor(df$group,levels=1:K)
df$species <- factor(df$species,levels=topSpecies)
dfc <- summarySE(df, measurevar="metric", groupvars=c("group","species"))

# plot
p <- ggplot(dfc, aes(x=species, y=metric, fill=group)) + 
    scale_fill_brewer(palette="Set1") +
    ggtitle(paste(gsub("_", " ", state),g1,"vs.",g2)) +
    ylab("Detection Rate") +
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=metric-se, ymax=metric+se),width=.2,position=position_dodge(.9)) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
    theme(plot.title=element_text(family="Times", face="bold", size=20),legend.position="right") +
    theme(legend.title=element_text(face="bold",size=15), legend.text=element_text(face="bold",size=15)) +
    theme(axis.title.x=element_text(face="bold",size=15), axis.text.x=element_text(face="bold",size=15)) +
    theme(axis.title.y=element_text(face="bold",size=15), axis.text.y=element_text(face="bold",size=15))
    
file <- paste("./plots/",state,"_",year,"_",g1,"vs",g2,"_",diff,"_Stats.jpeg",sep="")
ggsave(filename=file, plot=p, width=15,height=12)


outputFile <- paste("./results/",state,"_",g1,"vs",g2,"_",diff,"_stats.csv",sep="")
write.csv(output,outputFile,quote=FALSE,row.names=FALSE)





