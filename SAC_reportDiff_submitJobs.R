#! /usr/bin/env Rscript

states <- c("New_York","Florida","Texas","California") # c("BCR31")
diffs   <- c("most")
g1s    <- c(1)
g2s    <- c(3)

for (state in states) 
{
    stateAbb <- substr(state,1,1)
    for (diff in diffs) 
    {
        for (i in 1:length(g1s)) 
        {
            g1 <- g1s[i]
            g2 <- g2s[i]
            
            script <- paste("#!/bin/csh\n\n", 
                                "#$ -N ",stateAbb,".",substr(diff,1,1),".",g1,".",g2,
                                "\n\n# set working directory on all host to",
                                "\n# directory where the job was started",
                                "\n#$ -cwd",
                                "\n",
                                "\n# send all process STDOUT (fd 2) to this file",
                                "\n#$ -o ./tmp/",stateAbb,"_",substr(diff,1,1),"_",g1,"_",g2,"_o.txt",
                                "\n",
                                "\n# send all process STDERR (fd) to this file",
                                "\n#$ -e ./tmp/",stateAbb,"_",substr(diff,1,1),"_",g1,"_",g2,"_e.txt",
                                "\n", 
                                "\n#$ -q eecs",
                                "\n",
                                "\n# Commands \n",
                                "Rscript ./SAC_reportDiff.R ",state," ",diff," ",g1," ",g2,"\n", sep="")
                                #"Rscript ./SAC_reportDiff_BCR.R ",state," ",diff," ",g1," ",g2,"\n", sep="")
            
            file <- paste("./scripts/",state,"_",diff,"_",g1,"_",g2,".sh", sep="")
            write(script, file=file, append = FALSE)
            system(paste("qsub",file))
            Sys.sleep(1.0)
        }
    }
}
