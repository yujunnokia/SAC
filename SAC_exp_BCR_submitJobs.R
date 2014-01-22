# Submit jobs to cluster to run SAC exp on BCRs
# 
# Author: Jun Yu
# Version: Oct 2013
###############################################################################

states <- c("BCR30","BCR31","BCR32","BCR37")
Ks     <- 1:6

for (state in states) 
{
    for (K in Ks) 
    {
        script <- paste("#!/bin/csh\n\n", 
                        "#$ -N ",state,".",K,
                        "\n\n# set working directory on all host to",
                        "\n# directory where the job was started",
                        "\n#$ -cwd",
                        "\n",
                        "\n# send all process STDOUT (fd 2) to this file",
                        "\n#$ -o ./tmp/",state,"_",K,"_o.txt",
                        "\n",
                        "\n# send all process STDERR (fd) to this file",
                        "\n#$ -e ./tmp/",state,"_",K,"_e.txt",
                        "\n", 
						"\n#$ -q eecs",
						"\n",
                        "\n# Commands \n",
                        "Rscript ./SAC_exp_BCR.R ",state," ",K," \n", sep="")
        
        file <- paste("./scripts/",state,"_",K,".sh", sep="")
        write(script, file=file, append = FALSE)
        system(paste("qsub",file))
        Sys.sleep(1.0)
    }
}
