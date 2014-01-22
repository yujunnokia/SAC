#! /usr/bin/env Rscript

states <- c("New_York","Florida","Texas","California")
diffs   <- c("least","most")
g1s    <- c(1,1,2)
g2s    <- c(2,3,3)

for (state in states) 
{
    for (diff in diffs) 
    {
        for (i in 1:length(g1s)) 
        {
			g1 <- g1s[i]
			g2 <- g2s[i]
			
            system(paste("qdel ",substr(state,1,1),".",substr(diff,1,1),".",g1,".",g2,sep=""))
        }
    }
}
