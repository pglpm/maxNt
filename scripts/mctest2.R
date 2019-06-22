## this script computes the probability for the conditional frequencies
## frequency(symptom | snp)
## where each symptom can assume values listed in 'symptomvariants'
## and each snp can assume values (alleles) listed in 'snpvariants'
## such values can be tuples
library('ggplot2')
library('RColorBrewer')
library('cowplot')
library('png')
library('plot3D')
library('doParallel')
library('LaplacesDemon')
library('dplyr')

filename <- 'mcsamples.csv'

cores <- 1
thinning <- 1
mciterations <- 10000
mcdiscard <- 10

## setup parallel processing
`%dox%` <- `%do%`
if(cores>1){
`%dox%` <- `%dopar%`
cl <- makeCluster(cores, outfile="")
registerDoParallel(cl)
}

### Here is the set-up for the Markov-chain Monte Carlo sampling of the parameters
### Uses package LaplacesDemon

        PGF <- function(data){1}

        mydata <- list(y=1, PGF=PGF,
               parm.names=c('y'),
               mon.names=c('')
               )

        logprob <- function(parm,data){
            parm <- interval(parm,0,1)
            LP <- log(2*parm)
            list(LP=LP, Dev=-2*LP, Monitor=1, yhat=1, parm=parm)
}

        ##Initial.Values <- c(0,0) #log(mf/sum(mf)) #GIV(prob, mydata, n=1000, PGF=T)
        
        Initial.Values <- GIV(logprob, mydata, n=1000, PGF=F)

        usamples <- LaplacesDemon(logprob, mydata, Initial.Values,
                        Covar=NULL,
                        Thinning=thinning,
                        Iterations=mciterations, Status=mciterations,
                        #Chains=nchains,CPUs=nchains,LogFile=paste0(filename,'_LDlog'), #Packages=c('Matrix'),#Type="MPI",
                        ##Algorithm="RDMH"#, Specs=list(B=list(1:d,d1:d2,d3:dnp))
                        ##Algorithm="Slice", Specs=list(B=NULL, Bounds=c(0,1), m=Inf, Type="Continuous", w=0.001)
                        Algorithm="AFSS", Specs=list(A=mcdiscard, B=NULL, m=100, n=0, w=1)
                        )
    samples2 <- usamples$Posterior2
    samples1 <- usamples$Posterior1
    
    print('writing samples to file...')
#    write.csv(samples1,'testmcsamples.csv',row.names=F,col.names=F)
    write.table(samples1,file='testmcsamples2.csv',row.names=F,col.names=F, sep=',')
